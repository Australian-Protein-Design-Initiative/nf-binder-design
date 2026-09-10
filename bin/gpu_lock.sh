#!/usr/bin/env bash
# Atomic, crash-safe GPU slot allocation for single-node multi-GPU runs.
#
# Sourced from a process script; sets CUDA_VISIBLE_DEVICES to a GPU that no
# other pipeline task holds, and keeps holding that claim until the task exits.
#
# Why not read GPU usage from nvidia-smi
# -------------------------------------
# Because it cannot see the tasks. Apptainer is invoked with `--pid`, so each
# container gets its own PID namespace and `nvidia-smi --query-compute-apps`
# inside it reports "No running processes found" however loaded the card is.
# Falling back to `memory.used` does not work either: a task needs ~25 s of
# Python imports and checkpoint loading before it allocates any VRAM, so every
# task starting inside that window sees the same idle card and picks it.
#
# A static round-robin (CUDA_VISIBLE_DEVICES = task.index % nGPU) is not a fix
# either: it assigns by submission order, which drifts out of phase with
# completion order as soon as task durations vary. Measured over a 2x RTX 3060
# run, 41.4% of wall time had two tasks stacked on one card while the other sat
# completely idle, despite each card receiving the same number of tasks.
#
# The mechanism here
# ------------------
# One lock file per GPU slot, claimed with `flock -n`. flock is atomic, so two
# tasks cannot win the same slot no matter how closely they race, and the claim
# is a property of an open file descriptor: the kernel drops it when the task
# exits, including on SIGKILL or a machine crash. There are no stale locks to
# reap and no PID bookkeeping to get wrong.
#
# The descriptor is deliberately leaked into the caller's shell. That is what
# holds the slot for the lifetime of the task.
#
# Usage:
#   source /path/to/gpu_lock.sh
#   nfbd_acquire_gpu "0,1" "/path/to/lockdir" 1 3600

# Resolve a candidate spec ("0,1" / "all" / "auto") to a list of GPU indices.
nfbd_gpu_candidates() {
    local spec="${1:-all}"
    if [[ -z "${spec}" || "${spec}" == "all" || "${spec}" == "auto" ]]; then
        nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null | tr -d ' '
    else
        tr ',' '\n' <<< "${spec}" | sed '/^$/d' | tr -d ' '
    fi
}

# nfbd_acquire_gpu <candidates> <lock_dir> [slots_per_gpu] [timeout_seconds]
#
# On success: exports CUDA_VISIBLE_DEVICES and returns 0, holding the slot.
# On timeout: returns 1 without setting CUDA_VISIBLE_DEVICES.
nfbd_acquire_gpu() {
    local spec="${1:?candidates required}"
    local lock_dir="${2:?lock dir required}"
    local slots="${3:-1}"
    local timeout="${4:-3600}"

    if ! command -v flock >/dev/null 2>&1; then
        echo "gpu_lock: flock not found in this container; cannot allocate a GPU safely." >&2
        return 2
    fi

    local -a gpus
    mapfile -t gpus < <(nfbd_gpu_candidates "${spec}")
    if [[ ${#gpus[@]} -eq 0 ]]; then
        echo "gpu_lock: no GPUs matched '${spec}'" >&2
        return 2
    fi

    mkdir -p "${lock_dir}" || { echo "gpu_lock: cannot create ${lock_dir}" >&2; return 2; }

    # Build the slot list SLOT-MAJOR: slot 0 of every GPU, then slot 1 of every
    # GPU, and so on. Scanning it in order therefore fills every card's first
    # slot before doubling up on any of them, which is what makes placement
    # least-loaded-first.
    #
    # Building it GPU-major instead (all of GPU 0's slots, then all of GPU 1's)
    # looks equivalent and is not: with slots > 1 a task takes the next free
    # slot on the SAME card in preference to an idle one, which reproduces
    # exactly the stacked-on-one-card-while-the-other-idles behaviour this
    # script exists to prevent. Measured at 23.7% of wall time before the
    # ordering was fixed.
    #
    # The GPU order is rotated per task within each slot tier. Rotating the
    # whole list instead would break the tiering, since a rotation could put a
    # slot-1 entry ahead of an unclaimed slot-0.
    local -a slot_gpu slot_file
    local g s i
    local ngpu=${#gpus[@]}
    local offset=$(( ${RANDOM:-$$} % ngpu ))
    for (( s = 0; s < slots; s++ )); do
        for (( i = 0; i < ngpu; i++ )); do
            g=${gpus[$(( (i + offset) % ngpu ))]}
            slot_gpu+=("${g}")
            slot_file+=("${lock_dir}/gpu${g}.slot${s}.lock")
        done
    done
    local n=${#slot_gpu[@]}

    local deadline=$(( SECONDS + timeout ))
    local announced=0
    while :; do
        local idx fd
        for (( idx = 0; idx < n; idx++ )); do
            # Open first, then try the lock. If the lock fails the descriptor is
            # closed again, otherwise it stays open and the claim persists.
            exec {fd}>"${slot_file[$idx]}" || continue
            if flock -n "${fd}"; then
                export CUDA_VISIBLE_DEVICES="${slot_gpu[$idx]}"
                echo "gpu_lock: acquired GPU ${slot_gpu[$idx]} via $(basename "${slot_file[$idx]}")"
                return 0
            fi
            exec {fd}>&-
        done

        if (( SECONDS >= deadline )); then
            echo "gpu_lock: no GPU slot free after ${timeout}s (candidates: ${gpus[*]}, ${slots} slot(s) each)" >&2
            return 1
        fi
        if (( announced == 0 )); then
            echo "gpu_lock: all ${n} slot(s) busy, waiting..."
            announced=1
        fi
        # Jittered poll. `flock -n` on a handful of files costs microseconds,
        # so poll tightly -- this interval is dead time on a freed GPU. The
        # jitter stops waiters waking in lockstep and thrashing the same slot.
        # Kept in pure bash: an awk/date subprocess per poll would cost more
        # than the wait it is measuring.
        sleep "0.$(( 200 + RANDOM % 600 ))"
    done
}

# nfbd_record_gpu_trace <trace_dir> <process_name>
#
# Write one row describing the GPU(s) this task can see to
# <trace_dir>/<task_hash>.tsv. workflow.onComplete aggregates those rows into
# <outdir>/logs/gpu_trace_<datestamp>.txt, which Nextflow's own trace cannot
# provide: it records what a task cost, never which device it ran on.
#
# Call it after nfbd_acquire_gpu, and also where no lock is taken -- under SLURM
# the scheduler sets CUDA_VISIBLE_DEVICES itself, and that is exactly where
# knowing which model of card a task landed on matters most, because the nodes
# are not identical.
#
# One file per task rather than appends to one shared file. That keeps a resumed
# run's cached tasks in the output (their rows are already on disk) and lets a
# re-run task overwrite its own row and nothing else, with no locking needed.
#
# Diagnostic only, so every failure path returns 0: losing a provenance row must
# never lose a design.
nfbd_record_gpu_trace() {
    local trace_dir="${1:-}"
    local proc="${2:-unknown}"

    [[ -n "${trace_dir}" ]] || return 0
    command -v nvidia-smi >/dev/null 2>&1 || return 0
    mkdir -p "${trace_dir}" 2>/dev/null || return 0

    # The task hash is not available as `task.hash` here: the script text is
    # itself an input to the hash, so at render time Nextflow interpolates the
    # string "null" (verified, not assumed). The work directory encodes it
    # instead -- work/97/dc5a216c1963... is hash 97dc5a216c1963... -- and the
    # first 2+6 characters are what Nextflow prints and puts in the trace file's
    # `hash` column, so the two files join on it.
    local full short
    full="$(basename "$(dirname "${PWD}")")$(basename "${PWD}")" || return 0
    [[ ${#full} -ge 8 ]] || return 0
    short="${full:0:2}/${full:2:6}"

    # An unset CUDA_VISIBLE_DEVICES is written as '-' rather than left empty,
    # matching Nextflow's trace file. A trailing empty TSV field is silently
    # dropped by most readers, which would shift the column off the end of the
    # row for exactly the tasks where "no allocation happened" is the finding.
    #
    # -i accepts indices, UUIDs or PCI bus IDs, so this reports the right cards
    # whether the allocation came from the lock, from SLURM, or from nothing.
    # Note nvidia-smi does not itself honour CUDA_VISIBLE_DEVICES, which is why
    # the value is passed explicitly. Unset means the task can see every GPU on
    # the node, so report them all.
    # nvidia-smi can block indefinitely when the driver wedges, and a task that
    # hangs forever holding a Nextflow slot is worse than one that skips a
    # diagnostic row. Bound it where `timeout` is available.
    local -a smi=(nvidia-smi)
    if command -v timeout >/dev/null 2>&1; then
        smi=(timeout 10 nvidia-smi)
    fi

    local query='index,uuid,name,driver_version,memory.total'
    local out=''
    if [[ -n "${CUDA_VISIBLE_DEVICES:-}" ]]; then
        out=$("${smi[@]}" --query-gpu="${query}" --format=csv,noheader,nounits \
              -i "${CUDA_VISIBLE_DEVICES}" 2>/dev/null) || out=''
    else
        out=$("${smi[@]}" --query-gpu="${query}" --format=csv,noheader,nounits 2>/dev/null) || out=''
    fi
    [[ -n "${out}" ]] || return 0

    # Collapse nvidia-smi's one-line-per-GPU output into a single row, so the
    # trace stays one row per task. n_gpus and the comma-joined device columns
    # cover the multi-GPU-per-task case, which this pipeline does not currently
    # produce but is not prevented from producing.
    local fields=''
    fields=$(printf '%s\n' "${out}" | awk -F' *, *' -v OFS='\t' '
        {
            n++
            idx = idx sep $1; uuid = uuid sep $2; name = name sep $3; mem = mem sep $5
            drv = $4
            sep = ","
        }
        END { if (n) print n, idx, uuid, name, drv, mem }') || return 0
    [[ -n "${fields}" ]] || return 0

    # Write to a temporary file and rename. rename is atomic within a
    # directory, so a task killed mid-write leaves either the old row or the
    # new one, never a truncated line for the aggregator to trip over.
    local tmp="${trace_dir}/.${full}.$$.tmp"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        "${short}" \
        "${proc}" \
        "$(hostname 2>/dev/null || echo "${HOSTNAME:-unknown}")" \
        "${fields}" \
        "${CUDA_VISIBLE_DEVICES:--}" \
        > "${tmp}" 2>/dev/null || { rm -f "${tmp}" 2>/dev/null; return 0; }
    mv -f "${tmp}" "${trace_dir}/${full}.tsv" 2>/dev/null || rm -f "${tmp}" 2>/dev/null
    return 0
}
