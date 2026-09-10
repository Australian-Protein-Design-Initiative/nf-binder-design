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
