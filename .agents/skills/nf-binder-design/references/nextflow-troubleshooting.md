# Nextflow Troubleshooting — nf-binder-design

Practical gotchas, monitoring patterns, and failure diagnosis for running `nf-binder-design`.

---

## Running Multiple Methods Simultaneously

- **Nextflow cannot run multiple pipelines in the same project directory concurrently** — they share `.nextflow/cache/` and acquire an exclusive lock.
- Best practice is to keep each run in its own directory:
  - `runs/rfd/{run_name}`
  - `runs/bindcraft/{run_name}`
  - `runs/boltzgen/{run_name}`
- On a local GPU workstation, run a single pipeline instance at a time — running two pipelines simultaneously will cause GPU contention and OOM. Sequential runs or SLURM job scheduling avoids this.
- **Local multi-GPU execution is unreliable** — use `--gpu_devices 0` for single-GPU execution. If you must use multiple GPUs, a special dual-GPU Nextflow config is required (see `examples/` in the pipeline repo), and race conditions in GPU allocation are common. SLURM clusters handle multi-GPU scheduling correctly.

---

## GPU Memory & VRAM

- Monitor GPU usage using parsable flags (best for LLMs to read state):
  - **Local** (`-profile local`):
    ```bash
    nvidia-smi --query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used,memory.free --format=csv,largetable
    ```
  - **SLURM HPC**: attach to the compute node using `srun --overlap` (the `--overlap` flag shares the existing job allocation without blocking it):
    ```bash
    srun --overlap --jobid=<JOBID> nvidia-smi --query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used,memory.free --format=csv,largetable
    ```
  - **Live Watch**: Wrap in `watch -n 5` for a refreshing view:
    ```bash
    watch -n 5 "nvidia-smi --query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used,memory.free --format=csv,largetable"
    ```
- To see per-process memory usage:
  ```bash
  nvidia-smi --query-compute-apps=timestamp,gpu_name,pid,used_memory --format=csv
  ```

---

## Monitoring Pipeline Runs

### Always use `-with-trace` and `-with-report`

Every pipeline launch should include these flags:

```bash
DATESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p results/logs

nextflow run Australian-Protein-Design-Initiative/nf-binder-design --method rfd ... \
  -with-trace "results/logs/trace_${DATESTAMP}.txt" \
  -with-report "results/logs/report_${DATESTAMP}.html" \
  -profile local -resume
```

The **trace file** (TSV) is the primary monitoring tool — per-task status, hashes, and resource usage. The **report** (HTML) provides a visual summary after the run completes.

### Trace file columns

| Col | Field | Use for |
|-----|-------|---------|
| `$1` | `task_id` | Unique task number |
| `$2` | `hash` | Work dir: `work/$2*` |
| `$3` | `native_id` | OS process ID or SLURM job ID |
| `$4` | `name` | Stage name (e.g. `RFD:RFDIFFUSION (14)`) |
| `$5` | `status` | `COMPLETED`, `CACHED`, `FAILED`, or blank (running) |
| `$6` | `exit` | Exit code (non-zero = error) |
| `$8` | `duration` | Wall time |
| `$11` | `%cpu` | CPU usage |
| `$12` | `peak_rss` | Peak memory |

`-resume` marks reused tasks as `CACHED` (not `COMPLETED`). Match both when counting progress.

### One-shot summary

```bash
TRACE=$(ls -t results/logs/trace_*.txt 2>/dev/null | head -1)
echo "=== Completed/Cached ==="
tail -n +2 "$TRACE" | awk -F'\t' '$5=="COMPLETED" || $5=="CACHED"' \
  | awk -F'\t' '{print $4}' | sed 's/ (.*)//' | sort | uniq -c | sort -rn
echo
echo "=== Failed ==="
tail -n +2 "$TRACE" | awk -F'\t' '$5=="FAILED"' \
  | awk -F'\t' '{printf "%-45s hash: %-12s exit: %s\n", $4, $2, $6}'
```

### Watch in a separate terminal (auto-refreshes)

```bash
TRACE=$(ls -t results/logs/trace_*.txt | head -1)
watch -n 30 "tail -n +2 '$TRACE' | awk -F'\t' '\$5==\"COMPLETED\" || \$5==\"CACHED\"' | awk -F'\t' '{print \$4}' | sed 's/ (.*)//' | sort | uniq -c | sort -rn; echo; echo '--- Failed ---'; tail -n +2 '$TRACE' | awk -F'\t' '\$5==\"FAILED\"' | awk -F'\t' '{printf \"%-45s hash: %-12s exit: %s\n\", \$4, \$2, \$6}'"
```

### Live tail for failures

```bash
TRACE=$(ls -t results/logs/trace_*.txt | head -1)
tail -f "$TRACE" | awk -F'\t' '$5=="FAILED" {printf "FAILED: %-45s hash: %-12s exit: %s\n", $4, $2, $6}'
```

### Other useful monitoring commands

```bash
# GPU utilization
nvidia-smi --query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used,memory.free --format=csv,largetable

# Work directory size
du -sh work/

# Count output files (quick progress check without trace)
ls results/*/rfdiffusion/pdbs/ 2>/dev/null | wc -l      # RFD backbones
ls results/*/af2_initial_guess/pdbs/ 2>/dev/null | wc -l # AF2 structures
ls results/*/bindcraft/accepted/ 2>/dev/null | wc -l     # BindCraft accepted
```

---

## Diagnosing Failures

When a task fails, follow this escalation path:

**1. Check the trace file** for the failed task's hash:
```bash
TRACE=$(ls -t results/logs/trace_*.txt | head -1)
tail -n +2 "$TRACE" | awk -F'\t' '$5=="FAILED"' \
  | awk -F'\t' '{printf "%-45s hash: %-12s exit: %s\n", $4, $2, $6}'
```

**2. Read the task's command log** using the hash:
```bash
HASH="d3/1a460c"
WORK_DIR="work/$HASH"*   # glob expansion, fast even with large work dirs

cat "$WORK_DIR/.command.log"   # stdout/stderr
cat "$WORK_DIR/.command.err"
cat "$WORK_DIR/.command.sh"    # the command that ran
cat "$WORK_DIR/.command.trace" # resource usage
```

**3. Check `.nextflow.log` for tracebacks** if the work dir doesn't reveal the issue:
```bash
grep -B2 -A10 "Exception\|Traceback\|FAILED" .nextflow.log | grep -v "There's no process matching"
```

**4. Common error patterns in `.command.log`:**
```bash
grep -r "out of memory\|CUDA\|Cannot allocate" work/*/.command.log 2>/dev/null  # CUDA OOM
grep -rl "Traceback" work/*/.command.log 2>/dev/null                            # Python errors
grep -r "No space left\|Cannot write" work/*/.command.log 2>/dev/null           # Disk full
```

---

## Resuming Failed Runs

- Use `-resume` to skip completed tasks. Nextflow caches results in `work/` by content hash.
- If the run crashed (e.g. OOM), re-run the same command with `-resume` — it will pick up where it left off.
- If the input PDB was corrupted (e.g. a symlink loop), fix the file first, then run **without** `-resume` (the hash will differ and old cached results won't match).
- You can safely increase `--rfd_n_designs` / `--bindcraft_n_traj` / `--num_designs` and re-run with `-resume` to extend a successful run with more designs — cached tasks are reused.

---

## Pipeline Stage Streaming

- Nextflow stages run concurrently (RFdiffusion → Filter → MPNN → AF2 all overlap). You'll see progress on multiple stages at once.
- The display refreshes in-place using ANSI escape codes. Read the **last 10–15 lines** of the buffer — earlier lines may show stale status.

---

## Common Errors

### CUDA Out of Memory

```
RuntimeError: CUDA out of memory. Tried to allocate X GiB
```

- Trim your target or design smaller binders (fewer residues = less VRAM).
- Reduce `--rfd_batch_size` to 1.
- Use a GPU with more VRAM (A100 80 GB, H100) via a custom config.

### Clearing GPU Memory After Failed Runs

After an OOM crash, Python/JAX/PyTorch processes may hold GPU memory, causing subsequent runs to also fail. Before retrying:

```bash
# Find processes using GPUs
nvidia-smi --query-compute-apps=pid,used_memory --format=csv

# Kill lingering processes
kill -9 <pid>

# Verify GPU is clear
nvidia-smi --query-gpu=memory.free --format=csv
```

Always check `nvidia-smi` free memory before relaunching after a failure. If memory is still held by zombie processes, kill them first.

DO NOT indiscriminately kill other GPU processes that are not associated with the workflow.

### OOM Kill (RAM)

```
slurmstepd: error: Detected 1 oom-kill event(s)
```
or task log ends with `Killed`.

Increase the `memory` directive for the failing process in your config.

### Symlink Loops from Stale Work Directories

- Nextflow stages input files as symlinks in `work/`. If a run crashes, these can become self-referencing loops (e.g. `input/target.pdb -> target.pdb`).
- **Symptom**: File exists but has zero bytes, or tools report "empty PDB".
- **Fix**: Delete the symlink and restore the real file, or clean the work directory:
  ```bash
  rm -rf work/ .nextflow/
  ```
- **Prevention**: Run `ls -la input/` before launching to confirm files are not dangling symlinks.

### BoltzGen YAML: Relative Paths Fail

Always use **absolute paths** in BoltzGen YAML config files — BoltzGen runs inside a container with a different working directory:

```yaml
# Wrong
protein:
  filepath: input/target.pdb

# Correct
protein:
  filepath: /home/user/project/input/target.pdb
```

### Container Download Failures

Ensure `NXF_APPTAINER_CACHEDIR` points to a location with >60 GB free space and that the directory exists (see `references/setup-and-hpc.md`).

### `rg<=20` Filter Rejection Rate

`--rfd_filters "rg<=20"` filters by radius of gyration; most RFdiffusion outputs are extended. Expect only ~10–15% to pass. This is by design — it saves downstream compute. For more exploratory runs, relax to `"rg<=25"` or remove it.
