# Multiple GPUs on one machine

Nextflow has no built-in notion of a GPU as a consumable resource on the `local`
executor. If two tasks start at once on a two-GPU workstation, nothing in
Nextflow stops both from choosing the same card. This page describes how the
pipeline solves that, and how to configure it.

## Quick start

Pass the GPUs the pipeline may use, and select a local profile:

```bash
nextflow run main.nf --gpu_devices 0,1 -profile local -c nextflow.dual-gpu.config
```

That is the whole configuration. Tasks claim a GPU on startup, hold it for their
lifetime, and release it when they exit.

Leaving `--gpu_devices` unset (the default) disables GPU management entirely,
which is what you want under SLURM or any other scheduler that sets
`CUDA_VISIBLE_DEVICES` itself.

!!! tip "Spreading work across GPUs is only half the problem"

    Roughly 25 seconds of every RosettaFold3 task is Python imports and
    checkpoint loading before any GPU work starts, and that cost is paid once
    per task no matter which card it lands on. On a `-profile local` multi-GPU
    run, **raising batch sizes buys more throughput than any amount of
    scheduling.**

    - `--rf3_batch_size` groups already-generated sequences into one `rf3 fold`
      call. This is the one to reach for first: it amortises the startup cost
      with **no effect on design diversity** and no extra VRAM. At batch size 4
      an rfd3 task spends about 7 s of overhead per design instead of 28 s.
    - `--mpnn_batch_size` likewise generates several sequences per backbone in
      one call.
    - `--rfd3_batch_size` batches too, **but costs diversity**: every design in
      an RFdiffusion3 batch is sampled at the *same length*, so `--rfd3_batch_size 2`
      over 100 designs gives 50 independent length draws rather than 100. Raise
      it only if you are willing to trade binder length diversity for speed.

    The trade-off to know about: one failing design fails its whole batch.

## How allocation works

Each GPU task sources [`bin/gpu_lock.sh`](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/blob/main/bin/gpu_lock.sh)
and calls `nfbd_acquire_gpu`. That function keeps one lock file per GPU slot in
`params.gpu_lock_dir` (by default `<workDir>/.gpu_locks`) and claims one with
`flock`.

Two properties make this reliable where earlier approaches were not:

- **The claim is atomic.** `flock` either grants the lock or it does not, so two
  tasks starting in the same millisecond cannot both win the same card.
- **The kernel owns the release.** The lock lives on an open file descriptor
  held by the task shell. When the task exits the descriptor closes and the lock
  drops -- including when the task is killed, or the machine loses power. There
  are no stale lock files to clean up and no PID bookkeeping to go wrong.
- **Slots fill breadth-first.** The scan tries slot 0 of every GPU before slot 1
  of any of them, so a task always prefers an idle card over a second slot on a
  busy one. The GPU order is rotated per task so simultaneous starts do not all
  probe the same card first.

If every slot is busy the task polls until one frees, up to
`params.gpu_lock_timeout` seconds, then fails.

## What it is worth

Measured on a dual RTX 3060 workstation, running the same rfd3 workload twice
with **only the allocation mechanism changed** -- identical `maxForks`,
`memory`, `queueSize`, `pollInterval` and `submitRateLimit`:

| | static round-robin | flock allocation |
| --- | --- | --- |
| wall time | 967 s | 960 s |
| both GPUs working | 15.9% | **27.0%** |
| one card stacked while the other idles | 15.7% | **3.9%** |
| mean GPU utilisation | 21.3% | 23.1% |

Wall time is unchanged; what changes is that the work is spread properly. The
benefit grows with how contended the machine is -- on a full 800-design run the
static round-robin left one card stacked and the other idle for **41% of the
wall clock**.

Note the low absolute utilisation in both columns. That is not a scheduling
failure: roughly 25 s of every ~64 s task is Python imports and checkpoint
loading before any GPU work starts. Batching is the lever for that, not
scheduling -- see the throughput notes below.

## Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--gpu_devices` | `''` | GPUs this pipeline may claim: `0,1`, or `all`. Empty disables GPU management. |
| `--gpu_slots_per_device` | `1` | Concurrent tasks allowed per GPU. Override per process with `ext.gpu_slots`. |
| `--gpu_lock_timeout` | `14400` | Seconds a task waits for a free GPU before failing. |
| `--gpu_lock_dir` | `<workDir>/.gpu_locks` | Where lock files live. Must be shared by all tasks and visible inside containers. |
| `--gpu_trace_dir` | `<workDir>/.gpu_trace` | Where per-task GPU records are written before aggregation. Same constraints as `--gpu_lock_dir`. |
| `--gpu_trace_file` | `<outdir>/logs/gpu_trace_<datestamp>.txt` | Aggregated GPU trace for the run. |

Every process that runs a model claims a GPU, including MPNN. That matters
because MPNN's inference engine falls back to `torch.device("cuda")`, i.e. the
first visible device, so without a claim every MPNN task lands on GPU 0
regardless of what the heavier processes are doing.

For BoltzGen, note that `--devices` is a *count* of GPUs, not an index. A task
pinned to one card sees exactly one device, so leave `--devices` at 1 (or
unset); asking for 2 while the allocator has granted one will fail.

### What `gpu_slots_per_device` is for

It decides **how many tasks may share a card**, and it defaults to `1` so that a
run can never fail from two tasks competing for VRAM on the same device.

It is a VRAM cap, not the way to control concurrency. The lock decides *where* a
task runs; `maxForks` and the executor's memory budget decide *how many* run.
Conflating the two is the easiest way to make things slower -- if `maxForks`
exceeds the total slot count, the surplus tasks are dispatched only to block on
the lock while holding their memory reservation.

Raising it is worth doing for processes whose VRAM footprint is small relative
to the card. A task holds its GPU allocated but idle for roughly 25 s of Python
imports and checkpoint loading, and a second task can fill that gap. Because
that headroom is a property of the *process*, not the run, override it per
process rather than globally:

```groovy
process {
    // Measured peak VRAM ~3.4 GB on a 12 GB card, so two fit with room to spare.
    withName: RFDIFFUSION3 { ext.gpu_slots = 2; maxForks = 4 }
    withName: MPNN         { ext.gpu_slots = 2; maxForks = 4 }
}
```

Raise that process's `maxForks` alongside it, to `gpu_slots x nGPU`, or the
extra slots simply go unused.

One asymmetry to keep in mind: **`maxForks` is per process, but slots are
global.** Several GPU processes are usually eligible at once -- RFDIFFUSION3,
ROSETTAFOLD3 and BOLTZ_COMPARE_* all run in the same phase -- so what competes
for slots is the *sum* of their in-flight tasks, not any one `maxForks`. If that
sum exceeds the slot count the surplus blocks on the lock. In practice the
executor's `memory` budget usually binds first, which is another reason to
declare process memory accurately.

Check the VRAM headroom before raising it for a process. Host RAM is usually the
real ceiling anyway: a RosettaFold3 task peaks around **11.5 GB of host RAM**
against 3.4 GB of VRAM, so a 26 GB executor budget already caps concurrency at
two heavy tasks whatever the slot count says. Declare process `memory` honestly
-- under-declaring lets Nextflow dispatch more tasks than fit and the machine
swaps, which costs far more than an idle GPU.

## Recording which GPU ran what

Nextflow's `trace_*.txt` records what every task cost, but never which device it
ran on, and there is no way to add that as a trace column: trace observers run
in the head process, while the device is chosen inside the container. So each
task records it itself.

Every GPU task writes one row to `<outdir>/logs/gpu_trace_<datestamp>.txt`,
which shares its datestamp with the `trace_`, `report_` and `timeline_` files
from the same run:

```
timestamp             hash       process            hostname  n_gpus  gpu_index  gpu_uuid        gpu_name                 driver_version  memory_total_mib  cuda_visible_devices
2026-09-10T07:04:27Z  92/031dcb  RFD3:ROSETTAFOLD3  vrboxen   1       0          GPU-7f0d2332..  NVIDIA GeForce RTX 3060  580.173.02      12288             0
2026-09-10T07:04:27Z  92/235fbd  RFD3:MPNN          vrboxen   1       1          GPU-86d6fdb2..  NVIDIA GeForce RTX 3060  580.173.02      12288             1
```

The `hash` column is the same short hash Nextflow prints and puts in the trace
file, so the two join directly:

```bash
join -1 1 -2 2 -t $'\t' \
  <(tail -n +2 results/logs/trace_*.txt     | sort -k2,2 | awk -F'\t' -v OFS='\t' '{print $2,$4,$9}') \
  <(tail -n +2 results/logs/gpu_trace_*.txt | sort -k2,2)
```

Notes on how it behaves:

- **It records regardless of who chose the device.** The row is written whether
  the GPU came from the lock, from SLURM, or from nothing at all. On a
  heterogeneous cluster that is the point -- the model of card a task landed on
  explains a lot of otherwise puzzling variance in `realtime`.
- **`cuda_visible_devices` is `-` when nothing set it.** The task could see
  every GPU on the node and PyTorch picked one unsupervised, so `gpu_index`
  lists them all and `n_gpus` counts them all. A column of `-` under
  `-profile local` means `--gpu_devices` was never passed.
- **`-resume` keeps the history.** Records live in the work directory, one file
  per task hash, so a resumed run still reports the GPU its cached tasks ran on
  when they actually executed. A re-run task overwrites only its own record.
- **It cannot fail a task or a run.** This is diagnostic output, so every
  failure degrades to a missing row. The recording call is `|| true`, which in
  bash also suspends `set -e` for the whole function body, so nothing inside it
  can abort the script. Tested: `nvidia-smi` absent, returning garbage, exiting
  non-zero, or hanging forever; a bogus `CUDA_VISIBLE_DEVICES`; an unwritable,
  missing or non-directory trace path; and no `hostname` command. In each case
  the task ran to completion.

    A wedged driver is the one worth calling out. `nvidia-smi` can block
    indefinitely, and a task hanging forever while holding a Nextflow slot is
    worse than one that skips a row, so the call is wrapped in `timeout 10`.

    Aggregation is equally contained. An exception escaping
    `workflow.onComplete` is reported as `Failed to invoke workflow.onComplete
    event handler`, which makes a successful run look failed, so the whole
    handler is wrapped and degrades to a warning. Records that cannot be read
    are skipped individually, and rows whose column count does not match the
    header are dropped rather than written into the output.

!!! warning "One thing here is deliberately fatal"

    If `--gpu_devices` is set and `bin/gpu_lock.sh` cannot be sourced, the task
    fails. A GPU claim was asked for and could not be made, and running anyway
    would put tasks on the same card. Only the recording is optional; the
    allocation is not.

## Why not detect busy GPUs with `nvidia-smi`?

Because under containers it cannot see the tasks. Apptainer is invoked with
`--pid`, so each container has its own PID namespace and
`nvidia-smi --query-compute-apps` inside it reports `No running processes found`
however loaded the card is. Across one 800-design run, 250 of 250 task logs
reported no processes while the cards each held 3.2 GB.

Reading `memory.used` instead does not work either: a task spends roughly 25
seconds importing Python modules and loading its checkpoint before it allocates
any VRAM, so every task that starts inside that window sees the same idle card
and takes it.

A static round-robin (`CUDA_VISIBLE_DEVICES = task.index % nGPU`) is worth
understanding because it looks correct and is not. It assigns by *submission*
order, which drifts out of phase with *completion* order as soon as task
durations vary: tasks 1 and 2 take GPU 0 and 1, task 1 finishes and task 3 takes
GPU 0, task 3 finishes and task 4 takes GPU 1 -- where task 2 is still running.
Measured over a full run, that put two tasks on one card while the other sat
completely idle for 41% of the wall clock, even though each card received
exactly the same *number* of tasks.

## Throughput notes

Anything that delays dispatch is idle GPU time, so the local profiles avoid it:

- **No `submitRateLimit`.** Throttling launches caps dispatch regardless of how
  free the GPUs are -- at `'1/10sec'`, six a minute. It was only ever there to
  space out the old VRAM probe, and the lock has no race to space out.
- **Short `pollInterval`.** This is the delay between a task finishing and the
  next being dispatched, which lands directly on a card that just went free.
  The local profiles use 1-2 seconds rather than the 30 seconds some earlier
  configs used.
- **`maxForks` is flow control, not safety.** Exclusivity comes from the lock.
  `maxForks` only stops Nextflow dispatching a crowd of tasks that would sit
  blocked on the lock holding memory reservations. Setting it one above the GPU
  count lets the next task run its startup while the current one still owns the
  card, so the GPU is handed over warm.
- **Batch sizes matter more than any of this.** See the note at the top of the
  page -- per-task startup dominates, and batching is the only thing that
  amortises it.

## Using HyperQueue instead

[HyperQueue](https://it4innovations.github.io/hyperqueue/stable/) looks like the
right tool for this: it models GPUs as indexed resources, hands each task a
distinct one, and will not dispatch a task until a GPU is actually free -- so
unlike the lock, no Nextflow slot ever sits blocked. `conf/platforms/hyperqueue.config`
configures it.

It does allocate correctly. The problem is that **the allocation cannot reach
the workload**, because it is communicated through the environment and this
pipeline runs everything in containers. There are two independent barriers:

1. `apptainer.runOptions` includes `--cleanenv`, which strips
   `CUDA_VISIBLE_DEVICES` on the way into the container.
2. Nextflow's Apptainer launcher wraps the call in `env -`, wiping the
   environment down to a fixed allowlist (`PATH` and a few `APPTAINERENV_*`
   variables). So even exporting `APPTAINERENV_CUDA_VISIBLE_DEVICES` in a
   `beforeScript` does not survive to the `apptainer exec`.

The result is that HyperQueue assigns GPU 1 to a task, the container never hears
about it, the task sees both cards and PyTorch takes GPU 0 -- so everything
piles onto one card even though the scheduler did its job. Nothing in the
HyperQueue configuration can fix this; the value is discarded downstream of it.

Two notes on things that are *not* the problem, since both look suspicious:

- Nextflow's HyperQueue executor does translate the `accelerator` directive into
  `--resource gpus=N`, while HyperQueue's auto-detected resource is named
  `gpus/nvidia`. A task asking for `gpus` waits forever with no error. The
  configs here already avoid this by setting `accelerator = null` and requesting
  the resource through `clusterOptions` instead, which the executor appends
  verbatim.
- Nextflow logs `The support for HyperQueue is an experimental feature and it
  may change in a future release`. That is expected.

### Making HyperQueue work

Use HyperQueue for queueing and the lock for device selection. They compose
cleanly: HyperQueue ensures no more GPU tasks are dispatched than there are
GPUs, and `gpu_lock.sh` -- which runs *inside* the container and so is immune to
both barriers above -- decides which card each task takes.

```bash
hq server start &
hq worker start &

nextflow run main.nf -profile hyperqueue --gpu_devices 0,1
```

Passing `--gpu_devices` is what enables the in-container allocator. Without it,
tasks fall back to whatever PyTorch chooses, which is GPU 0 for all of them.

The lock-based default needs no daemon, no version coordination between server
and worker, and no second scheduler to reason about, which is why it is the
default.
