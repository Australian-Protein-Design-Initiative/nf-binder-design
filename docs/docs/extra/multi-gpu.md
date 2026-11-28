# Multiple GPUs on one machine

Nextflow has no built-in notion of a GPU as a consumable resource on the `local`
executor. If two tasks start at once on a two-GPU workstation, Nextflow does
not prevent both from choosing the same card. This page describes how the
pipeline coordinates GPU assignment and how to configure it.

## Quick start

Add the `multi_gpu` profile alongside the executor profile:

```bash
nextflow run main.nf -profile local,multi_gpu
```

That is enough on its own. The profile detects the cards with `nvidia-smi`,
turns on GPU management, and caps every GPU process at one task per card.
Tasks allocate an available GPU on startup, retain it during execution, and
release it upon exit.

To restrict a run to particular cards, name them; the profile sizes itself to
whatever you pass:

```bash
nextflow run main.nf -profile local,multi_gpu --gpu_devices 0,1
```

!!! warning "Do not combine `multi_gpu` with a scheduler profile"

    Use it only with `-profile local`. Under SLURM the scheduler assigns the
    card, and the profile's device count is resolved from `nvidia-smi` wherever
    the config is parsed, which on a cluster is the submit node and not the node
    the task runs on.

Leaving `--gpu_devices` unset without the profile (the default) disables GPU
management entirely, which is appropriate under SLURM or any other scheduler
that sets `CUDA_VISIBLE_DEVICES` externally. On a multi-GPU workstation that is
almost never what you want, so a task that finds several cards visible and none
assigned to it prints a warning naming `--gpu_devices` in its task log.

The older `-c nextflow.dual-gpu.config` overlays in `examples/` predate this
profile and still work; they differ only in also setting per-process `memory`
and `cpus` tuned to a specific workstation.

A single-GPU machine also benefits from this configuration. With
`--gpu_devices 0`, the card admits one GPU task at a time, preventing multiple
models from loading concurrently. Without this constraint, runs that execute
multiple fold engines concurrently can trigger out-of-memory errors even when
each engine fits within GPU memory on its own.

!!! tip "Throughput bottlenecks during task initialization"

    Approximately 25 seconds of every RosettaFold3 task is spent importing Python
    modules and loading model checkpoints before GPU execution begins. This
    initialization occurs once per task regardless of card assignment. On a
    `-profile local` multi-GPU run, **increasing batch sizes yields greater
    throughput improvements than scheduling adjustments alone.**

    - `--rf3_batch_size` groups generated sequences into a single `rf3 fold`
      execution. This amortises startup overhead without altering design
      diversity or requiring extra VRAM. At batch size 4, an rfd3 task spends
      about 7 seconds of overhead per design instead of 28 seconds.
    - `--mpnn_batch_size` generates several sequences per backbone in a single call.
    - `--rfd3_batch_size` also batches generation, but restricts diversity:
      every design in an RFdiffusion3 batch is sampled at the *same length*, so
      `--rfd3_batch_size 2` across 100 designs yields 50 independent length draws
      instead of 100. Increase this parameter only when uniform length is acceptable.

    Batch failure trade-off: a failure in one design causes the entire batch task
    to fail.

## Allocation mechanism

Each GPU task sources [`bin/gpu_lock.sh`](https://github.com/Australian-Protein-Design-Initiative/nf-binder-design/blob/main/bin/gpu_lock.sh)
and calls `nfbd_acquire_gpu`. That function maintains one lock file per GPU slot in
`params.gpu_lock_dir` (default: `<workDir>/.gpu_locks`) and claims a slot using
`flock`.

Key design properties:

- Atomic locking: `flock` either grants the lock or rejects it immediately,
  preventing two tasks starting at the exact same millisecond from acquiring the
  same device.
- Kernel-managed release: the lock is tied to an open file descriptor held by
  the task process. When the process terminates, the operating system closes
  the descriptor and releases the lock automatically, even if the task is killed
  or the host system loses power. There are no stale lock files or PID tracking
  tables to clean up.
- Breadth-first slot assignment: the acquisition routine checks slot 0 across
  all GPUs before evaluating slot 1 on any device, preferring an idle GPU over a
  second process on an active card. GPU search order is rotated per task to
  prevent concurrent starts from contending on the first device.

When all slots are occupied, the task polls until a slot becomes available or
until `params.gpu_lock_timeout` seconds elapse, after which the task fails.

The lock directory and lock files are created world-writable because tasks may run
under different user IDs: certain container images run as root, whereas others
run as the invoking user via `containerOptions` (such as `conf/platforms/arm64.config`
under Docker). The first task to execute creates the files. If a subsequent task
cannot access a lock file for writing, it reports an error and exits immediately
rather than waiting for the timeout. In that scenario, delete the lock directory
and restart the run; locks maintain no state across tasks, so deleting the
directory does not lose data.

## Benchmark results

Measured on a dual RTX 3060 workstation, running the same rfd3 workload twice
with identical settings (`maxForks`, `memory`, `queueSize`, `pollInterval`,
`submitRateLimit`) and varying only the allocation mechanism:

| Metric | Static round-robin | flock allocation |
| --- | --- | --- |
| Wall time | 967 s | 960 s |
| Both GPUs active | 15.9% | **27.0%** |
| One GPU active while other idles | 15.7% | **3.9%** |
| Mean GPU utilisation | 21.3% | 23.1% |

While overall wall time remains similar, work distribution across devices improves
substantially. The difference becomes pronounced on heavily contended systems:
on an 800-design run, static round-robin left one card saturated and the other
idle for 41% of total elapsed time.

Overall utilisation numbers reflect Python initialization and checkpoint loading
overhead (roughly 25 seconds of every ~64-second task). Addressing this latency
requires larger batch sizes rather than scheduling changes.

## Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--gpu_devices` | `''` (`all` under `-profile multi_gpu`) | GPUs available to the pipeline: `0,1`, or `all`. Empty disables GPU management. |
| `--gpu_slots_per_device` | `1` | Concurrent tasks allowed per GPU. Can be overridden per process with `ext.gpu_slots`. |
| `--gpu_lock_timeout` | `14400` | Seconds a task waits for an available GPU before failing. |
| `--gpu_lock_dir` | `<workDir>/.gpu_locks` | Directory containing lock files. Must be shared across tasks and mounted inside containers. |
| `--gpu_trace_dir` | `<workDir>/.gpu_trace` | Directory for per-task GPU records before aggregation. Must be shared across tasks and mounted inside containers. |
| `--gpu_trace_file` | `<outdir>/logs/gpu_trace_<datestamp>.txt` | Aggregated GPU trace for the run. |

Every process that executes a neural network model claims a GPU, including MPNN.
This prevents contention: MPNN falls back to `torch.device("cuda")` (the first
visible device), so without explicit allocation every MPNN task executes on GPU 0
regardless of other active workloads.

For BoltzGen, `--devices` specifies a GPU count rather than a device index.
Because each task is pinned to a single card, keep `--devices` set to 1 (or unset);
specifying 2 while only one device is allocated causes an error.

### Managing concurrency with `gpu_slots_per_device`

`--gpu_slots_per_device` defines the maximum number of concurrent tasks assigned to
a single GPU card. It defaults to `1` to avoid VRAM exhaustion from concurrent tasks.

This setting functions as a VRAM guard rather than an overall concurrency limiter.
The file lock determines *where* a task executes, whereas `maxForks` and executor
memory configurations determine *how many* tasks run in total. If `maxForks`
exceeds the available slots, surplus tasks launch and block waiting for locks
while continuing to consume host RAM.

Increasing slot counts is effective for processes with modest VRAM usage relative
to available card capacity. Because task initialization leaves the GPU largely idle
for the first ~25 seconds, a second task can run concurrently without saturating
compute resources. Configure this per process rather than globally:

```groovy
process {
    // Peak VRAM ~3.4 GB on a 12 GB card; two instances run comfortably.
    withName: RFDIFFUSION3 { ext.gpu_slots = 2; maxForks = 4 }
    withName: MPNN         { ext.gpu_slots = 2; maxForks = 4 }
}
```

Increase the corresponding `maxForks` setting to `gpu_slots * nGPU` to ensure the
additional slots are utilised.

Note that `maxForks` is evaluated per process, whereas lock slots are global across
all processes. When multiple GPU processes run concurrently (such as `RFDIFFUSION3`,
`ROSETTAFOLD3`, and `BOLTZ_COMPARE_*`), total slot demand equals the sum of their
in-flight tasks. Set process `memory` directives accurately: under-declaring memory
allows Nextflow to dispatch more tasks than host RAM can support, resulting in swap
thrashing that degrades performance far more than an idle GPU.

## GPU allocation tracing

Nextflow's `trace_*.txt` reports resource usage per task, but does not identify the
assigned GPU device because trace observers execute in the parent process while device
assignment occurs inside the container. To record this information, each task logs
its assigned device directly.

Every GPU task appends one line to `<outdir>/logs/gpu_trace_<datestamp>.txt`, sharing
the datestamp with the run's `trace_`, `report_`, and `timeline_` files:

```
timestamp             hash       process            hostname  n_gpus  gpu_index  gpu_uuid        gpu_name                 driver_version  memory_total_mib  cuda_visible_devices
2026-09-10T07:04:27Z  92/031dcb  RFD3:ROSETTAFOLD3  vrboxen   1       0          GPU-7f0d2332..  NVIDIA GeForce RTX 3060  580.173.02      12288             0
2026-09-10T07:04:27Z  92/235fbd  RFD3:MPNN          vrboxen   1       1          GPU-86d6fdb2..  NVIDIA GeForce RTX 3060  580.173.02      12288             1
```

The `hash` column matches Nextflow's task hash, allowing direct joins with standard
trace files:

```bash
join -1 1 -2 2 -t $'\t' \
  <(tail -n +2 results/logs/trace_*.txt     | sort -k2,2 | awk -F'\t' -v OFS='\t' '{print $2,$4,$9}') \
  <(tail -n +2 results/logs/gpu_trace_*.txt | sort -k2,2)
```

Tracing details:

- Universal logging: a trace record is generated regardless of allocation method,
  whether from the lock script, SLURM environment variables, or default selection.
  On heterogeneous clusters, recording device models helps explain variance in task runtimes.
- Unset environment variables: `cuda_visible_devices` records `-` when unset. In this
  case, the task accesses all host GPUs and PyTorch chooses a default; `gpu_index`
  and `n_gpus` report all detected devices. A column of `-` under `-profile local`
  indicates that `--gpu_devices` was not supplied.
- Resume support: records reside in the work directory (one file per task hash),
  ensuring resumed runs retain the original execution device for cached tasks. Re-run
  tasks overwrite only their corresponding record.
- Fault tolerance: tracing commands never cause a task or workflow to fail. The
  logging call appends `|| true` to prevent shell exit under `set -e`. The invocation
  is wrapped in `timeout 10` to avoid hangs if the driver or `nvidia-smi` stalls.
  Similarly, aggregation inside `workflow.onComplete` catches exceptions and emits a
  warning instead of failing the pipeline. Incomplete or malformed records are discarded.

!!! warning "Lock acquisition failure is fatal"

    If `--gpu_devices` is provided but `bin/gpu_lock.sh` cannot be sourced, the task
    terminates immediately. When GPU locking is enabled, failing to acquire a lock
    would allow uncoordinated GPU sharing. Device allocation is mandatory when
    configured, even though trace logging remains optional.

## Limitations of runtime GPU detection with `nvidia-smi`

Querying active tasks via `nvidia-smi` inside containers is ineffective: Apptainer
launches containers with `--pid`, isolating the container PID namespace. As a result,
`nvidia-smi --query-compute-apps` reports no running processes even under heavy load.

Checking `memory.used` dynamically is similarly unreliable: tasks spend approximately
25 seconds loading Python libraries and weights before allocating VRAM. Multiple tasks
launching within that initialization window observe identical idle VRAM readings and
contend for the same card.

Static round-robin allocation (`CUDA_VISIBLE_DEVICES = task.index % nGPU`) also fails
in practice. Because task execution times vary, submission order diverges quickly from
completion order: tasks 1 and 2 take GPUs 0 and 1; task 1 finishes and task 3 claims
GPU 0; task 3 finishes and task 4 claims GPU 1 while task 2 is still running on it.
Across full runs, this imbalance leaves cards idle while others run multiple concurrent jobs.

## Throughput considerations

To maximize GPU utilisation on local systems, local execution profiles use the
following settings:

- Unconstrained dispatch: `submitRateLimit` is disabled. Launch throttling is
  unnecessary because lock acquisition manages concurrency directly.
- Low polling interval: `pollInterval` is set to 1-2 seconds (rather than 30 seconds),
  allowing tasks to launch immediately when a card becomes available.
- Flow control via `maxForks`: exclusivity is guaranteed by `flock`. Setting
  `maxForks` to one above total GPU count allows the subsequent task to start its
  initialization while the previous task completes GPU compute.
- Batch size optimization: per-task initialization dominates execution time; grouping
  multiple designs per task amortises this overhead substantially.

## HyperQueue integration

[HyperQueue](https://it4innovations.github.io/hyperqueue/stable/) provides resource-aware
job queuing, tracking GPUs as indexed resources and dispatching tasks only when devices
become free. `conf/platforms/hyperqueue.config` configures this integration.

However, container execution presents a challenge for direct HyperQueue GPU assignment:
environment variables set by HyperQueue are stripped by container isolation.
Specifically:

1. `apptainer.runOptions` specifies `--cleanenv`, which removes `CUDA_VISIBLE_DEVICES`.
2. Nextflow's Apptainer integration executes commands via `env -`, preserving only a
   strict allowlist of environment variables. Setting `APPTAINERENV_CUDA_VISIBLE_DEVICES`
   in a `beforeScript` does not persist into the container invocation.

Without additional coordination, the container cannot read HyperQueue's device assignment,
and PyTorch defaults to GPU 0 for all concurrent tasks.

### HyperQueue with containerized lock allocation

HyperQueue and the file lock mechanism can be combined effectively: HyperQueue manages
task dispatch and queue concurrency, while `gpu_lock.sh` executes inside the container
to determine device assignment:

```bash
hq server start &
hq worker start &

nextflow run main.nf -profile hyperqueue --gpu_devices 0,1
```

Supplying `--gpu_devices` enables the in-container allocation script.

The default lock-based setup requires no background daemons, worker coordination, or
external scheduler configuration, making it suitable for standard local execution.
