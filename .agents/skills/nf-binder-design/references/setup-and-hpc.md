# Setup and HPC Configuration

## Table of Contents

- [Install checklist](#install-checklist)
- [Choosing a profile](#choosing-a-profile)
- [Check first](#check-first)
- [Prerequisites](#prerequisites)
- [Installing the Pipeline](#installing-the-pipeline)
- [Container Setup](#container-setup)
- [Local Workstation](#local-workstation)
- [HPC Clusters with SLURM](#hpc-clusters-with-slurm)
- [Environment Variables](#environment-variables)
- [Platform-Specific Profiles](#platform-specific-profiles)
- [Custom HPC configuration](#custom-hpc-configuration)
- [SBATCH Script Template](#sbatch-script-template)
- [Utility Scripts](#utility-scripts)
- [Troubleshooting](#troubleshooting)

---

## Install checklist

Walk through these steps in order when helping a user set up from scratch. Skip steps that already pass.

### 1. Check Java

Nextflow requires **Java 17 or later** ([Seqera installation guide](https://docs.seqera.io/nextflow/install)). Java 11 is no longer supported in recent Nextflow releases.

```bash
java -version
```

If Java is missing or below 17, install a current LTS release (Temurin recommended). Example with [SDKMAN](https://sdkman.io/):

```bash
curl -s https://get.sdkman.io | bash
# open a new shell, then:
sdk install java 17.0.10-tem
java -version
```

On HPC, check modules first: `module avail java` or `module avail OpenJDK`.

### 2. Check or install Nextflow

```bash
nextflow info
```

If that fails, on HPC try `module avail nextflow` / `module load nextflow`, then `nextflow info` again.

Otherwise install the self-installing package ([Seqera docs](https://docs.seqera.io/nextflow/install#self-install)):

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/
export PATH="$PATH:$HOME/.local/bin"
nextflow info
```

Add the `export PATH=...` line to `~/.bashrc` or `~/.zshrc` if needed.

### 3. Check or install the pipeline

See [Check first](#check-first) below, then:

```bash
nextflow pull Australian-Protein-Design-Initiative/nf-binder-design
nextflow run Australian-Protein-Design-Initiative/nf-binder-design --help
```

### 4. Smoke test a method

Confirm the pipeline resolves and prints method help:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd --help
```

### 5. Choose and verify the execution profile

Use [Choosing a profile](#choosing-a-profile) below. Run a minimal test with the correct `-profile` (and `-c` / `--slurm_account` / site flags as needed) before a production-scale job.

---

## Choosing a profile

| Your environment | Profile | Notes |
|------------------|---------|-------|
| Single machine with local NVIDIA GPU(s) | `-profile local` | One pipeline instance at a time on a workstation; see dual-GPU caveats below |
| SLURM cluster (generic) | `-profile slurm` | Add `--slurm_account=YOUR_ACCOUNT` if required |
| Monash M3 | `-profile slurm,m3` or `-profile m3` | BDI partitions: `-profile slurm,m3_bdi` |
| MLeRP | `-profile slurm,mlerp` | |
| NCI Gadi (PBS Pro, not SLURM) | `-profile nci_gadi` | Uses PBS executor; set `PROJECT` env var |
| Pawsey Setonix | `-profile pawsey_setonix` | AMD MI250X; see Pawsey docs in repo |
| HyperQueue (`hq`) | `-profile hyperqueue` | Or `-profile hq,hyperqueue` |
| Your cluster is not listed | `-profile slurm -c conf/platforms/your_site.config` | Copy the closest file from `conf/platforms/` and adapt — see [Custom HPC configuration](#custom-hpc-configuration) |

**Verify:** re-run `nextflow info` and a small test job with your chosen profile before scaling up.

**Common mistakes:**
- Using `-profile local` on an HPC login node (submits all GPU work on the login node)
- Using `--gpu_devices=0,1` with SLURM (the scheduler assigns GPUs; omit this flag)
- Using `-profile slurm` on Gadi (use `nci_gadi` instead)

---

## Check first

Before installing, check if the pipeline is already present at `~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design`.

If it exists:

- Notify the user that the pipeline is already installed.
- Ask if they would like to update it. This path is a git working copy — update with `git pull` inside that directory, or run `nextflow pull Australian-Protein-Design-Initiative/nf-binder-design`.
- Note: you can checkout tags for prior versions or switch to `develop` or other branches if requested.

---

## Prerequisites

- **Java** 17+ ([Seqera requirements](https://docs.seqera.io/nextflow/install#requirements))
- **Nextflow** 23.04+ (verify with `nextflow info`)
- **Apptainer** (for containers; usually pre-installed on HPC clusters)
- **NVIDIA GPU** with CUDA support (AMD GPUs on Pawsey Setonix — see `pawsey_setonix` profile)
- **~60 GB+ storage** for containers

The [Install checklist](#install-checklist) above covers Java and Nextflow installation in full.

## Installing the Pipeline

### Option 1: Nextflow Pull (recommended for users)

```bash
nextflow pull Australian-Protein-Design-Initiative/nf-binder-design

# Pull a specific version:
nextflow pull -r 0.2.0 Australian-Protein-Design-Initiative/nf-binder-design

# Test:
nextflow run Australian-Protein-Design-Initiative/nf-binder-design --help
```

In this case, the pipeline code is cached in `~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/`.

### Option 2: Git Clone (for development)

```bash
git clone https://github.com/Australian-Protein-Design-Initiative/nf-binder-design
cd nf-binder-design
nextflow run main.nf --help
```

---

## Container Setup

Containers are automatically downloaded when running the pipeline. Apptainer is the default container runtime (pre-installed on most HPC clusters).

---

## Local Workstation

Use `-profile local` for single-node execution:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --input_pdb target.pdb \
  ... \
  -profile local \
  -resume
```

### Multi-GPU workstations

Pass `--gpu_devices=0,1` (or `all`) and base the configuration on `examples/*/nextflow.dual-gpu.config`, e.g.:

```bash
nextflow run main.nf --method rfd3 --gpu_devices 0,1 \
  -profile local -c nextflow.dual-gpu.config -resume
```

Each GPU task sources `bin/gpu_lock.sh` and claims a `flock`-backed slot for its lifetime. The claim is atomic and the kernel releases it when the task exits (including on SIGKILL), so tasks cannot collide on a card and there are no stale locks to clean up. Slots fill breadth-first, so an idle card is always preferred over a second slot on a busy one.

Do **not** pass `--gpu_devices` under SLURM: the scheduler assigns GPUs and sets `CUDA_VISIBLE_DEVICES` itself. Leaving it unset (the default) disables the in-pipeline allocator.

Key points when writing a local multi-GPU config:

- **Do not set `submitRateLimit`.** It throttles dispatch regardless of GPU availability. It only ever existed to space out the old `nvidia-smi` VRAM probe, which no longer exists.
- **Keep `pollInterval` short** (1-2 s). It is the delay between a task finishing and the next being dispatched onto the card that just freed.
- **`gpu_slots_per_device` defaults to 1**, one task per card. Override per process with `ext.gpu_slots` for processes with VRAM headroom (RFDIFFUSION3 and MPNN peak around 3.4 GB and 3 GB), and raise that process's `maxForks` to `gpu_slots × nGPU` to match. Raising `maxForks` beyond the total slot count just makes tasks block on the lock while holding a memory reservation.
- **Declare process `memory` honestly.** Host RAM is usually the binding constraint, not VRAM: a RosettaFold3 task peaks around 11.5 GB of host RAM against 3.4 GB of VRAM. Under-declaring lets Nextflow overcommit and the machine swaps.

**Throughput: batch sizes matter more than GPU scheduling.** Roughly 25 s of every RosettaFold3 task is Python imports and checkpoint loading before any GPU work starts, paid once per task. `--rf3_batch_size` and `--mpnn_batch_size` amortise that with no cost to design diversity. `--rfd3_batch_size` also amortises it, **but every design in an RFdiffusion3 batch is sampled at the same length**, so raising it trades binder length diversity for speed - `--rfd3_batch_size 2` over 100 designs gives 50 independent length draws, not 100.

**Which GPU ran what** is recorded to `<outdir>/logs/gpu_trace_<datestamp>.txt`, one row per GPU task: hostname, `n_gpus`, and the GPU index, UUID, model, driver version and memory. It shares its datestamp with `trace_<datestamp>.txt` and joins to it on the `hash` column. This happens automatically, under SLURM as well as locally. It is strictly diagnostic and cannot fail a task: a missing, broken or hung `nvidia-smi`, an unwritable trace directory, or an unreadable record all degrade to a missing row. The one deliberate exception is that if `--gpu_devices` is set and `bin/gpu_lock.sh` cannot be sourced, the task fails, because a requested GPU claim could not be made.

See `docs/docs/extra/multi-gpu.md` for the full explanation, including why `nvidia-smi` cannot be used to detect busy GPUs under containers, and how this compares to HyperQueue.

---

## HPC Clusters with SLURM

Use `-profile slurm` for SLURM execution (this is the default if no profile is specified):

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --slurm_account=ab12 \
  ... \
  -profile slurm \
  -resume
```

If `--slurm_account` is omitted, the pipeline will use the default SLURM account for the user. The account used must have access to the GPU partition(s) specified in the config file(s).

---

## Environment Variables

Set these before running on HPC clusters (if not already in your environment):

```bash
# Apptainer cache — set to a location with sufficient space (not /home)
export APPTAINER_CACHEDIR=/scratch/myproject/apptainer_cache
export NXF_APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR}

# Temporary directory — must have sufficient space
export TMPDIR=/scratch/myproject/tmp
export NXF_TEMP=$TMPDIR
mkdir -p $TMPDIR
```

---

## Platform-Specific Profiles

Site-specific configs in `conf/platforms/` can be activated with `-profile`:

| Profile | Site | Typical `-profile` |
|---------|------|-------------------|
| `m3` | Monash M3 HPC cluster | `slurm,m3` |
| `m3_bdi` | Monash M3 with BDI partitions | `slurm,m3_bdi` |
| `mlerp` | MLeRP HPC cluster | `slurm,mlerp` |
| `nci_gadi` | NCI Gadi (PBS Pro, not SLURM) | `nci_gadi` |
| `pawsey_setonix` | Pawsey Setonix (AMD MI250X GPUs) | `pawsey_setonix` |
| `spartan_a100` | Spartan HPC (University of Melbourne, A100) | `slurm,spartan_a100` |
| `spartan_l40s` | Spartan HPC (University of Melbourne, L40S) | `slurm,spartan_l40s` |
| `hyperqueue` | Generic HyperQueue cluster | `hyperqueue` |

Usage examples:

```bash
-profile slurm,m3          # SLURM + Monash M3
-profile nci_gadi           # NCI Gadi (PBS)
-c conf/platforms/m3.config   # equivalent to loading m3 config directly
```

Copy and adapt these configs for your own HPC cluster. Pull requests for additional clusters are welcome.

---

## Custom HPC configuration

If the user's cluster is not in the table above:

1. Ask for their scheduler (SLURM, PBS Pro, HyperQueue), GPU partition/queue names, account/project ID, and scratch filesystem path.
2. Copy the closest match from `conf/platforms/` (e.g. `m3.config` for SLURM + GPU, `nci_gadi.config` for PBS).
3. Edit:
   - `clusterOptions` / queue names and `--gres=gpu:1` (or site equivalent)
   - `--account=` or project ID
   - Apptainer bind mounts for scratch (`/scratch`, `/project`, etc.)
   - `withName:` resource blocks for GPU processes (`RFDIFFUSION`, `BINDCRAFT`, `BOLTZGEN_DESIGN`, etc.)
4. Lightweight processes (`UNIQUE_ID`, `GET_CONTIGS`, `FILTER_DESIGNS`, …) should use `executor = 'local'`.
5. Test with a minimal run before production scale:

```bash
nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --input_pdb 'input/target.pdb' \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs 2 \
  -c conf/platforms/your_site.config \
  -profile slurm \
  --slurm_account=ab12 \
  -resume
```

See `examples/pdl1-rfd/run-m3-full.sh` and docs [M3 HPC examples](https://australian-protein-design-initiative.github.io/nf-binder-design/extra/m3-hpc-examples/) for worked examples.

---

## SBATCH Script Template

For SLURM-based HPC clusters, wrap the Nextflow command in an SBATCH script that runs on a login/head node:

```bash
#!/bin/bash
#SBATCH --account=ab12
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=nf-binder-design
#SBATCH --output=nf-binder-design-%j.log
#SBATCH --error=nf-binder-design-%j.err

export NXF_APPTAINER_CACHEDIR=/scratch/ab12/${USER}/apptainer_cache
export APPTAINER_CACHEDIR=$NXF_APPTAINER_CACHEDIR

mkdir -p results/logs

nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --method rfd \
  --slurm_account=ab12 \
  --input_pdb input/target.pdb \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=100 \
  --rfd_batch_size=5 \
  -resume \
  -profile slurm
```

Submit: `sbatch run.sh`

Monitor: Check `nf-binder-design-*.log` or `.nextflow.log`

---

## Utility Scripts

Scripts in `bin/` can be run standalone with [uv](https://docs.astral.sh/uv/):

**From a git clone** (repo root):

```bash
cd /path/to/nf-binder-design
uv run bin/get_contigs.py --help
```

**After `nextflow pull`** (Nextflow assets cache):

```bash
uv run ~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/bin/get_contigs.py --help
```

Use scripts from the same pipeline version you run with Nextflow.

| Script | Purpose |
|--------|---------|
| `af2_combine_scores.py` | Combine AF2 scores (mid-run progress) |
| `get_contigs.py` | Extract contigs from a PDB structure |
| `trim_to_contigs.py` | Trim a PDB to specified contigs |
| `renumber_chains.py` | Renumber chain residues |
| `pdb_to_fasta.py` | Extract FASTA from PDB |
| `filter_designs.py` | Design filter plugin system |
| `calculate_shape_scores.py` | Shape-based scoring metrics |
| `merge_scores.py` | Merge scoring tables |

---

## Troubleshooting

→ See `references/nextflow-troubleshooting.md` for monitoring commands, failure diagnosis, and common error patterns (CUDA OOM, OOM kill, symlink loops, disk space, BoltzGen YAML paths).
