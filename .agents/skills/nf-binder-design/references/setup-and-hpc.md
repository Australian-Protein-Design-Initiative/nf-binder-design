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

For multi-GPU workstations, base the configuration on `examples/*/nextflow.dual-gpu.config` (e.g. `examples/pdl1-rfd3/nextflow.dual-gpu.config`). Pass it with `-c nextflow.dual-gpu.config` and set `--gpu_devices=0,1` (or similar). The config must:

- Set a slower than default `submitRateLimit` (e.g. `'1/10sec'` or `'1/2sec'`) to reduce race conditions in the busy-GPU detection preamble for each task.
- Set `maxForks` on GPU processes to the number of available GPUs (see the example configs for `RFDIFFUSION3`, `ROSETTAFOLD3`, `BOLTZ_COMPARE_*`, etc.).

**However**, local multi-GPU execution remains unreliable — race conditions in GPU allocation are common even with tuning. **Prefer single-GPU execution** (`--gpu_devices 0`) with sequential runs unless you are on a SLURM cluster where the scheduler handles GPU assignment. Do not use `--gpu_devices=0,1` with SLURM — the scheduler assigns GPUs.

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

DATESTAMP=$(date +%Y%m%d_%H%M%S)
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
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt \
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
