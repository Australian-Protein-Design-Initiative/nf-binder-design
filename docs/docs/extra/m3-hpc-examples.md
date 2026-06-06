# M3 HPC cluster examples

Here are some examples for specific workflows on the Monash [M3 HPC cluster](https://docs.erc.monash.edu/Compute/HPC/M3/).

These can be adapted for other HPC clusters the use SLURM.

## Setup

Load the Nextflow module:

```bash
module load nextflow/24.04.3

nextflow info
```

Pull the workflow:

```bash
nextflow pull Australian-Protein-Design-Initiative/nf-binder-design
```

(this creates a local cache of the workflow in `~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/`)

You can also pull a specific version like: 

```bash
nextflow pull -r 0.2.0 Australian-Protein-Design-Initiative/nf-binder-design
```

Add the `-r` flag (eg `-r 0.2.0`) to examples below when using a specific version of the workflow.

## RFdiffusion Example

Create an SBATCH script named `run.sh` like:

```bash
#!/bin/bash
#SBATCH --account=ab12
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=nf-binder-design
#SBATCH --output=nf-binder-design-%j.log
#SBATCH --error=nf-binder-design-%j.err

# Change this to your project ID
export PROJECT_ID=ab12

# Set the apptainer cache directory to a location in your /scratch2 directory (not in /home)
export NXF_APPTAINER_CACHEDIR=/scratch2/${PROJECT_ID}$/${USER}/apptainer_cache
export APPTAINER_CACHEDIR=$NXF_APPTAINER_CACHEDIR

# When apptainer downloads container images, it often stages them into /tmp (TMDIR) by default.
# Occasionally /tmp on the node fills up, so we set a custom tmp directory to prevent this.
export APPTAINER_TMPDIR=/scratch2/${PROJECT_ID}/${USER}/tmp
export TMPDIR=${APPTAINER_TMPDIR}
mkdir -p ${APPTAINER_TMPDIR}

DATESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p results/logs

module load nextflow/24.04.3

nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --slurm_account ${PROJECT_ID} \
  --method rfd \
  --input_pdb input/PDL1.pdb \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=4 \
  --rfd_batch_size=1 \
  --pmpnn_seqs_per_struct=2 \
  --pmpnn_relax_cycles=1 \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt \
  -resume \
  -profile slurm,m3

# or, if you have access to the `bdi` partition, use:
# -profile slurm,m3_bdi
```

Get the input files for this example (from the `examples/pdl1-rfd` directory):

```bash
mkdir -p input
cp ~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/examples/pdl1-rfd/input/PDL1.pdb input/
```

Submit the job to the queue:
```bash
sbatch run.sh
```

Monitor the output of `nf-binder-design-%j.log`, or `.nextflow.log`.

> If using an interactive `smux` or low resource CPU-only Strudel session, you can run without `sbatch` like: `./run.sh` - Nextflow will still submit jobs to the queue. This can be convenient for debugging.

## BindCraft Example

Create an SBATCH script named `run.sh` like:

```bash
#!/bin/bash
#SBATCH --account=ab12
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=nf-binder-design
#SBATCH --output=nf-binder-design-%j.log
#SBATCH --error=nf-binder-design-%j.err

# Change this to your project ID
export PROJECT_ID=ab12

# Set the apptainer cache directory to a location in your /scratch2 directory (not in /home)
export NXF_APPTAINER_CACHEDIR=/scratch2/${PROJECT_ID}/${USER}/apptainer_cache
export APPTAINER_CACHEDIR=$NXF_APPTAINER_CACHEDIR

# When apptainer downloads container images, it often stages them into /tmp (TMDIR) by default.
# Occasionally /tmp on the node fills up, so we set a custom tmp directory to prevent this.
export APPTAINER_TMPDIR=/scratch2/${PROJECT_ID}/${USER}/tmp
export TMPDIR=${APPTAINER_TMPDIR}
mkdir -p ${APPTAINER_TMPDIR}

DATESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p results/logs

module load nextflow/24.04.3

nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --slurm_account ${PROJECT_ID} \
  --method bindcraft \
  --input_pdb input/PDL1.pdb \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-120" \
  --bindcraft_n_traj 4 \
  --bindcraft_batch_size 1 \
  --bindcraft_advanced_settings_preset "default_4stage_multimer" \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt \
  -resume \
  -profile slurm,m3

# or, if you have access to the `bdi` partition, use:
# -profile slurm,m3_bdi
```

Get the input files for this example (from the `examples/pdl1-bindcraft` directory):

```bash
mkdir -p input
cp ~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/examples/pdl1-bindcraft/input/PDL1.pdb input/
```

Submit the job to the queue:
```bash
sbatch run.sh
```

Monitor the output of `nf-binder-design-%j.log`, or `.nextflow.log`.

> If using an interactive `smux` or low resource CPU-only Strudel session, you can run without `sbatch` like: `./run.sh` - Nextflow will still submit jobs to the queue. This can be convenient for debugging.

## Troubleshooting failures

If you encounter an error, the `.nextflow.log` file will report the failing task and the path into the `./work` directory, like:

```
Work dir:
  /home/harshil/repos/nf-core/fetchngs/work/ab/123456789aabbccddeeff123456789
```

This directory contains:

- `command.log`: contains both stdout and stderr from the task
- `exitcode`: created when the job ends, with exit code
- `command.run`: wrapper script used to run the job (handles environment setup and job submission)
- `command.sh`: command used for this task
- `command.trace`: logs of compute resource usage
- any input files used for the task (often symlinked from an upstream task)
- any output files generated before the error occurred

You can diagnose the problem by inspecting the `command.log` file, the options used in `.command.sh` and the input and output files.

### Common errors

----

**Problem/message:**

```
RuntimeError: CUDA out of memory. Tried to allocate 8.94 GiB (GPU 0; 15.90 GiB total capacity; 8.94 GiB already allocated; 6.34 GiB free; 0 bytes cached)
```
_(or some variation involving `CUDA` and a failure to allocate memory)_

The task requires more GPU memory (VRAM) than was available to it.

**Solution:**

- If possible, trim your target or generate smaller binders – fewer residues means less VRAM usage.
- Create a custom `nextflow.config` based on `~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/conf/platforms/m3.config` to use a GPU with more VRAM (such as an A100 80GB or H100 80GB).

----

**Problem/message:**

`.command.log` ends with:
```
slurmstepd: error: Detected 1 oom-kill event(s) in step 1090990.batch cgroup.
```

or simply ends with:
```
Killed
```

The task exceeded the memory (RAM) allocated to the task.

**Solution:**

- Create a custom `nextflow.config` based on `~/.nextflow/assets/Australian-Protein-Design-Initiative/nf-binder-design/conf/platforms/m3.config` and increase the `memory` directive for the task.

**Problem/message:**

```
no space left on device
```
or:
```
Disk quota exceeded
```

The filesystem where containers are being downloaded or `results` and `work` are being written is full / over quota.

**Solution:**

- Check that the filesystem where `work` and `results` are being written is not over quota.

- Ensure the environment variables `NXF_APPTAINER_CACHEDIR` and `APPTAINER_TMPDIR` are set in your run.sh script or `~/.bashrc` to a location with sufficient space (ie your `/scratch2` directory), eg:

```bash
export PROJECT_ID=ab12
export NXF_APPTAINER_CACHEDIR=/scratch2/${PROJECT_ID}/${USER}/apptainer_cache
export APPTAINER_TMPDIR=/scratch2/${PROJECT_ID}/${USER}/tmp
mkdir -p ${APPTAINER_TMPDIR}
```

- Clean up files in your home directory if it's full - see [M3 docs: "Run over your storage quota?"](https://docs.erc.monash.edu/Compute/HPC/M3/Files/RunOverQuota/).
