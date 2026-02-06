# M3 HPC cluster examples

Here are some examples for specific workflows on the [M3 HPC cluster](https://docs.erc.monash.edu/Compute/HPC/M3/).

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

## RFdiffusion Workflows

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

# Set the apptainer cache directory to a location in your /scratch2 directory (not in /home)
export NXF_APPTAINER_CACHEDIR=/scratch2/ab12/${USER}/apptainer_cache
export APPTAINER_CACHEDIR=$NXF_APPTAINER_CACHEDIR

DATESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p results/logs

nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
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
