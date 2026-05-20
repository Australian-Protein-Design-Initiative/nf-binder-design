# Pawsey Setonix

> ⚠️ Note - the Pawsey Setonix supercomputer uses AMD MI250X GPUs - these do not support CUDA.
> As a result, only the `rfd` workflow is currently supported on Setonix.
> Other workflows may be ported and supported in the future. however this requires tools and containers to be modified to support ROCm, which is not trivial in many cases.

Run **nf-binder-design** on [Pawsey Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925914/Nextflow) using the `pawsey_setonix` profile. GPU processes use AMD MI250X (ROCm); see `conf/platforms/pawsey_setonix.config` for per-process container and SLURM settings.

## Workflow nodes

SSH to a workflow node: `setonix-workflow.pawsey.org.au` (`setonix-07`, `setonix-08`). These nodes are intended for long-running workflows like Nextflow that submit work to the SLURM queue.

Follow [How to Run Workflows on the Workflow Nodes](https://pawsey.atlassian.net/wiki/spaces/US/pages/286097469/How+to+Run+Workflows+on+the+Workflow+Nodes) and start a `screen` or `tmux` session before launching the pipeline.

## Setup

Load Singularity and Nextflow (versions on Setonix may differ; check `module avail`):

```bash
module load singularity/4.1.0-mpi-gpu
module load nextflow/25.04.6
```

Do **not** `module load rocm` in the shell that launches Nextflow. Host ROCm paths in `LD_LIBRARY_PATH` can break containers that are not ROCm images (even `/bin/bash` inside the container).

`unset SBATCH_EXPORT` so SLURM variables from an interactive or batch login shell are not passed into child jobs in a way that confuses the executor.

### Singularity cache and temp directories

The profile uses `$MYSOFTWARE/.nextflow_singularity` as the default image cache. Set these explicitly if you prefer a fixed path (replace project and username):

```bash
export SINGULARITY_CACHEDIR=/software/projects/${PAWSEY_PROJECT}/${USER}/.nextflow_singularity
export NXF_SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR}

export SINGULARITY_TMPDIR=/scratch/${PAWSEY_PROJECT}/${USER}/tmp
export NXF_SINGULARITY_TMPDIR=${SINGULARITY_TMPDIR}
mkdir -p "${SINGULARITY_TMPDIR}"
```

Use **scratch** for `SINGULARITY_TMPDIR` so large image layers and build temps do not fill project quota under `/software/projects`.

### GPU SLURM account

GPU jobs on Setonix must use the project account with a **`-gpu` suffix** (for example `pawsey1343-gpu`, not `pawsey1343`). Pass it to Nextflow so each GPU process gets `--account=...` in `clusterOptions`:

```bash
--slurm_account=pawsey1343-gpu
```

## Example run script

From your run directory (with `input/*.pdb` and a local or pulled copy of the pipeline), use a script like:

```bash
#!/bin/bash

module load singularity/4.1.0-mpi-gpu
module load nextflow/25.04.6

unset SBATCH_EXPORT

export SINGULARITY_CACHEDIR=/software/projects/${PAWSEY_PROJECT}/${USER}/.nextflow_singularity
export NXF_SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR}
export SINGULARITY_TMPDIR=/scratch/${PAWSEY_PROJECT}/${USER}/tmp
export NXF_SINGULARITY_TMPDIR=${SINGULARITY_TMPDIR}
mkdir -p "${SINGULARITY_TMPDIR}"

DATESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p results/logs

# Local git checkout (development or pinned revision)
PIPELINE_DIR=./nf-binder-design/

nextflow run ${PIPELINE_DIR}/main.nf \
  --slurm_account=${PAWSEY_PROJECT}-gpu \
  --method rfd \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 45-65]" \
  --hotspot_res "A56" \
  --rfd_n_designs=1 \
  -profile pawsey_setonix \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
```

### Running from the Nextflow hub

Instead of `PIPELINE_DIR`, pull the release and run by name:

```bash
nextflow pull Australian-Protein-Design-Initiative/nf-binder-design || true

nextflow run Australian-Protein-Design-Initiative/nf-binder-design \
  --slurm_account=${PAWSEY_PROJECT}-gpu \
  --method rfd \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=4 \
  -profile pawsey_setonix \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
```

## References

- [Setonix GPU partition quick start](https://pawsey.atlassian.net/wiki/spaces/US/pages/51928618/Setonix+GPU+Partition+Quick+Start)
- [Singularity on Pawsey](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925894/Singularity)
- [Job scheduling (SLURM partitions)](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925964/Job+Scheduling)
- [Nextflow on Pawsey Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925914/Nextflow)
