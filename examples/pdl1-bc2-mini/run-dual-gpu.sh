#!/bin/bash
# Four BindCraft2 trajectories against PD-L1, two at a time across two GPUs.
#
# Stage budgets are left at the BindCraft2 defaults so the accept rate is
# representative. PD-L1 is an easy target, so four trajectories has a fair
# chance of yielding an accepted design - but it is still a small sample, and
# zero accepted is a plausible outcome rather than a pipeline failure.

set -euo pipefail

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method bc2 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-65" \
  --bc2_n_traj 4 \
  --bc2_batch_size 1 \
  --gpu_devices "0,1" \
  -c nextflow.dual-gpu.config \
  -profile local \
  -resume
