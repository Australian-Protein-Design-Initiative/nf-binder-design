#!/bin/bash
# The same VHH campaign as run.sh, two trajectories at a time across two GPUs.

set -euo pipefail

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method bc2 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A54,A56,A66,A115" \
  --bc2_modality VHH \
  --bc2_humanize \
  --bc2_n_traj 4 \
  --bc2_batch_size 1 \
  --gpu_devices "0,1" \
  -c nextflow.dual-gpu.config \
  -profile local \
  -resume
