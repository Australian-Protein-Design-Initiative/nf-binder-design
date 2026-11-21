#!/bin/bash
# Four BindCraft2 trajectories against PD-L1, run one at a time on a single GPU.
#
# Same design settings as run-dual-gpu.sh, so the results are comparable; it
# just takes about twice as long. Use --gpu_devices to pick a different card.

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
  -profile local \
  -resume
