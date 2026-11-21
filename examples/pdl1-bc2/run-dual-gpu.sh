#!/bin/bash

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method bc2 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-120" \
  --bc2_n_traj 8 \
  --bc2_batch_size 1 \
  --gpu_devices "0,1" \
  -c nextflow.dual-gpu.config \
  -profile local \
  -resume
