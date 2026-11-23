#!/bin/bash

# Same design settings as run-local-simple.sh, one generation task per GPU.
# Four designs at two per task is two tasks, so both cards stay busy.
# Don't use --gpu_devices with SLURM - the scheduler assigns the device.

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.docker.config \
  -c nextflow.dual-gpu.config \
  --method rfd4_proteina \
  --design_name pdl1 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --contigs "50-120, /0, A18-132" \
  --hotspot_res "A54,A56,A115,A123,A125" \
  --rfd4_n_designs=4 \
  --rfd4_designs_per_task=2 \
  --rfd4_recycle_strategy 'medium/front_loaded' \
  --gpu_devices "0,1" \
  -profile local \
  -resume
