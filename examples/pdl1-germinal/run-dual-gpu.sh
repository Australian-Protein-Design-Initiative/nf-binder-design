#!/bin/bash

# This example is designed to run on a local machine with 2 GPUs.
# We use the --gpu_devices=0,1 flag to specify the GPUs to use.
# Don't use --gpu_devices=0,1 with SLURM - the scheduler should handle this for you.

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p results/logs

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.dual-gpu.config \
  --method germinal \
  --germinal_config configs/pdl1_vhh.yaml \
  --germinal_pdb_dir pdbs \
  --germinal_experiment_name pdl1_vhh \
  --germinal_n_traj 4 \
  --germinal_batch_size 1 \
  --gpu_devices=0,1 \
  --outdir results \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
