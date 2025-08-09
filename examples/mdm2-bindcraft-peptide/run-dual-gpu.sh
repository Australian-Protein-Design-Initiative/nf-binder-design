#!/bin/bash

# This example is designed to run on a local machine with 2 GPUs.
# We use the --gpu_devices=0,1 flag to specify the GPUs to use.
# Don't use --gpu_devices=0,1 with SLURM - the scheduler should handle this for you.

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/bindcraft.nf  \
  -c nextflow.dual-gpu.config \
  -params-file params.json \
  --gpu_devices=0,1 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
