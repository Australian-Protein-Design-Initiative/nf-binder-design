#!/bin/bash

# This example is designed to run on a local machine with 2 GPUs.
# We use the --gpu_devices=0,1 flag to specify the GPUs to use.
# Don't use --gpu_devices=0,1 with SLURM - the scheduler should handle this for you.

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.dual-gpu.config \
  --method bindcraft \
  --input_pdb 'input/6aru_cropped.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A412" \
  --binder_length_range "55-120" \
  --bindcraft_n_traj 2 \
  --bindcraft_batch_size 1 \
  --bindcraft_advanced_settings_preset "default_4stage_multimer" \
  --gpu_devices=0,1 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt

#  --binder_length_range "65-120" \
#  --contigs "[F2-23/F84-175/F205-267/0 G91-171/G209-263/0]"
