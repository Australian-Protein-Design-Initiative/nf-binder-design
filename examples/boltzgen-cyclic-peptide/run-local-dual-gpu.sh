#!/bin/bash

# This example is designed to run on a local machine with 2 GPUs.
# We use the --devices=2 and --num_workers=2 options to specify the number GPUs to use,
# and --batch_size=2 so that designs are even distributed across both GPUs.
# The same settings should work on SLURM if each node has 2 GPUs.
#
# We cannot target to specific GPUs in this example - to do this you should set
# CUDA_VISIBLE_DEVICES in the nextflow.config apptainer settings

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method boltzgen \
  --config_yaml 2vsm-cyclic.yaml \
  --outdir results \
  --design_name 2vsm-cyclic \
  --protocol peptide-anything \
  --num_designs 10 \
  --batch_size 2 \
  --budget 5 \
  --devices 2 \
  --num_workers 2 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt

