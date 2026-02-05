#!/bin/bash

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method boltzgen \
  --config_yaml 2vsm-cyclic.yaml \
  --outdir results \
  --design_name 2vsm-cyclic \
  --protocol peptide-anything \
  --num_designs 10 \
  --batch_size 1 \
  --budget 5 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
