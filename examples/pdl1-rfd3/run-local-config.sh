#!/bin/bash

# Example using a user-supplied rfd3 JSON config file

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method rfd3 \
  --rfd3_config pdl1_rfd3.json \
  --outdir results \
  --rfd3_n_designs=2 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
