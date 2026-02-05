#!/bin/bash

# Based on settings from: https://doi.org/10.1101/2025.07.23.666285

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

# TODO: We need support for filter sets and use: peptide_filters_relaxed
nextflow run ${PIPELINE_DIR}/main.nf \
  --method bindcraft \
  -params-file params.json \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
