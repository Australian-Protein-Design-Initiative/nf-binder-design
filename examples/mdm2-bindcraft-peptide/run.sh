#!/bin/bash

# Based on settings from: https://doi.org/10.1101/2025.07.23.666285

PIPELINE_DIR=../../

# TODO: We need support for filter sets and use: peptide_filters_relaxed
nextflow run ${PIPELINE_DIR}/main.nf \
  --method bindcraft \
  -params-file params.json \
  -profile local \
  -resume
