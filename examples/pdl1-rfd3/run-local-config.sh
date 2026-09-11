#!/bin/bash

# Example using a user-supplied rfd3 JSON config file

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method rfd3 \
  --rfd3_config pdl1_rfd3.json \
  --outdir results \
  --rfd3_n_designs=2 \
  -profile local \
  -resume
