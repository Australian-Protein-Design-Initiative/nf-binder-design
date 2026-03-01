#!/bin/bash

# Simple example using --contigs / --hotspot_res params (auto-generates rfd3 JSON config)

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method rfd3 \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 50-120]" \
  --hotspot_res "A56" \
  --rfd3_n_designs=8 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
