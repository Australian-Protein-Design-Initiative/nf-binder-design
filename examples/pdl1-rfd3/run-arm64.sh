#!/bin/bash

# ARM64 / NVIDIA GB10 (Blackwell) test run.
# Simple example using --contigs / --hotspot_res params (auto-generates rfd3 JSON config)

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow-arm64.config \
  --method rfd3 \
  --design_name pdl1 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --contigs "[A18-132/0 50-120]" \
  --hotspot_res "A56" \
  --rfd3_n_designs=3 \
  --rfd3_filters="rg<=20" \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
