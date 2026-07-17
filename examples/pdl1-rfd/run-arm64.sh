#!/bin/bash

# ARM64 / NVIDIA GB10 (Blackwell) test run.

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow-arm64.config \
  --method rfd \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=4 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
