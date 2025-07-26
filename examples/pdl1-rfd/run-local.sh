#!/bin/bash

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf  \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "[A56]" \
  --rfd_n_designs=4 \
  --rfd_batch_size=1 \
  --pmpnn_seqs_per_struct=2 \
  --pmpnn_relax_cycles=1 \
  --rfd_filters="rg<=20" \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
