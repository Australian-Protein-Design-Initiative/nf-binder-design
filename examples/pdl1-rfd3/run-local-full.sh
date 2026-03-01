#!/bin/bash

# Simple example using --contigs / --hotspot_res params (auto-generates rfd3 JSON config)

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method rfd3 \
  --design_name pdl1 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --contigs "[A18-132/0 50-120]" \
  --hotspot_res "A56" \
  --rfd3_n_designs=4 \
  --rfd3_filters="rg<=20" \
  --rfd3_is_non_loopy=false \
  --refold_with 'boltz' \
  --refold_filter_sort='pair_pae_min' \
  --refold_max=2 \
  --refold_use_msa_server=true \
  --refold_target_fasta='input/full/3BIK_B.fasta' \
  --refold_target_templates='input/full/' \
  --output_rmsd_aligned=true \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
