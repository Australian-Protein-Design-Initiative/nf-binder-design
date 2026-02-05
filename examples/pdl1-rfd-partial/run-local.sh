#!/bin/bash

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method rfd_partial \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_partial_per_binder=1 \
  --rfd_batch_size=1 \
  --rfd_partial_T "5,15" \
  --pmpnn_seqs_per_struct=1 \
  --pmpnn_relax_cycles=3 \
  --rfd_filters="rg<=20" \
  --refold_af2ig_filters="pae_interaction<=10;plddt_binder>=80" \
  --af2ig_recycle=3 \
  --refold_max=100 \
  --refold_use_msa_server=true \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
