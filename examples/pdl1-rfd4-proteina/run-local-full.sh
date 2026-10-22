#!/bin/bash

# Extended example: more designs, structural filtering, Boltz-2 refolding and FoldSeek search.

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.docker.config \
  --method rfd4_proteina \
  --rfd4_config pdl1_rfd4.yaml \
  --outdir results \
  --rfd4_n_designs=20 \
  --rfd4_designs_per_task=5 \
  --rfd4_batch_size=1 \
  --rfd4_recycle_strategy='thinking_medium' \
  --rfd4_filters="rg<=20" \
  --full_refold_with 'boltz' \
  --full_refold_filter_sort='pair_pae_min' \
  --full_refold_max=5 \
  --full_refold_create_target_msa \
  --full_refold_use_msa_server \
  --do_foldseek \
  --foldseek_database=CATH50 \
  --foldseek_databases_path="$(pwd)/../databases/foldseek" \
  -profile local \
  -resume
