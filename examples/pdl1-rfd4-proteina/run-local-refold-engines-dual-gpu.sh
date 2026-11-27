#!/bin/bash

# Refold scoring with several engines: fold each design's binder against the target with
# each one, compare the result with the design, score the rebuilt interface, and rank on
# the cross-engine consensus. Every engine gives the same columns, prefixed with its name.
#

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.dual-gpu.config \
  --method rfd4_proteina \
  --rfd4_config pdl1_rfd4.yaml \
  --outdir results-refold-engines \
  --rfd4_n_designs=8 \
  --rfd4_designs_per_task=4 \
  --rfd4_recycle_strategy 'medium/front_loaded' \
  --refold_with 'esmfold2,rf3' \
  --refold_filter 'rfd4_ppi_pass>=1' \
  --refold_sort_by='-rfd4_oracle_ppi_overall_iptm' \
  --refold_max=2 \
  --n_predictions=1 \
  --use_msa_server \
  --consensus_metric=ipsae \
  -profile local \
  -resume
