#!/bin/bash
set -euo pipefail
if [ "$(id -gn)" != "alphafold" ]; then exec sg alphafold -c "$0 $*"; fi

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.m3.config \
  --method fold_pulldown \
  --targets input/targets.fasta \
  --binders input/binders.fasta \
  --outdir results \
  --methods af2,boltz,rf3,protenix \
  --msa_method mmseqs2_colabfold \
  --use_remote_server true \
  --create_target_msa true \
  --create_binder_msa false \
  --n_predictions 5 \
  --af2_keep_models all \
  -profile local -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt \
  "$@"
