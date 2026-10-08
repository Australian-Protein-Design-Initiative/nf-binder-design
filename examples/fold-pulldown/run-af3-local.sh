#!/bin/bash
set -euo pipefail
# AlphaFold3 + Protenix target x binder pulldown (ColabFold remote MSAs for targets).

# AlphaFold3 weights are not bundled - fetch them once with
# ../../models/download_af3_weights.sh (default location), or set AF3_MODEL_DIR.
AF3_MODEL_DIR=${AF3_MODEL_DIR:-$(cd ../.. && pwd)/models/alphafold3}

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.m3.config \
  --method fold_pulldown \
  --targets input/targets.fasta \
  --binders input/binders.fasta \
  --outdir results-af3 \
  --methods af3,protenix \
  --af3_model_dir "${AF3_MODEL_DIR}" \
  --msa_method mmseqs2_colabfold \
  --use_remote_server true \
  --create_target_msa true \
  --create_binder_msa false \
  --n_predictions 5 \
  -profile local -resume \
  -with-report "results-af3/logs/report_${DATESTAMP}.html" \
  -with-trace "results-af3/logs/trace_${DATESTAMP}.txt" \
  "$@"
