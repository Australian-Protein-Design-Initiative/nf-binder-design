#!/bin/bash
set -euo pipefail
# AlphaFold3 + Protenix target x binder pulldown (ColabFold remote MSAs for targets).

# The site configs under conf/platforms/ still use top-level `def`, which Nextflow
# >=26's default (strict) parser rejects, so ask for the v1 parser. Pinning an older
# Nextflow instead does not work any more: the workflow scripts no longer compile
# under 24.10.0, which rejects the `def` bindings in a workflow body with
# "Variable `log` already defined in the process scope".
export NXF_SYNTAX_PARSER=v1

# AlphaFold3 weights are not bundled - fetch them once with
# ../../models/download_af3_weights.sh (default location), or set AF3_MODEL_DIR.
AF3_MODEL_DIR=${AF3_MODEL_DIR:-$(cd ../.. && pwd)/models/alphafold3}

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)
DEFAULT_SLURM_ACCOUNT=$(sacctmgr --parsable2 show user -s "${USER}" | tail -1 | cut -f 2 -d \|)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.m3.config \
  --method fold_pulldown \
  --slurm_account "${DEFAULT_SLURM_ACCOUNT}" \
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
  -profile slurm,m3 -resume \
  -with-report "results-af3/logs/report_${DATESTAMP}.html" \
  -with-trace "results-af3/logs/trace_${DATESTAMP}.txt" \
  "$@"
