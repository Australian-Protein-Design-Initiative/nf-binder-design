#!/bin/bash
set -euo pipefail

# Pin Nextflow 24.10.0: site configs under conf/platforms/ still use top-level
# `def`, which Nextflow >=26's default (strict) parser rejects. Use
# NXF_SYNTAX_PARSER=v1 with Nextflow 26, or pin <26 as here.
export NXF_VER=24.10.0

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)
DEFAULT_SLURM_ACCOUNT=$(sacctmgr --parsable2 show user -s ${USER} | tail -1 | cut -f 2 -d \|)

# Mosaic Multispecifics binders x PD-L1 + IL-7Ra.
# ColabFold remote MSAs for targets; binders stay query-only. Every engine,
# AF2 included, gets the target's ColabFold a3m (see README).
nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.m3.config \
  --method fold_pulldown \
  --slurm_account ${DEFAULT_SLURM_ACCOUNT} \
  --targets input/targets.fasta \
  --binders input/binders.fasta \
  --outdir results \
  --methods af2,boltz,rf3,protenix,openfold3,esmfold2,esmfold2_fast \
  --msa_method mmseqs2_colabfold \
  --use_remote_server true \
  --create_target_msa true \
  --create_binder_msa false \
  --n_predictions 5 \
  --af2_keep_models all \
  -profile slurm,m3 -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt \
  "$@"
