#!/bin/bash
set -euo pipefail

# The site configs under conf/platforms/ still use top-level `def`, which Nextflow
# >=26's default (strict) parser rejects, so ask for the v1 parser. Pinning an older
# Nextflow instead does not work any more: the workflow scripts no longer compile
# under 24.10.0, which rejects the `def` bindings in a workflow body with
# "Variable `log` already defined in the process scope".
export NXF_SYNTAX_PARSER=v1

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)
DEFAULT_SLURM_ACCOUNT=$(sacctmgr --parsable2 show user -s ${USER} | tail -1 | cut -f 2 -d \|)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.m3.config \
  --method fold \
  --slurm_account ${DEFAULT_SLURM_ACCOUNT} \
  --input 'input/pdl1.fasta' \
  --outdir results \
  --methods af2,boltz,rf3,protenix,openfold3,esmfold2,esmfold2_fast \
  --msa_method jackhmmer_af2 \
  --n_predictions 1 \
  --af2_keep_models best \
  -profile slurm,m3 -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt \
  "$@"
