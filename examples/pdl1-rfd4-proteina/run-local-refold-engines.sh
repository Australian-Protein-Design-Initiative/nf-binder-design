#!/bin/bash

# Refold scoring with several engines: fold each design's binder against the target with
# each one, compare the result with the design, score the rebuilt interface, and rank on
# the cross-engine consensus. Every engine gives the same columns, prefixed with its name.
#
# Note: chai1 and openfold3 are amd64 only. On an arm64 host use boltz, esmfold2,
# esmfold2_fast or rf3, with -c ../../conf/platforms/arm64.config.

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.docker.config \
  --method rfd4_proteina \
  --rfd4_config pdl1_rfd4.yaml \
  --outdir results-refold-engines \
  --rfd4_n_designs=20 \
  --rfd4_designs_per_task=5 \
  --refold_with 'boltz,chai1,openfold3' \
  --refold_max=10 \
  --n_predictions=3 \
  --use_msa_server \
  --consensus_metric=ipsae \
  -profile local \
  -resume
