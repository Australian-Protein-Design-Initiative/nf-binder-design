#!/bin/bash

# Oracle refold scoring: co-fold each design's binder sequence against the target with
# several predictors and rank on their consensus, instead of Boltz-2 alone.
#
# --refold_oracles is independent of --full_refold_with. Both are used here: Boltz-2 gives
# the RMSD-to-design and binder-monomer metrics, the oracle set gives per-tool iptm/ipSAE
# and a cross-tool consensus_z, and all of it lands in one combined_scores.tsv.
#
# Note: the oracle containers are amd64 only, so this script does not run on an arm64 host.

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.docker.config \
  --method rfd4_proteina \
  --rfd4_config pdl1_rfd4.yaml \
  --outdir results-oracles \
  --rfd4_n_designs=20 \
  --rfd4_designs_per_task=5 \
  --rfd4_filters="rg<=20" \
  --full_refold_with 'boltz' \
  --full_refold_max=10 \
  --full_refold_use_msa_server \
  --refold_oracles 'boltz,chai1,openfold3' \
  --refold_oracles_max=10 \
  --n_predictions=3 \
  --use_msa_server \
  --consensus_metric=ipsae \
  -profile local \
  -resume
