#!/bin/bash

# Full RFD4-Proteina path on an arm64 host (eg DGX Spark GB10): generate, filter, refold with
# Boltz-2, and combine. conf/platforms/arm64.config swaps the amd64-only images for their arm64
# builds, without which the post-generation steps fail with `exec format error`.
#
# Kept small on purpose. Raise --rfd4_n_designs and drop the quickstart settings for a real run.

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c ${PIPELINE_DIR}/conf/platforms/arm64.config \
  -c nextflow.docker.config \
  --method rfd4_proteina \
  --rfd4_config pdl1_rfd4.yaml \
  --outdir results-arm64 \
  --rfd4_n_designs=2 \
  --rfd4_designs_per_task=2 \
  --rfd4_num_workers=0 \
  --rfd4_inference_experiment=quickstart \
  --rfd4_compile_lengths=false \
  --rfd4_filters="rg<=25" \
  --full_refold_with 'boltz' \
  --full_refold_max=2 \
  --full_refold_use_msa_server \
  -profile local \
  -resume
