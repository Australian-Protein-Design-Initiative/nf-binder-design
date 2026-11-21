#!/bin/bash

# Full RFD4-Proteina path on an arm64 host (eg DGX Spark GB10): generate, filter, refold, and
# combine. conf/platforms/arm64.config swaps the amd64-only images for their arm64 builds, without
# which the post-generation steps fail with `exec format error`.
#
# boltz, esmfold2, esmfold2_fast and rf3 have arm64 images, so any of them can go in
# --refold_with here. RF3's own flags are --refold_rf3_* in this workflow, since the
# plain --rf3_* names belong to the design-time RF3 step in --method rfd3.
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
  --refold_with 'boltz,esmfold2_fast' \
  --refold_max=2 \
  --use_msa_server \
  -profile local \
  -resume
