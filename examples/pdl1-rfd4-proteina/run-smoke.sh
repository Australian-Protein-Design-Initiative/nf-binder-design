#!/bin/bash

# Fastest possible check that the container and GPU work: one design, no torch.compile warm-up.

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.docker.config \
  --method rfd4_proteina \
  --rfd4_config pdl1_rfd4.yaml \
  --outdir results-smoke \
  --rfd4_n_designs=1 \
  --rfd4_inference_experiment=quickstart \
  --rfd4_compile_lengths=false \
  -profile local \
  -resume
