#!/bin/bash

# Example using a user-supplied RFD4-Proteina condition spec, which gives access to the full
# conditioning vocabulary (burial, hydrogen bonds, pi-stacking, CATH labels, and so on).

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.docker.config \
  --method rfd4_proteina \
  --rfd4_config pdl1_rfd4.yaml \
  --outdir results \
  --rfd4_n_designs=4 \
  --rfd4_designs_per_task=2 \
  -profile local \
  -resume
