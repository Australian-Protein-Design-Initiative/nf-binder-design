#!/bin/bash

# Simple example using --contigs / --hotspot_res params (auto-generates the condition spec)

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.docker.config \
  --method rfd4_proteina \
  --design_name pdl1 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --contigs "50-120, /0, A18-132" \
  --hotspot_res "A54,A56,A115,A123,A125" \
  --rfd4_n_designs=4 \
  --rfd4_designs_per_task=2 \
  -profile local \
  -resume
