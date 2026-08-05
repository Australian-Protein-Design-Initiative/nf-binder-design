#!/bin/bash

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method rfd \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=4 \
  -profile local \
  -resume
