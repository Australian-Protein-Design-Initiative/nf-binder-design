#!/bin/bash

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method boltzgen \
  --config_yaml pfoa.yaml \
  --outdir results \
  --design_name pfoa \
  --protocol protein-small_molecule \
  --num_designs 2 \
  --batch_size 1 \
  --budget 1 \
  -profile local \
  -resume

