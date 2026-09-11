#!/bin/bash

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method boltzgen \
  --config_yaml pfoa.yaml \
  --outdir results \
  --design_name pfoa \
  --protocol protein-small_molecule \
  --num_designs 4 \
  --batch_size 2 \
  --budget 10 \
  --devices 2 \
  --num_workers 2 \
  -profile local \
  -resume

