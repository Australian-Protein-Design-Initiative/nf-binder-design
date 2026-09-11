#!/bin/bash

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method boltzgen \
  --config_yaml 1g13prot.yaml \
  --outdir results \
  --design_name 1g13prot \
  --protocol protein-anything \
  --num_designs 4 \
  --batch_size 2 \
  --budget 2 \
  --devices 2 \
  --num_workers 2 \
  -profile local \
  -resume

