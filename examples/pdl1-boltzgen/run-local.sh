#!/bin/bash

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method boltzgen \
  --config_yaml pdl1-binder.yaml \
  --outdir results \
  --design_name pdl1-binder \
  --protocol protein-anything \
  --num_designs 4 \
  --batch_size 2 \
  --budget 2 \
  --do_foldseek \
  --foldseek_database CATH50 \
  --foldseek_databases_path "$(pwd)/../databases/foldseek" \
  -profile local \
  -resume
