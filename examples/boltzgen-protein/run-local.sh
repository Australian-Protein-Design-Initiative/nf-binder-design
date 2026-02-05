#!/bin/bash

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method boltzgen \
  --config_yaml 1g13prot.yaml \
  --outdir results \
  --design_name 1g13prot \
  --protocol protein-anything \
  --num_designs 4 \
  --batch_size 1 \
  --budget 2 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt

