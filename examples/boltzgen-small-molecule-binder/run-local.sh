#!/bin/bash

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/boltzgen.nf \
  --config_yaml pfoa.yaml \
  --outdir results \
  --design_name pfoa \
  --protocol protein-small_molecule \
  --num_designs 2 \
  --batch_size 1 \
  --budget 1 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt

