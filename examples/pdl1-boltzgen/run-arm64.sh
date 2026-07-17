#!/bin/bash

# ARM64 / NVIDIA GB10 (Blackwell) test run.
# Skips --do_foldseek: no arm64 foldseek container yet and no local databases directory.
# See nextflow-arm64.config for a fork()+native-thread-pool deadlock workaround
# (forces --num_workers 0 via ext.args - the CLI --num_workers 0 below is a no-op
# due to a Groovy falsy-zero bug in the module, kept only for documentation).

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow-arm64.config \
  --method boltzgen \
  --config_yaml pdl1-binder.yaml \
  --outdir results \
  --design_name pdl1-binder \
  --protocol protein-anything \
  --num_designs 4 \
  --batch_size 2 \
  --budget 2 \
  --num_workers 0 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
