#!/bin/bash
PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/boltzgen_filter.nf \
  --run results/boltzgen/merged \
  --budget 2 \
  --refolding_rmsd_threshold 3.0 \
  --filter_biased=false \
  --additional_filters 'ALA_fraction<0.3' 'filter_rmsd_design<2.5' \
  --metrics_override plip_hbonds_refolded=4 \
  --alpha 0.2 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt

