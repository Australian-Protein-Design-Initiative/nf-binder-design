#!/bin/bash
set -euo pipefail
# Re-exec under the `alphafold` group so the AF2 database under /mnt/datasets
# (group=alphafold, mode 750) is readable. Harmless if you already have access.
if [ "$(id -gn)" != "alphafold" ]; then exec sg alphafold -c "$0 $*"; fi

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/fold.nf \
  -c nextflow.m3.config \
  --input 'input/pdl1.fasta' \
  --outdir results \
  --methods af2,boltz,rf3,protenix \
  --msa_method jackhmmer_af2 \
  -profile local -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
