#!/bin/bash

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p results/logs

nextflow run ${PIPELINE_DIR}/main.nf \
  --method germinal \
  --germinal_config configs/pdl1_vhh.yaml \
  --germinal_pdb_dir pdbs \
  --germinal_experiment_name pdl1_vhh \
  --germinal_n_traj 2 \
  --germinal_batch_size 1 \
  --outdir results \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
