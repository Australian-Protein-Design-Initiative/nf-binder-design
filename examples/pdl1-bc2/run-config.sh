#!/bin/bash

# Config mode: the campaign JSON drives the run. input/PDL1.pdb is referenced by
# a relative path inside campaign.json, detected, and staged into each task.

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method bc2 \
  --bc2_config campaign.json \
  --outdir results-config \
  --bc2_n_traj 4 \
  --bc2_batch_size 1 \
  --bc2_set 'filters.i_pTM.threshold=0.7' \
  -profile local \
  -resume
