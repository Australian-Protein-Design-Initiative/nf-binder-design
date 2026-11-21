#!/bin/bash

# A humanised VHH (nanobody) campaign against the same target. --bc2_modality
# supplies the scaffold and its CDR mutation mask from inside the container, so
# nothing extra needs staging.

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method bc2 \
  --bc2_config campaign.json \
  --outdir results-vhh \
  --bc2_modality VHH \
  --bc2_humanize \
  --bc2_n_traj 4 \
  --bc2_batch_size 1 \
  -profile local \
  -resume
