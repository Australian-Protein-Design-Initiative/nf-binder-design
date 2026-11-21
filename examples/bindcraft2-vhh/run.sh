#!/bin/bash
# Humanised VHH (nanobody) binders against PD-L1, one trajectory at a time on a
# single GPU.
#
# --bc2_modality VHH switches BindCraft2 from de novo backbone generation to
# fold conditioning on its shipped VHH scaffold, redesigning the CDR loops
# against the target. The scaffold and its CDR mutation mask both live inside
# the container, so nothing extra is staged.
#
# --bc2_humanize adds human-like sequence preferences and an MHC-II anchor
# score objective, which is what you usually want for a therapeutic nanobody.

set -euo pipefail

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method bc2 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A54,A56,A66,A115" \
  --bc2_modality VHH \
  --bc2_humanize \
  --bc2_n_traj 4 \
  --bc2_batch_size 1 \
  -profile local \
  -resume
