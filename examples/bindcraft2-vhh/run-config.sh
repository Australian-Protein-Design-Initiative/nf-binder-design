#!/bin/bash
# The same campaign as run.sh, driven entirely by campaign.json.
#
# Two things differ. The target is BindCraft2's shipped "hPDL1" preset rather
# than a local PDB, so no structure is staged at all - the pipeline logs "no
# files to stage (shipped presets only)". And the modality and the humanize
# property are set in the JSON instead of on the command line, which is
# equivalent to --bc2_modality VHH --bc2_humanize.
#
# Run `bindcraft design --list-targets` in the container for the other shipped
# targets (hPD1, hIL2R, hIL7RA, mPDL1, dynorphin_a).

set -euo pipefail

PIPELINE_DIR=../../

nextflow run ${PIPELINE_DIR}/main.nf \
  --method bc2 \
  --bc2_config campaign.json \
  --outdir results-config \
  --bc2_n_traj 4 \
  --bc2_batch_size 1 \
  -profile local \
  -resume
