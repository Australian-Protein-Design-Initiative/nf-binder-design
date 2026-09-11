#!/bin/bash
set -euo pipefail

# Pin Nextflow 24.10.0: site configs under conf/platforms/ still use top-level
# `def`, which Nextflow >=26's default (strict) parser rejects. Use
# NXF_SYNTAX_PARSER=v1 with Nextflow 26, or pin <26 as here.
export NXF_VER=24.10.0

PIPELINE_DIR=../..

DEFAULT_SLURM_ACCOUNT=$(sacctmgr --parsable2 show user -s ${USER} | tail -1 | cut -f 2 -d \|)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method bindcraft \
  --slurm_account $DEFAULT_SLURM_ACCOUNT \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-120" \
  --bindcraft_n_traj 4 \
  --bindcraft_batch_size 1 \
  --bindcraft_advanced_settings_preset "default_4stage_multimer" \
  -profile slurm,m3 \
  -resume

# Alternatively, instead of --target_chains you can specify RFDiffusion-style 
# contigs defining the regions to use, eg:
#   --contigs "[A18-132/0]"

# --do_foldseek \
