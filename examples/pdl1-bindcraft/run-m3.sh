#!/bin/bash

PIPELINE_DIR=../../software/nf-binder-design

DATESTAMP=$(date +%Y%m%d_%H%M%S)

DEFAULT_SLURM_ACCOUNT=$(sacctmgr --parsable2 show user -s ${USER} | tail -1 | cut -f 2 -d \|)

nextflow run ${PIPELINE_DIR}/bindcraft.nf \
  --slurm_account $DEFAULT_SLURM_ACCOUNT \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --target_chains "A" \
  --hotspot_res "A56" \
  --binder_length_range "55-120" \
  --bindcraft_n_traj 4 \
  --bindcraft_batch_size 1 \
  --bindcraft_advanced_settings_preset "default_4stage_multimer" \
  # --do_foldseek \
  -profile slurm,m3 \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt

# Alternatively, instead of --target_chains you can specify RFDiffusion-style 
# contigs defining the regions to use, eg:
#   --contigs "[A18-132/0]"
