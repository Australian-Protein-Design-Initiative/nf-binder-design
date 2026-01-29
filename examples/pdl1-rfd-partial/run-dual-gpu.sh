#!/bin/bash

# This example is designed to run on a local machine with 2 GPUs.
# We use the --gpu_devices=0,1 flag to specify the GPUs to use.
# Don't use --gpu_devices=0,1 with SLURM - the scheduler should handle this for you.

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/partial.nf  \
  -c nextflow.dual-gpu.config \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_partial_per_binder=1 \
  --rfd_batch_size=1 \
  --rfd_partial_T "5,15" \
  --pmpnn_seqs_per_struct=1 \
  --pmpnn_relax_cycles=3 \
  --rfd_filters="rg<=20" \
  --refold_af2ig_filters="pae_interaction<=10;plddt_binder>=80" \
  --af2ig_recycle=3 \
  --refold_max=100 \
  --refold_use_msa_server=true \
  --gpu_devices=0,1 \
  --enable_gpu_stats=true \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
