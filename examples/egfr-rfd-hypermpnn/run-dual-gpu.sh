#!/bin/bash

# This example is designed to run on a local machine with 2 GPUs.
# We use the --gpu_devices=0,1 flag to specify the GPUs to use.
# Don't use --gpu_devices=0,1 with SLURM - the scheduler should handle this for you.

# Here we use the HyperMPNN weights for ProteinMPNN 
# (run `./models/download_hypermpnn_weights.sh` first to download them)
#
# The input file 6aru.pdb contains the full structure of the EGFR/Fab complex 
# - we use the RFDiffusion --contig syntax to 'crop' to only the EGFR target domain.

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf  \
  -c nextflow.dual-gpu.config \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A310-481/0 65-120]" \
  --hotspot_res "[A412]" \
  --rfd_n_designs=4 \
  --rfd_batch_size=1 \
  --pmpnn_seqs_per_struct=2 \
  --pmpnn_relax_cycles=1 \
  --rfd_filters="rg<=20" \
  --pmpnn_weigths="$(realpath ${PIPELINE_DIR}/models/HyperMPNN/retrained_models/v48_020_epoch300_hyper.pt)" \
  --gpu_devices=0,1 \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
