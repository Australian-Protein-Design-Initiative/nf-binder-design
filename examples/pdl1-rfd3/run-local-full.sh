#!/bin/bash

# Simple example using --contigs / --hotspot_res params (auto-generates rfd3 JSON config)

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  --method rfd3 \
  --design_name pdl1 \
  --input_pdb 'input/PDL1.pdb' \
  --outdir results \
  --contigs "[A18-132/0 50-120]" \
  --hotspot_res "A54,A56,A115,A123,A125" \
  --rfd3_n_designs=4 \
  --rfd3_filters="rg<=20" \
  --rfd3_is_non_loopy=false \
  --rfd3_hotspot_subsample=0.4 \
  --mpnn_preset 'soluble' \
  --mpnn_weights_noise '030' \
  --rf3_target_template='input/full/3BIK.pdb' \
  --rf3_template_selection='A/*/18-26,A/*/42-53,A/*/132-137,A/*/184-190,A/*/226-229' \
  --rf3_use_msa_server \
  --rf3_create_target_msa \
  --full_refold_with 'boltz' \
  --full_refold_filter_sort='pair_pae_min' \
  --full_refold_max=2 \
  --full_refold_create_target_msa \
  --full_refold_use_msa_server \
  --full_refold_target_fasta='input/full/3BIK_B.fasta' \
  --full_refold_target_templates='input/full/' \
  --output_rmsd_aligned \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt

#  --rf3_create_target_msa=true \
#  --rf3_use_msa_server=true \