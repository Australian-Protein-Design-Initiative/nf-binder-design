#!/bin/bash
#
# In this example we specify additional options:
#   --pmpnn_seqs_per_struct=2
#   --pmpnn_relax_cycles=3
# to generate two sequences per backbone design with ProteinMPNN with 3 FastRelax cycles
#   --rfd_filters="rg<=20"
#   --af2ig_recycle=3
# to exclude backbone designs >20 A radius of gyration early in the pipeline
# run AlphaFold2 initial guess with 3 recycles

# We refold designs (using Boltz-2):
#  --refold_af2ig_filters="pae_interaction<=10;plddt_binder>=80"
# only for designs with AF2 initial guess scores: pae_interaction <= 10 and binder pLDDT >= 80
# We refold a maximum of 100 designs:
#  --refold_max=100
# using the public ColabFold MMSeqs2 server to generate the MSA for the target sequence:
#   --refold_use_msa_server=true
# Rather than refold using the truncated target (PDL1.pdb) we use the full target sequence:
#  --refold_target_fasta='input/full/3BIK_B.fasta'
# and use full length target template(s) to improve the target prediction:
#  --refold_target_templates='input/full/'

PIPELINE_DIR=../../

DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf  \
  --input_pdb 'input/*.pdb' \
  --outdir results \
  --contigs "[A18-132/0 65-120]" \
  --hotspot_res "A56" \
  --rfd_n_designs=4 \
  --rfd_batch_size=1 \
  --pmpnn_seqs_per_struct=2 \
  --pmpnn_relax_cycles=3 \
  --rfd_filters="rg<=20" \
  --refold_af2ig_filters="pae_interaction<=10;plddt_binder>=80" \
  --af2ig_recycle=3 \
  --refold_max=100 \
  --refold_use_msa_server=true \
  --refold_target_fasta='input/full/3BIK_B.fasta' \
  --refold_target_templates='input/full/' \
  --output_rmsd_aligned=true \
  -profile local \
  -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
