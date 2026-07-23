#!/bin/bash
set -euo pipefail
# AF2 DBs at /mnt/datasets/alphafold are group=alphafold, mode 750. Re-exec under that
# group so the sbatch-submitted jobs inherit the GID
if [ "$(id -gn)" != "alphafold" ]; then exec sg alphafold -c "$0 $*"; fi

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)
DEFAULT_SLURM_ACCOUNT=$(sacctmgr --parsable2 show user -s ${USER} | tail -1 | cut -f 2 -d \|)

# Multimer complex folding: input/complex.fasta has 2 records -> chains A, B.
# jackhmmer_af2 is required for paired multimer (its rich UniProt/UniRef headers
# carry the taxonomy bin/msa_taxonomy.py turns into each engine's paired MSA;
# ColabFold headers are taxonomy-less). AF2 uses the 2021 DB snapshot for its
# native multimer pairing (see nextflow.m3.config).
nextflow run ${PIPELINE_DIR}/fold.nf \
  -c nextflow.m3.config \
  --slurm_account ${DEFAULT_SLURM_ACCOUNT} \
  --input 'input/complex.fasta' \
  --outdir results \
  --methods af2,boltz,rf3,protenix \
  --msa_method jackhmmer_af2 \
  --n_predictions 1 \
  --af2_keep_models best \
  -profile slurm,m3 -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt
