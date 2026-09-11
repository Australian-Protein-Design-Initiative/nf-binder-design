#!/bin/bash
set -euo pipefail
# AF2 DBs at /mnt/datasets/alphafold are group=alphafold, mode 750.
if [ "$(id -gn)" != "alphafold" ]; then exec sg alphafold -c "$0 $*"; fi

# Pin Nextflow 24.10.0: conf/platforms/m3.config uses `def random_choice(...)`,
# which Nextflow >=26 fails to parse ("Unexpected input: '('").
export NXF_VER=24.10.0

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)
DEFAULT_SLURM_ACCOUNT=$(sacctmgr --parsable2 show user -s ${USER} | tail -1 | cut -f 2 -d \|)

# Mosaic Multispecifics binders x PD-L1 + IL-7Ra.
# ColabFold remote MSA for targets; all fold engines. AF2 uses query-only target
# MSA under mmseqs2_colabfold (see README); Boltz/RF3/Protenix get the ColabFold a3ms.
nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.m3.config \
  --method fold_pulldown \
  --slurm_account ${DEFAULT_SLURM_ACCOUNT} \
  --targets input/targets.fasta \
  --binders input/binders.fasta \
  --outdir results \
  --methods af2,af2_mono,boltz,rf3,protenix \
  --msa_method mmseqs2_colabfold \
  --use_remote_server true \
  --create_target_msa true \
  --create_binder_msa false \
  --n_predictions 5 \
  --af2_keep_models all \
  -profile slurm,m3 -resume \
  -with-report results/logs/report_${DATESTAMP}.html \
  -with-trace results/logs/trace_${DATESTAMP}.txt \
  "$@"
