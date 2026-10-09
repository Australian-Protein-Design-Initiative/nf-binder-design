#!/bin/bash
#
# As run.sh, but with AlphaFold3 instead of Protenix for structure prediction
# (configs/pdl1_vhh_af3.yaml sets structure_model: "af3").
#
# AlphaFold3's weights and databases are not in the container and must be bound
# in -- see docs/workflows/germinal.md. Point these at your own copies:
#
#   AF3_WEIGHTS    directory holding exactly one af3.bin.zst
#                  (../../models/download_af3_weights.sh will fetch it)
#   AF3_DATABASES  the AlphaFold3 public databases (~630 GB)
#                  on M3, use /mnt/datasets/alphafold3/3.0.0
#
# e.g. AF3_WEIGHTS=/data/af3_weights AF3_DATABASES=/data/af3_databases ./run-af3.sh

set -euo pipefail

PIPELINE_DIR=../../

: "${AF3_WEIGHTS:?set AF3_WEIGHTS to the directory holding af3.bin.zst}"
: "${AF3_DATABASES:?set AF3_DATABASES to the AlphaFold3 public databases directory}"
export AF3_WEIGHTS AF3_DATABASES

mkdir -p results/logs

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.af3.config \
  --method germinal \
  --germinal_config configs/pdl1_vhh_af3.yaml \
  --germinal_pdb_dir pdbs \
  --germinal_experiment_name pdl1_vhh_af3 \
  --germinal_n_traj 2 \
  --germinal_batch_size 1 \
  --outdir results \
  -profile local \
  -resume
