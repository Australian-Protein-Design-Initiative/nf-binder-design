#!/bin/bash
set -euo pipefail

# PD-L1 VHH design with Germinal on M3, with a Boltz-2 refold of the accepted designs.
#
# The BDI GPU partition (bdi/bdiq) needs a SLURM account with QOS access; the default
# below is whatever sacctmgr reports for you, override with
# SLURM_ACCOUNT=yt41 ./run-m3-bdi.sh
#
# The site configs under conf/platforms/ still use top-level `def`, which Nextflow
# >=26's default (strict) parser rejects, so ask for the v1 parser. Pinning an older
# Nextflow instead does not work any more: the workflow scripts no longer compile
# under 24.10.0, which rejects the `def` bindings in a workflow body with
# "Variable `log` already defined in the process scope".
export NXF_SYNTAX_PARSER=v1

PIPELINE_DIR=../..

SLURM_ACCOUNT="${SLURM_ACCOUNT:-$(sacctmgr --parsable2 show user -s "${USER}" | tail -1 | cut -f 2 -d \|)}"

mkdir -p results/logs

# --refold_with boltz is the held-out second opinion: Germinal has already cofolded
# every accepted design with its own AF2 pass, so the refold folds the binder against
# the target again with an independent predictor. --refold_max 2 keeps the refold short;
# drop it to refold every accepted design.
nextflow run ${PIPELINE_DIR}/main.nf \
  --method germinal \
  --slurm_account "${SLURM_ACCOUNT}" \
  --germinal_config configs/pdl1_vhh.yaml \
  --germinal_pdb_dir pdbs \
  --germinal_experiment_name pdl1_vhh \
  --germinal_n_traj 4 \
  --germinal_batch_size 1 \
  --refold_with boltz \
  --refold_max 2 \
  --outdir results \
  -profile m3_bdi \
  -resume
