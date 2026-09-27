#!/bin/bash
set -euo pipefail
# AlphaFold3 + Protenix on the PD-L1 monomer (shared jackhmmer MSA).

# jackhmmer_af2 MSAs read the AF2 DBs at /mnt/datasets/alphafold (group=alphafold,
# mode 750). Re-exec under that group so submitted jobs inherit the GID.
if [ "$(id -gn)" != "alphafold" ]; then exec sg alphafold -c "$0 $*"; fi

# AlphaFold3 weights are not bundled - fetch them once with
# ../../models/download_af3_weights.sh (default location), or set AF3_MODEL_DIR.
AF3_MODEL_DIR=${AF3_MODEL_DIR:-$(cd ../.. && pwd)/models/alphafold3}

PIPELINE_DIR=../..
DATESTAMP=$(date +%Y%m%d_%H%M%S)

nextflow run ${PIPELINE_DIR}/main.nf \
  -c nextflow.m3.config \
  --method fold \
  --input 'input/pdl1.fasta' \
  --outdir results-af3 \
  --methods af3,protenix \
  --af3_model_dir "${AF3_MODEL_DIR}" \
  --msa_method jackhmmer_af2 \
  --n_predictions 1 \
  -profile local -resume \
  -with-report "results-af3/logs/report_${DATESTAMP}.html" \
  -with-trace "results-af3/logs/trace_${DATESTAMP}.txt" \
  "$@"
