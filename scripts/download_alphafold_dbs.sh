#!/usr/bin/env bash
# Download AlphaFold genetic databases (and model params) using DeepMind's
# official scripts from https://github.com/google-deepmind/alphafold
#
# Usage:
#   scripts/download_alphafold_dbs.sh /path/to/alphafold_dbs [full_dbs|reduced_dbs]
#
# Requirements: aria2c, rsync, git, ~2.6 TB free for full_dbs (~556 GB download).
# Point fold.nf at the download directory with --af2_db_path.
set -euo pipefail

if [[ $# -lt 1 ]]; then
  cat <<'EOF'
Usage: scripts/download_alphafold_dbs.sh <DOWNLOAD_DIR> [full_dbs|reduced_dbs]

  DOWNLOAD_DIR   Destination for the AlphaFold database tree (not inside a
                 git clone of alphafold - large DBs would slow Docker builds).
  MODE           full_dbs (default) or reduced_dbs.

Requires: aria2c, rsync, git. Full DBs are ~556 GB download / ~2.62 TB unzipped.
See: https://github.com/google-deepmind/alphafold#genetic-databases
EOF
  exit 1
fi

DOWNLOAD_MODE="${2:-full_dbs}"
AF_REPO="${AF_REPO:-https://github.com/google-deepmind/alphafold.git}"
AF_REF="${AF_REF:-main}"

if [[ "${DOWNLOAD_MODE}" != "full_dbs" && "${DOWNLOAD_MODE}" != "reduced_dbs" ]]; then
  echo "Error: MODE must be full_dbs or reduced_dbs (got '${DOWNLOAD_MODE}')" >&2
  exit 1
fi

for cmd in aria2c rsync git; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Error: ${cmd} not found in PATH" >&2
    exit 1
  fi
done

mkdir -p "$1"
DOWNLOAD_DIR="$(cd "$1" && pwd)"
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/af2-db-scripts.XXXXXX")"
cleanup() { rm -rf "${WORKDIR}"; }
trap cleanup EXIT

echo "Cloning AlphaFold scripts (${AF_REPO}@${AF_REF})..."
git clone --depth 1 --branch "${AF_REF}" "${AF_REPO}" "${WORKDIR}/alphafold"

echo "Downloading AlphaFold data into ${DOWNLOAD_DIR} (${DOWNLOAD_MODE})..."
bash "${WORKDIR}/alphafold/scripts/download_all_data.sh" "${DOWNLOAD_DIR}" "${DOWNLOAD_MODE}"

echo
echo "Done. Expected layout under ${DOWNLOAD_DIR}:"
echo "  bfd/ or small_bfd/  mgnify/  params/  pdb70/  pdb_mmcif/"
echo "  uniref30/  uniref90/  (+ uniprot/ pdb_seqres/ for multimer)"
echo
echo "Use with fold.nf:"
echo "  --msa_method jackhmmer_af2 --af2_db_path ${DOWNLOAD_DIR} --af2_db_preset ${DOWNLOAD_MODE}"
echo
echo "Ensure the directory is readable by pipeline jobs (e.g. chmod -R a+rX '${DOWNLOAD_DIR}')."
echo "On Apptainer/Singularity, bind-mount ${DOWNLOAD_DIR} into the container."
