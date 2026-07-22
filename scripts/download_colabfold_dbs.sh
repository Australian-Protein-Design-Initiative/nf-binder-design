#!/usr/bin/env bash
# Download and index ColabFold MMseqs2 databases for local MSA search.
# Uses ColabFold's setup_databases.sh (https://colabfold.mmseqs.com/ /
# https://github.com/sokrypton/ColabFold).
#
# Usage:
#   scripts/download_colabfold_dbs.sh /path/to/colabfold_dbs
#
# Optional environment variables:
#   UNIREF30DB              UniRef30 archive stem (default: uniref30_2302)
#   CFDB                    Env DB archive stem (default: colabfold_envdb_202108)
#   MMSEQS_NO_INDEX=1       Skip mmseqs createindex (faster setup; slower search)
#   DOWNLOADS_ONLY=1        Download archives only (no unpack / index)
#   FAST_PREBUILT_DATABASES Prebuilt .db.tar.gz (default: 1; ColabFold default)
#   GPU=1                   Build GPU-capable indexes (needs GPU MMseqs2)
#   SKIP_TEMPLATES=1        Skip PDB mmCIF / Foldseek template downloads
#
# Requirements: mmseqs in PATH, plus aria2c or curl or wget.
# Indexing the full ColabFold env DB typically wants hundreds of GB of RAM
# (ColabFold docs cite ~768 GB with preloaded indexes for single-query search).
set -euo pipefail

if [[ $# -lt 1 ]]; then
  cat <<'EOF'
Usage: scripts/download_colabfold_dbs.sh <DOWNLOAD_DIR>

  DOWNLOAD_DIR   Destination for ColabFold / MMseqs2 databases.

Requires: mmseqs, and one of aria2c / curl / wget.
Optional: SKIP_TEMPLATES=1 MMSEQS_NO_INDEX=1 DOWNLOADS_ONLY=1 GPU=1

See: https://colabfold.mmseqs.com/
     https://github.com/sokrypton/ColabFold/blob/main/setup_databases.sh
EOF
  exit 1
fi

UNIREF30DB="${UNIREF30DB:-uniref30_2302}"
CFDB="${CFDB:-colabfold_envdb_202108}"
SETUP_URL="${SETUP_URL:-https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh}"

if ! command -v mmseqs >/dev/null 2>&1; then
  echo "Error: mmseqs not found in PATH (install MMseqs2 first)" >&2
  exit 1
fi

mkdir -p "$1"
DOWNLOAD_DIR="$(cd "$1" && pwd)"
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/colabfold-db-scripts.XXXXXX")"
cleanup() { rm -rf "${WORKDIR}"; }
trap cleanup EXIT

echo "Fetching ColabFold setup_databases.sh..."
if command -v curl >/dev/null 2>&1; then
  curl -fsSL "${SETUP_URL}" -o "${WORKDIR}/setup_databases.sh"
elif command -v wget >/dev/null 2>&1; then
  wget -q -O "${WORKDIR}/setup_databases.sh" "${SETUP_URL}"
else
  echo "Error: need curl or wget to fetch setup_databases.sh" >&2
  exit 1
fi
chmod +x "${WORKDIR}/setup_databases.sh"

echo "Running setup_databases.sh in ${DOWNLOAD_DIR}..."
echo "  UNIREF30DB=${UNIREF30DB}  CFDB=${CFDB}"
echo "  FAST_PREBUILT_DATABASES=${FAST_PREBUILT_DATABASES:-1}"
echo "  MMSEQS_NO_INDEX=${MMSEQS_NO_INDEX:-}  SKIP_TEMPLATES=${SKIP_TEMPLATES:-}"
(
  cd "${DOWNLOAD_DIR}"
  export UNIREF30DB CFDB
  # Prefer prebuilt expandable-profile DBs (default in current ColabFold script).
  export FAST_PREBUILT_DATABASES="${FAST_PREBUILT_DATABASES:-1}"
  # Optional passthroughs (empty = unset behaviour from upstream script).
  [[ -n "${MMSEQS_NO_INDEX:-}" ]] && export MMSEQS_NO_INDEX
  [[ -n "${DOWNLOADS_ONLY:-}" ]] && export DOWNLOADS_ONLY
  [[ -n "${GPU:-}" ]] && export GPU
  [[ -n "${SKIP_TEMPLATES:-}" ]] && export SKIP_TEMPLATES
  bash "${WORKDIR}/setup_databases.sh" "${DOWNLOAD_DIR}"
)

# Arrange prefixes so fold.nf / boltz_pulldown --uniref30 and --colabfold_envdb
# can point at dedicated directories (MMSEQS_COLABFOLDSEARCH globs uniref30_* and
# colabfold_envdb* inside each path).
UNIREF_DIR="${DOWNLOAD_DIR}/uniref30"
ENVDB_DIR="${DOWNLOAD_DIR}/colabfold_envdb"
mkdir -p "${UNIREF_DIR}" "${ENVDB_DIR}"

shopt -s nullglob
moved_any=0
for f in "${DOWNLOAD_DIR}/${UNIREF30DB}"*; do
  base="$(basename "${f}")"
  if [[ ! -e "${UNIREF_DIR}/${base}" ]]; then
    mv "${f}" "${UNIREF_DIR}/"
    moved_any=1
  fi
done
for f in "${DOWNLOAD_DIR}/${CFDB}"*; do
  base="$(basename "${f}")"
  if [[ ! -e "${ENVDB_DIR}/${base}" ]]; then
    mv "${f}" "${ENVDB_DIR}/"
    moved_any=1
  fi
done
shopt -u nullglob

# setup_databases names DBs as ${UNIREF30DB}_db; pipeline globs uniref30_*.
# Ensure the expected prefix exists (uniref30_2302_db already matches).
if [[ ! -e "${UNIREF_DIR}/${UNIREF30DB}_db" && ! -e "${UNIREF_DIR}/${UNIREF30DB}" ]]; then
  echo "Warning: expected ${UNIREF30DB}_db under ${UNIREF_DIR} not found" >&2
fi
if [[ ! -e "${ENVDB_DIR}/${CFDB}_db" && ! -e "${ENVDB_DIR}/${CFDB}" ]]; then
  echo "Warning: expected ${CFDB}_db under ${ENVDB_DIR} not found" >&2
fi

echo
if [[ "${moved_any}" -eq 1 ]]; then
  echo "Organised databases into:"
  echo "  ${UNIREF_DIR}"
  echo "  ${ENVDB_DIR}"
fi
echo
echo "Use with fold.nf (local ColabFold MSA):"
echo "  --msa_method mmseqs2_colabfold \\"
echo "  --uniref30 ${UNIREF_DIR} \\"
echo "  --colabfold_envdb ${ENVDB_DIR}"
echo
echo "Or use the remote ColabFold API instead (no local DBs):"
echo "  --msa_method mmseqs2_colabfold --use_remote_server true"
echo
echo "On Apptainer/Singularity, bind-mount ${DOWNLOAD_DIR} into the container."
