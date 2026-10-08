#!/bin/bash
#
# Download the AlphaFold3 model parameters (af3.bin.zst) for --methods af3.
#
# The weights are NOT distributed with nf-binder-design or its containers - they
# are subject to Google DeepMind's AlphaFold3 Model Parameters Terms of Use and
# must be obtained by each user/organisation.
#
# Usage:
#   ./download_af3_weights.sh [-o DIR] [--accept-terms] [--decompress]
#
#   -o, --output DIR   Target directory (default: models/alphafold3 next to this
#                      script, which is the pipeline's default --af3_model_dir)
#   --accept-terms     Non-interactive acceptance of the terms of use
#   --decompress       Store af3.bin instead of af3.bin.zst (AF3 reads either;
#                      the compressed form is smaller and loads fine)

set -euo pipefail

AF3_WEIGHTS_URL="https://storage.googleapis.com/alphafold3/af3.bin.zst"
TERMS_URL="https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md"
POLICY_URL="https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_PROHIBITED_USE_POLICY.md"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/alphafold3"
ACCEPT_TERMS=false
DECOMPRESS=false

usage() {
    sed -n '2,17p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output)
            OUT_DIR="$2"
            shift 2
            ;;
        --accept-terms)
            ACCEPT_TERMS=true
            shift
            ;;
        --decompress)
            DECOMPRESS=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

cat >&2 <<EOF

################################################################################
#                                                                              #
#               ALPHAFOLD3 MODEL PARAMETERS - TERMS OF USE                     #
#                                                                              #
################################################################################

The AlphaFold3 model parameters are made available by Google DeepMind under the
AlphaFold3 Model Parameters Terms of Use and Prohibited Use Policy:

  ${TERMS_URL}
  ${POLICY_URL}

Read these in full before continuing. In summary (this is NOT a substitute for
the terms themselves):

  * NON-COMMERCIAL USE ONLY, by or on behalf of a non-commercial organisation
    (universities, non-profit research institutes, government bodies, etc).
    No commercial activities, including research done on behalf of a
    commercial organisation.
  * DO NOT SHARE OR REDISTRIBUTE the parameters outside your organisation
    (do not bake them into public containers, public buckets or git repos).
  * Do not use AlphaFold3 outputs to train other structure prediction models.
  * Outputs are subject to the AlphaFold3 Output Terms of Use.

Weights will be downloaded from:
  ${AF3_WEIGHTS_URL}
into:
  ${OUT_DIR}

EOF

if [[ "${ACCEPT_TERMS}" != "true" ]]; then
    if [[ ! -t 0 ]]; then
        echo "ERROR: not running interactively - re-run with --accept-terms once you have read and agree to the terms." >&2
        exit 1
    fi
    read -r -p "Have you read and do you agree to the AlphaFold3 Model Parameters Terms of Use? Type 'yes' to continue: " answer
    if [[ "${answer}" != "yes" ]]; then
        echo "Terms not accepted - aborting." >&2
        exit 1
    fi
fi

mkdir -p "${OUT_DIR}"

# AF3's --model_dir loader refuses to start if more than one model file matches,
# so never leave a second (different) weights file alongside ours.
existing=$(find "${OUT_DIR}" -maxdepth 1 -type f \( -name '*.bin' -o -name '*.bin.zst' \) \
    ! -name 'af3.bin.zst' ! -name 'af3.bin' | head -n 1)
if [[ -n "${existing}" ]]; then
    echo "ERROR: ${OUT_DIR} already contains another model file (${existing}). AlphaFold3 requires exactly one model file per --model_dir." >&2
    exit 1
fi

if [[ -s "${OUT_DIR}/af3.bin" || -s "${OUT_DIR}/af3.bin.zst" ]]; then
    echo "AlphaFold3 weights already present in ${OUT_DIR} - nothing to do." >&2
else
    tmp="${OUT_DIR}/af3.bin.zst.part"
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 -C - -o "${tmp}" "${AF3_WEIGHTS_URL}"
    elif command -v wget >/dev/null 2>&1; then
        wget -c -O "${tmp}" "${AF3_WEIGHTS_URL}"
    else
        echo "ERROR: need curl or wget to download the weights." >&2
        exit 1
    fi
    mv "${tmp}" "${OUT_DIR}/af3.bin.zst"

    if [[ "${DECOMPRESS}" == "true" ]]; then
        if ! command -v zstd >/dev/null 2>&1; then
            echo "ERROR: --decompress needs zstd (or omit it - AlphaFold3 reads af3.bin.zst directly)." >&2
            exit 1
        fi
        zstd -d --rm "${OUT_DIR}/af3.bin.zst" -o "${OUT_DIR}/af3.bin"
    fi
fi

chmod -R go-rwx "${OUT_DIR}" 2>/dev/null || true

echo "" >&2
echo "AlphaFold3 weights are in: ${OUT_DIR}" >&2
ls -lh "${OUT_DIR}" >&2
echo "" >&2
if [[ "${OUT_DIR}" == "${SCRIPT_DIR}/alphafold3" ]]; then
    echo "This is the default location - no extra pipeline flags needed for --methods af3." >&2
else
    echo "Run the pipeline with:  --af3_model_dir ${OUT_DIR}" >&2
fi
