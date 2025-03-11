#!/bin/bash

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <results_dir>" >&2
    echo
    echo "This is a helper script to aggregate scores from the ./results directory"
    echo "It's almost equivalent to the final COMBINE_SCORES task of the pipeline,"
    echo "but can be run mid-workflow run to get a quick overview of the scores so far."
    echo
    exit 1
fi

# Get the directory where this script lives
script_dir="$(dirname "$0")"

results_dir="$1"
pdb_dir="${results_dir}/af2_initial_guess/pdbs"
# scores_dir="${results_dir}/af2_initial_guess"  # DEPRECATED old path
scores_dir="${results_dir}/af2_initial_guess/scores"

# Ensure directories exist
if [ ! -d "$pdb_dir" ] || [ ! -d "$scores_dir" ]; then
    echo "Error: Required directories not found in results path" >&2
    echo "Expected: ${pdb_dir} and ${scores_dir}" >&2
    exit 1
fi

# Create temp files that will be cleaned up on exit
tmp_af2_scores=$(mktemp)
tmp_shape_scores=$(mktemp)
tmp_combined_scores=$(mktemp)

cleanup() {
    rm -f "$tmp_af2_scores" "$tmp_shape_scores" "$tmp_combined_scores"
}
trap cleanup EXIT

# Run AF2 score aggregation
"${script_dir}/bin/af2_combine_scores.py" "$scores_dir" --output "$tmp_af2_scores"

# Calculate shape scores
"${script_dir}/bin/calculate_shape_scores.py" --chain A ${pdb_dir}/*.pdb >"$tmp_shape_scores"

# Merge scores
"${script_dir}/bin/merge_scores.py" "$tmp_af2_scores" "$tmp_shape_scores" | tee "$tmp_combined_scores"

# Generate FASTA output to stdout
# "${script_dir}/bin/pdb_to_fasta.py" \
# --scores-table "$tmp_combined_scores" \
# --scores pae_interaction,rg,length \
# --chain A \
# "${pdb_dir}"/*.pdb