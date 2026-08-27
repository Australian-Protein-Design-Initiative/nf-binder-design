// Generic per-structure confidence parser for fold.nf's RF3 and Protenix
// engines: optionally run bin/ipsae.py on the full PAE JSON + structure, then
// bin/fold/parse_fold_confidence.py flattens the summary JSON (plus ipSAE TSV)
// into a single normalized TSV row on stdout. The RF3/Protenix subworkflows fan
// this out one call per sample and collectFile the rows into
// <tool>_fold_scores.tsv. AF2 has its own FOLD_SCORE_AF2; Boltz has
// FOLD_PARSE_BOLTZ_CONFIDENCE. CPU-only, local executor.
process FOLD_PARSE_CONFIDENCE {
    tag "${meta.id} ${tool} ${model}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), val(tool), val(model), val(original_file), val(predictions_file), path(json_file), path(ipsae_pae), path(ipsae_structure), val(do_ipsae)

    output:
    stdout

    script:
    def ipsae_fmt = tool == 'rf3' ? 'rf3' : 'af3'
    """
    set -euo pipefail
    ipsae_args=()
    if [[ "${do_ipsae}" == "true" ]]; then
        python3 ${projectDir}/bin/ipsae.py --format ${ipsae_fmt} \\
            "${ipsae_pae}" "${ipsae_structure}" 10 10 \\
            || echo "ipsae.py failed for ${ipsae_structure}" >&2
        tsv=\$(ls -1 *_10_10_ipsae.tsv 2>/dev/null | head -n 1 || true)
        if [[ -n "\${tsv}" ]]; then
            ipsae_args=(--ipsae-tsv "\${tsv}")
        fi
    fi
    python3 ${projectDir}/bin/fold/parse_fold_confidence.py \\
        --tool "${tool}" \\
        --id "${meta.id}" \\
        --model "${model}" \\
        --original-file "${original_file}" \\
        --predictions-file "${predictions_file}" \\
        --json "${json_file}" \\
        "\${ipsae_args[@]}"
    """
}
