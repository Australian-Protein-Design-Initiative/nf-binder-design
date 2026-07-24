// Generic per-structure confidence parser for fold.nf's RF3 and Protenix
// engines: bin/parse_fold_confidence.py flattens one summary JSON into a single
// normalized TSV row (header + row) on stdout. The RF3/Protenix subworkflows fan
// this out one call per sample and collectFile the rows into
// <tool>_fold_scores.tsv. AF2 has its own FOLD_SCORE_AF2 (adds ipSAE); Boltz has
// FOLD_PARSE_BOLTZ_CONFIDENCE. CPU-only, local executor.
process FOLD_PARSE_CONFIDENCE {
    tag "${meta.id} ${tool} ${model}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), val(tool), val(model), val(original_file), val(predictions_file), path(json_file)

    output:
    stdout

    script:
    """
    python3 ${projectDir}/bin/parse_fold_confidence.py \
        --tool "${tool}" \
        --id "${meta.id}" \
        --model "${model}" \
        --original-file "${original_file}" \
        --predictions-file "${predictions_file}" \
        --json "${json_file}"
    """
}
