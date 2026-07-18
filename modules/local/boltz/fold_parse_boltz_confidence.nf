// fold.nf-local confidence parser: reuses bin/parse_boltz_confidence.py (the
// same script boltz_pulldown.nf's PARSE_BOLTZ_CONFIDENCE_JSON uses) without
// the target/binder metadata columns that script's binder-design callers add
// - fold.nf's meta has no target/binder split (it's a single fold target).
process FOLD_PARSE_BOLTZ_CONFIDENCE {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(json_file), path(ipsae_tsv)

    output:
    stdout

    script:
    """
    python3 ${projectDir}/bin/parse_boltz_confidence.py \
        --json "${json_file}" \
        --id "${meta.id}" \
        --merge-ipsae "${ipsae_tsv}"
    """
}
