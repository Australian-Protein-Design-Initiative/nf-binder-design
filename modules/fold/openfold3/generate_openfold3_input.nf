// Monomer OpenFold3 query JSON for fold.nf: one protein chain whose MSA dir holds
// the shared a3m as colabfold_main.a3m. Seeds live in the per-task runner YAML,
// so this runs once per input rather than once per batch.
process GENERATE_OPENFOLD3_INPUT {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3m)

    output:
    tuple val(meta), path(fasta), path('msa_*'), path('openfold3_query.json'), emit: with_json

    script:
    """
    python ${projectDir}/bin/fold/make_openfold3_input.py \\
        --fasta ${fasta} \\
        --name '${meta.id}' \\
        --a3m ${a3m} \\
        -o openfold3_query.json
    """
}
