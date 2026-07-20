process GENERATE_RF3_FOLD_INPUT {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3m)

    output:
    tuple val(meta), path(fasta), path(a3m), path('rf3_fold.json'), emit: with_json

    script:
    """
    python ${projectDir}/bin/make_rf3_fold_spec.py \
        --fasta ${fasta} \
        --name '${meta.id}' \
        --a3m ${a3m} \
        -o rf3_fold.json
    """
}
