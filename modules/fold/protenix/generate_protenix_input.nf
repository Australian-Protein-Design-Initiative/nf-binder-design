process GENERATE_PROTENIX_INPUT {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3m)

    output:
    tuple val(meta), path(fasta), path(a3m), path('protenix_input.json'), emit: with_json

    script:
    """
    python ${projectDir}/bin/make_protenix_input.py \
        --fasta ${fasta} \
        --name '${meta.id}' \
        --a3m ${a3m} \
        -o protenix_input.json
    """
}
