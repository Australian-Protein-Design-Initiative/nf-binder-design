process GENERATE_RF3_FOLD_INPUT {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3m)
    path templates

    output:
    tuple val(meta), path(fasta), path(a3m), path('rf3_fold.json'), path('rf3_templates'), emit: with_json

    script:
    def template_chains = (meta.template_chains ?: []) as List
    def templates_arg = params.templates \
        ? "--templates-dir ${templates}" + (template_chains ? " --template-chains ${template_chains.join(' ')}" : '') \
        : ''
    """
    python ${projectDir}/bin/fold/make_rf3_fold_spec.py \
        --fasta ${fasta} \
        --name '${meta.id}' \
        --a3m ${a3m} \
        ${templates_arg} \
        -o rf3_fold.json
    """
}
