// Multimer RF3 input JSON for fold.nf: one {seq, chain_id, msa_path} component
// per chain, each pointing at that chain's TaxID=-annotated a3m (rendered by
// bin/msa_taxonomy.py --tool rf3). RF3/atomworks pairs chains internally by
// matching TaxID=<n>. a3ms arrive as a chain-ordered list (FOLD_MSA); each is
// staged alongside the JSON so RF3_FOLD finds them by basename.
process GENERATE_RF3_FOLD_INPUT_COMPLEX {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3ms)

    output:
    tuple val(meta), path(fasta), path(a3ms), path('rf3_fold.json'), emit: with_json

    script:
    def files = (a3ms instanceof List) ? a3ms : [a3ms]
    def a3m_arg = files.collect { it.name }.join(' ')
    """
    python ${projectDir}/bin/make_rf3_fold_spec.py \
        --fasta ${fasta} \
        --name '${meta.id}' \
        --a3m ${a3m_arg} \
        -o rf3_fold.json
    """
}
