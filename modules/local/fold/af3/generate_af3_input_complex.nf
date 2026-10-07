// Multimer AlphaFold3 input JSON for fold.nf. Consumes the same per-chain bundle
// as Protenix (*.protenix_paired.a3m + *.protenix_unpaired.a3m), but only uses
// the unpaired files: the Protenix paired render's UniRef100_ACC_SPECIES headers
// do not match AF3's species regex, so each unpaired a3m (original headers) is
// re-rendered with msa_taxonomy.py --tool af3 into tr|ACC|ACC_SPECIES headers.
// Bundle order is chain order (FOLD_MSA sorts by chain_index; FOLD_PULLDOWN_MSA
// emits target then binder), so it is kept rather than sorting by name, which
// would swap chains in a pulldown whenever the binder id sorts first.
process GENERATE_AF3_INPUT_COMPLEX {
    tag "${meta.id}${meta.fold_batch ? " batch${meta.fold_batch}" : ''}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3ms)
    path templates

    output:
    tuple val(meta), path(fasta), path('chain_*_{unpaired,paired}.a3m'), path('af3_input.json'), emit: with_json

    script:
    def files = (a3ms instanceof List) ? a3ms : [a3ms]
    def unpaired = files.findAll { it.name.endsWith('.protenix_unpaired.a3m') }
    def unpaired_arg = unpaired.collect { it.name }.join(' ')
    def template_chains = (meta.template_chains ?: []) as List
    def templates_arg = params.templates \
        ? "--templates-dir ${templates}" + (template_chains ? " --template-chains ${template_chains.join(' ')}" : '') \
        : ''
    def paired_arg = unpaired.collect { it.name.replaceFirst(/\.protenix_unpaired\.a3m$/, '.af3_paired.a3m') }.join(' ')
    """
    set -euo pipefail
    for f in ${unpaired_arg}; do
        python ${projectDir}/bin/fold/msa_taxonomy.py \\
            --a3m "\${f}" \\
            --tool af3 \\
            --chain-id "\${f%.protenix_unpaired.a3m}" \\
            --out "\${f%.protenix_unpaired.a3m}.af3_paired.a3m"
    done

    python ${projectDir}/bin/fold/make_af3_input.py \\
        --fasta ${fasta} \\
        --name '${meta.id}' \\
        --a3m ${unpaired_arg} \\
        --paired-a3m ${paired_arg} \\
        --seed ${meta.af3_seed} \\
        ${templates_arg} \\
        -o af3_input.json
    """
}
