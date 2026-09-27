// Multimer OpenFold3 query JSON for fold.nf. Consumes the same per-chain bundle
// as Protenix/AF3 but only its *.protenix_unpaired.a3m files (original headers):
// each becomes the chain's colabfold_main.a3m, and is re-rendered with
// msa_taxonomy.py --tool openfold3 into the uniprot_hits.a3m OpenFold3 pairs on.
// Bundle order is chain order, so it is kept rather than sorted by name (see
// generate_af3_input_complex.nf).
process GENERATE_OPENFOLD3_INPUT_COMPLEX {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3ms)

    output:
    tuple val(meta), path(fasta), path('msa_*'), path('openfold3_query.json'), emit: with_json

    script:
    def files = (a3ms instanceof List) ? a3ms : [a3ms]
    def unpaired = files.findAll { it.name.endsWith('.protenix_unpaired.a3m') }
    def unpaired_arg = unpaired.collect { it.name }.join(' ')
    def pairing_arg = unpaired.collect { it.name.replaceFirst(/\.protenix_unpaired\.a3m$/, '.of3_pairing.a3m') }.join(' ')
    """
    set -euo pipefail
    for f in ${unpaired_arg}; do
        python ${projectDir}/bin/fold/msa_taxonomy.py \\
            --a3m "\${f}" \\
            --tool openfold3 \\
            --chain-id "\${f%.protenix_unpaired.a3m}" \\
            --out "\${f%.protenix_unpaired.a3m}.of3_pairing.a3m"
    done

    python ${projectDir}/bin/fold/make_openfold3_input.py \\
        --fasta ${fasta} \\
        --name '${meta.id}' \\
        --a3m ${unpaired_arg} \\
        --pairing-a3m ${pairing_arg} \\
        -o openfold3_query.json
    """
}
