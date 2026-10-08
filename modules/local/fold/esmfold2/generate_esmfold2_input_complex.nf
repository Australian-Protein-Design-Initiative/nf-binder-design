// Multimer ESMFold2 per-chain a3ms for fold.nf. Consumes the same per-chain
// bundle as Protenix/AF3/OpenFold3 but only its *.protenix_unpaired.a3m files
// (original headers), re-rendering each with msa_taxonomy.py --tool esmfold2 so
// the hit headers carry the key=<taxid> tokens ESMFold2 pairs on.
//
// ESMFold2 has no query JSON: it takes one MSA object per chain and builds the
// paired + block-diagonal layout itself, so this stage only rewrites headers.
// Bundle order is chain order (see generate_af3_input_complex.nf), so the outputs
// are named chain_A/chain_B/... by position and the predict task sorts on that
// rather than trusting the staging order.
process GENERATE_ESMFOLD2_INPUT_COMPLEX {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3ms)

    output:
    tuple val(meta), path(fasta), path('chain_*.esmfold2.a3m'), emit: with_msa

    script:
    def files = (a3ms instanceof List) ? a3ms : [a3ms]
    def unpaired = files.findAll { it.name.endsWith('.protenix_unpaired.a3m') }
    def chain_ids = ('A'..'Z').toList()
    def renders = unpaired.withIndex().collect { f, i ->
        "python ${projectDir}/bin/fold/msa_taxonomy.py " +
            "--a3m '${f.name}' --tool esmfold2 --chain-id '${chain_ids[i]}' " +
            "--out 'chain_${chain_ids[i]}.esmfold2.a3m'"
    }.join('\n    ')
    """
    set -euo pipefail
    ${renders}
    """
}
