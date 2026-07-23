// Split a multi-record complex FASTA into one single-record FASTA per chain,
// for fold.nf's multimer per-chain MSA search (plans/fold-nf-multimer-paired-msa.md
// §4.1). Each output is named <id>.chain_<CHAIN>.fasta so its stem is unique
// (the AF2 jackhmmer MSA stage keys its output dir on the FASTA stem).
process SPLIT_COMPLEX_FASTA {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.chain_*.fasta"), emit: chains

    script:
    """
    python ${projectDir}/bin/split_complex_fasta.py \
        --fasta ${fasta} \
        --name '${meta.id}'
    """
}
