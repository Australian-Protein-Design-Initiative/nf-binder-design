// Multimer Protenix input JSON for fold.nf: one proteinChain per chain, each
// with its own unpairedMsaPath + pairedMsaPath (rendered by bin/fold/msa_taxonomy.py
// --tool protenix). Protenix pairs chains internally by species *mnemonic*, not
// numeric TaxID=. The a3m bundle arrives as a combined list of per-chain
// *.protenix_paired.a3m + *.protenix_unpaired.a3m, each group already in chain
// order (FOLD_MSA sorts by chain_index; FOLD_PULLDOWN_MSA emits target then
// binder). Split by suffix but do NOT sort by name: pulldown files are named
// after the target/binder ids, so a binder id sorting first swapped the MSAs.
process GENERATE_PROTENIX_INPUT_COMPLEX {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    tuple val(meta), path(fasta), path(a3ms)

    output:
    tuple val(meta), path(fasta), path(a3ms), path('protenix_input.json'), emit: with_json

    script:
    def files = (a3ms instanceof List) ? a3ms : [a3ms]
    def paired = files.findAll { it.name.endsWith('.protenix_paired.a3m') }
    def unpaired = files.findAll { it.name.endsWith('.protenix_unpaired.a3m') }
    def paired_arg = paired.collect { it.name }.join(' ')
    def unpaired_arg = unpaired.collect { it.name }.join(' ')
    """
    python ${projectDir}/bin/fold/make_protenix_input.py \
        --fasta ${fasta} \
        --name '${meta.id}' \
        --a3m ${unpaired_arg} \
        --paired-a3m ${paired_arg} \
        -o protenix_input.json
    """
}
