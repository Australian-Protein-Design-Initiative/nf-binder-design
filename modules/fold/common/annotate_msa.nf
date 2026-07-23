// Render one per-chain unpaired a3m into each engine's native paired-MSA format
// via the shared bin/msa_taxonomy.py (the single source of taxonomy truth; see
// plans/fold-nf-multimer-paired-msa.md §4). Runs once per chain; the grouped
// per-complex bundle is assembled downstream in FOLD_MSA.
//
// Renders all three tools' files unconditionally (cheap) so a mixed --methods
// run needs only one ANNOTATE_MSA task per chain. meta.id is the per-chain id
// (<complex>.chain_<CHAIN>), so every output name is unique across chains.
process ANNOTATE_MSA {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    // Publish the rendered per-chain paired/unpaired MSAs so pairing depth is
    // inspectable (the msa_taxonomy.py stderr log also reports it per chain).
    publishDir "${params.outdir}/fold/msa/paired", pattern: "*.{a3m,csv}", mode: 'copy'

    input:
    tuple val(meta), path(fasta), path(a3m)

    output:
    tuple val(meta), path(rf3_a3m), path(px_paired), path(px_unpaired), path(boltz_csv), emit: rendered

    script:
    def stem = meta.id
    def chain = meta.chain_id ?: '?'
    rf3_a3m = "${stem}.rf3.a3m"
    px_paired = "${stem}.protenix_paired.a3m"
    px_unpaired = "${stem}.protenix_unpaired.a3m"
    boltz_csv = "${stem}.boltz.csv"
    """
    python ${projectDir}/bin/msa_taxonomy.py --a3m ${a3m} --tool rf3 \
        --out ${rf3_a3m} --chain-id '${chain}'
    python ${projectDir}/bin/msa_taxonomy.py --a3m ${a3m} --tool protenix \
        --paired-out ${px_paired} --unpaired-out ${px_unpaired} --chain-id '${chain}'
    python ${projectDir}/bin/msa_taxonomy.py --a3m ${a3m} --tool boltz \
        --out ${boltz_csv} --chain-id '${chain}'
    """
}
