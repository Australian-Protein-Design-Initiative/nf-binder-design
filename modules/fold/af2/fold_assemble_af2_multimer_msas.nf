// Assemble AF2 multimer precomputed-MSA tree + features.pkl for one
// target+binder pair. Target MSA: jackhmmer dir and/or ColabFold a3m for
// chain A; binder B is always query-only. Pass assets/dummy_files/empty
// for unused slots. Uses the AF2 container so features.pkl can be built
// with alphafold.data (predict loads that pickle, not the raw msas/).
process FOLD_ASSEMBLE_AF2_MULTIMER_MSAS {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/ghcr.io-australian-protein-design-initiative-containers-alphafold2-2.3.2-custom.img'

    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/af2/msas",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString()
            return rel.startsWith("${meta.id}/") ? rel : null
        }
    )

    input:
    tuple val(meta), path(pair_fasta), path(target_msa_dir), path(target_a3m)

    output:
    tuple val(meta), path(pair_fasta), path("${meta.id}"), emit: msas

    script:
    """
    python ${projectDir}/bin/fold/assemble_af2_multimer_msas.py \\
        --pair-fasta ${pair_fasta} \\
        --pair-id '${meta.id}' \\
        --target-msa-dir ${target_msa_dir} \\
        --target-a3m ${target_a3m} \\
        -o .

    python ${projectDir}/bin/fold/af2_multimer_features_from_msas.py \\
        --fasta ${pair_fasta} \\
        --msas-dir '${meta.id}'
    """
}
