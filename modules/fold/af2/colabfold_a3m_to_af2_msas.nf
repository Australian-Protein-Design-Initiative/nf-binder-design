process COLABFOLD_A3M_TO_AF2_MSAS {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/ghcr.io-australian-protein-design-initiative-containers-alphafold2-2.3.2-custom.img'

    // AF2-only bridge (features.pkl + provenance a3m copy). Shared ColabFold
    // a3m is already published by MMSEQS_COLABFOLDSEARCH under
    // fold/msa/mmseqs2_colabfold/.
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/af2/msas",
        mode: 'copy',
        pattern: "${meta.id}"
    )

    input:
    tuple val(meta), path(fasta), path(a3m)

    output:
    tuple val(meta), path(fasta), path("${meta.id}"), emit: msas

    script:
    """
    # The ALPHAFOLD2 predict module only reads features.pkl under
    # --use_precomputed_msas=true (see bin/fold/colabfold_a3m_to_af2_msas.py's
    # docstring for the ground-truth investigation) - the msas/ copy of the
    # raw a3m here is for provenance/debugging only, not read by AF2 itself.
    mkdir -p "${meta.id}/msas"
    cp "${a3m}" "${meta.id}/msas/colabfold.a3m"

    python ${projectDir}/bin/fold/colabfold_a3m_to_af2_msas.py \
        --fasta ${fasta} \
        --a3m ${a3m} \
        --output-dir "${meta.id}"
    """
}
