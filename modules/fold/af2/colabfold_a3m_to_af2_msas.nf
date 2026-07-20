process COLABFOLD_A3M_TO_AF2_MSAS {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif'

    // Published alongside the native jackhmmer_af2 MSA dirs (same publishDir
    // pattern/saveAs caveat as ALPHAFOLD2_JACKHMMER_MSA: pattern operates on
    // the top-level "${meta.id}" directory item, not files nested inside it).
    publishDir(
        path: "${params.outdir}/af2/msas",
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
    # --use_precomputed_msas=true (see bin/colabfold_a3m_to_af2_msas.py's
    # docstring for the ground-truth investigation) - the msas/ copy of the
    # raw a3m here is for provenance/debugging only, not read by AF2 itself.
    mkdir -p "${meta.id}/msas"
    cp "${a3m}" "${meta.id}/msas/colabfold.a3m"

    python ${projectDir}/bin/colabfold_a3m_to_af2_msas.py \
        --fasta ${fasta} \
        --a3m ${a3m} \
        --output-dir "${meta.id}"
    """
}
