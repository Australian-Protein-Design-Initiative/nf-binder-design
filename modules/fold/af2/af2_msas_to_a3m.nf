process AF2_MSAS_TO_A3M {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif'

    // Shared jackhmmer-derived a3m for Boltz/RF3/Protenix under the same tree as
    // the native AF2 jackhmmer MSA dirs.
    publishDir(
        path: "${params.outdir}/fold/msa/jackhmmer_af2",
        mode: 'copy',
        pattern: '*.a3m'
    )

    input:
    tuple val(meta), path(fasta), path(af2_msa_dir)

    output:
    tuple val(meta), path(fasta), path("${meta.id}.a3m"), emit: a3m

    script:
    """
    python ${projectDir}/bin/af2_msas_to_a3m.py \
        --msas-dir "${af2_msa_dir}/msas" \
        --output "${meta.id}.a3m"
    """
}
