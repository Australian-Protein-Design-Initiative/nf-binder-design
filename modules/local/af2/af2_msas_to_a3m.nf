process AF2_MSAS_TO_A3M {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif'

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
