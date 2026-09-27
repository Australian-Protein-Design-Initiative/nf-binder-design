process AF2_MSAS_TO_A3M {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/ghcr.io-australian-protein-design-initiative-containers-alphafold2-2.3.2-custom.img'

    // Shared jackhmmer-derived a3m for Boltz/RF3/Protenix under the same tree as
    // the native AF2 jackhmmer MSA dirs.
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/msa/jackhmmer_af2",
        mode: 'copy',
        pattern: '*.a3m'
    )

    input:
    tuple val(meta), path(fasta), path(af2_msa_dir)

    output:
    tuple val(meta), path(fasta), path("${meta.id}.a3m"), emit: a3m

    script:
    """
    set -euo pipefail
    MSAS="${af2_msa_dir}/msas"
    if [[ ! -f "\${MSAS}/bfd_uniref_hits.a3m" && ! -f "\${MSAS}/uniref90_hits.sto" ]]; then
        if [[ -d "\${MSAS}/A" ]]; then
            MSAS="\${MSAS}/A"
        fi
    fi
    python ${projectDir}/bin/fold/af2_msas_to_a3m.py \
        --msas-dir "\${MSAS}" \
        --output "${meta.id}.a3m"
    """
}
