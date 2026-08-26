// Assemble AF2 multimer precomputed-MSA tree for one target+binder pair.
// Target MSA: jackhmmer dir and/or ColabFold a3m for chain A; binder B is
// always query-only. Pass assets/dummy_files/empty for unused slots.
process FOLD_ASSEMBLE_AF2_MULTIMER_MSAS {
    tag "${meta.id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

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
    python3 ${projectDir}/bin/fold/assemble_af2_multimer_msas.py \\
        --pair-fasta ${pair_fasta} \\
        --pair-id '${meta.id}' \\
        --target-msa-dir ${target_msa_dir} \\
        --target-a3m ${target_a3m} \\
        -o .
    """
}
