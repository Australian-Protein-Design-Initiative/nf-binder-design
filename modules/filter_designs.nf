process FILTER_DESIGNS {
    tag "Filtering ${pdb.name}"
    container 'ghcr.io/australian-protein-design-initiative/containers/mdanalysis:2.8.0'

    publishDir "${params.outdir}/filtering/backbones/scores", mode: 'copy', pattern: '*.scores.tsv'
    publishDir "${params.outdir}/filtering/backbones/accepted", mode: 'copy', pattern: 'accepted/*.pdb'
    publishDir "${params.outdir}/filtering/backbones/rejected", mode: 'copy', pattern: 'rejected/*.pdb'

    input:
    path pdb
    val filters
    val binder_chains

    output:
    path "accepted/${pdb.name}", emit: accepted, optional: true
    path "rejected/${pdb.name}", emit: rejected, optional: true
    path "${pdb.baseName}.scores.tsv", emit: scores

    script:
    def filter_args_str = filters ? filters.split(';').collect { "--filter \"${it}\"" }.join(' ') : ''
    """
    filter_designs.py \\
        ${pdb} \\
        ${filter_args_str} \\
        --binder-chains ${binder_chains} \\
        --collect-in . \\
        --output ${pdb.baseName}.scores.tsv
    """
}
