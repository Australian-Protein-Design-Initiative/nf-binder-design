process TRIM_TO_CONTIGS {

    container 'ghcr.io/australian-protein-design-initiative/containers/mdanalysis:2.8.0'

    input:
    path input_pdb
    val contigs

    output:
    path "*.pdb", emit: pdb

    script:
    """
    ${baseDir}/bin/trim_to_contigs.py \\
        ${input_pdb} \\
        "${contigs}" \\
        --output "${input_pdb.simpleName}_contigs.pdb"
    """
}