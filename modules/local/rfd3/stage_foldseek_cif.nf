process RFD3_STAGE_FOLDSEEK_CIF {
    executor 'local'

    input:
    tuple val(meta), path(cif)

    output:
    path "${meta.id}_design.cif"

    script:
    """
    cp ${cif} ${meta.id}_design.cif
    """
}
