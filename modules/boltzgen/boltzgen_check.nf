process BOLTZGEN_CHECK {

    container '/home/perry/projects/nf-binder-design/repos/boltzgen/boltzgen-file-index-offset.sif'

    input:
    path 'input/*'

    output:
    path "${params.design_name}.cif"

    script:
    """
    boltzgen check ${params.design_name}.yaml
    """
}
