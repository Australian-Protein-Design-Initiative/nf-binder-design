process GENERATE_RFD4_CONFIG {
    executor 'local'

    publishDir path: "${params.outdir}/rfd4_proteina", pattern: '*.json', mode: 'copy'

    input:
    val design_name
    path input_structure
    val contigs
    val hotspot_res
    val partial_trajectory_fraction

    output:
    path "${design_name}.json", emit: config_spec

    script:
    def hotspot_arg = hotspot_res ? "--hotspot-res '${hotspot_res}'" : ''
    def partial_arg = partial_trajectory_fraction ? "--partial-trajectory-fraction ${partial_trajectory_fraction}" : ''
    """
    ${projectDir}/bin/rfd4_proteina/rfd4_config.py generate \\
        --design-name "${design_name}" \\
        --input "${input_structure}" \\
        --contigs '${contigs}' \\
        ${hotspot_arg} \\
        ${partial_arg} \\
        -o ${design_name}.json
    """
}
