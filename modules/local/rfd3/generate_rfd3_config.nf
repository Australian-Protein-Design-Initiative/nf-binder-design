process GENERATE_RFD3_CONFIG {
    executor 'local'

    input:
    val design_name
    path input_pdb
    val contigs
    val hotspot_res
    val partial_t
    val allow_loopy

    output:
    path "${design_name}.json", emit: config_json

    script:
    def hotspot_arg = hotspot_res ? "--hotspot-res '${hotspot_res}'" : ''
    def partial_t_arg = partial_t ? "--partial-t ${partial_t}" : ''
    def allow_loopy_arg = allow_loopy ? '--allow-loopy' : ''
    """
    ${projectDir}/bin/rfd3/stage_rfd3_config.py generate \
        --design-name "${design_name}" \
        --input-pdb "${input_pdb}" \
        --contigs '${contigs}' \
        ${hotspot_arg} \
        ${partial_t_arg} \
        ${allow_loopy_arg} \
        -o ${design_name}.json
    """
}
