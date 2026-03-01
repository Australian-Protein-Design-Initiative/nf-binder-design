process GENERATE_RFD3_CONFIG {
    executor 'local'

    input:
    val design_name
    path input_pdb
    val contigs
    val hotspot_res
    val partial_t
    val is_non_loopy  // null = omit from config, true = add true, false = add false

    output:
    path "${design_name}.json", emit: config_json

    script:
    def hotspot_arg = hotspot_res ? "--hotspot-res '${hotspot_res}'" : ''
    def partial_t_arg = partial_t ? "--partial-t ${partial_t}" : ''
    def is_non_loopy_arg = (is_non_loopy == null) ? '' : (is_non_loopy ? '--is-non-loopy' : '--allow-loopy')
    """
    ${projectDir}/bin/rfd3/stage_rfd3_config.py generate \
        --design-name "${design_name}" \
        --input-pdb "${input_pdb}" \
        --contigs '${contigs}' \
        ${hotspot_arg} \
        ${partial_t_arg} \
        ${is_non_loopy_arg} \
        -o ${design_name}.json
    """
}
