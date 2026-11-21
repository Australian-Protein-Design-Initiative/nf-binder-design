process BC2_GENERATE_CONFIG {
    tag "${design_name}"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.6'

    input:
    path input_pdb
    val design_name
    val target_chains
    val hotspot_res
    val binder_length_range

    output:
    path 'campaign.json', emit: config_json

    script:
    def chains_arg = target_chains ? "--target-chains '${target_chains}'" : ''
    def hotspot_arg = hotspot_res ? "--hotspot-res '${hotspot_res}'" : ''
    def lengths_arg = binder_length_range ? "--binder-length-range '${binder_length_range}'" : ''
    """
    ${projectDir}/bin/bindcraft2/stage_bc2_config.py generate \\
        --design-name "${design_name}" \\
        --input-pdb "${input_pdb}" \\
        ${chains_arg} \\
        ${hotspot_arg} \\
        ${lengths_arg} \\
        -o campaign.json
    """
}
