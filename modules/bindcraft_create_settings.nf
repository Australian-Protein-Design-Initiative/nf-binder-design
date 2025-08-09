process BINDCRAFT_CREATE_SETTINGS {
    container 'ghcr.io/australian-protein-design-initiative/containers/mdanalysis:2.8.0'

    input:
    tuple val(batch_id), val(n_designs)
    path input_pdb
    val hotspot_res
    val target_chains
    val binder_length_range
    val design_name
    val hotspot_subsample

    output:
    path 'settings.json', emit: settings_json
    val batch_id,         emit: batch_id

    script:
    """
    ${baseDir}/bin/create_bindcraft_settings.py \\
        --input_pdb ./${input_pdb} \\
        --hotspot_res "${hotspot_res}" \\
        --target_chains "${target_chains}" \\
        --binder_length_range "${binder_length_range}" \\
        --design_name "${design_name}_${batch_id}" \\
        --n_designs ${n_designs} \\
        --output_dir "results" \\
        --output_json "settings.json" \\
        --hotspot_subsample ${hotspot_subsample}
    """
}
