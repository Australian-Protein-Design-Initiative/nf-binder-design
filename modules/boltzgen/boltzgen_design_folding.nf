process BOLTZGEN_DESIGN_FOLDING {
    tag "batch_${start_index}"

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:71cf788'

    //publishDir path: "${params.outdir}/boltzgen/batch_${start_index}", pattern: '**', mode: 'copy'
    publishDir path: "${params.outdir}/boltzgen/batches/design_folding", pattern: '**', mode: 'copy'

    input:
    path batch_dir
    path config_yaml
    path input_files
    val design_name
    val protocol
    val start_index
    val devices
    val num_workers

    output:
    path ("batch_${start_index}"), type: 'dir', emit: batch_dir

    script:
    def config_basename = config_yaml.name
    """
    # Copy batch directory structure
    cp -r ${batch_dir}/* batch_${start_index}/ || true
    
    # Stage config.yaml (skip if already exists with same name)
    if [ ! -f ${config_basename} ]; then
        cp ${config_yaml} ${config_basename}
    fi
    
    # Stage input files at correct relative paths
    ${projectDir}/bin/boltzgen/stage_boltzgen_inputs.py ${config_basename} input_files --config-dir .
    
    # Run boltzgen design folding step
    # HF_HOME is set to /models/boltzgen in container with pre-cached weights
    boltzgen run ${config_basename} \
        --output batch_${start_index}/ \
        --protocol ${protocol} \
        --steps design_folding \
        --devices ${devices} \
        --num_workers ${num_workers} \
        --cache /models/boltzgen \
        ${task.ext.args ?: ''}
    
    # Rename files to add start_index offset
    if [ -d batch_${start_index}/intermediate_designs_inverse_folded/refold_design_cif ]; then
        ${projectDir}/bin/boltzgen/rename_boltzgen_files.py \
            batch_${start_index}/intermediate_designs_inverse_folded/refold_design_cif \
            ${design_name} \
            ${start_index}
    fi
    if [ -d batch_${start_index}/intermediate_designs_inverse_folded/fold_out_design_npz ]; then
        ${projectDir}/bin/boltzgen/rename_boltzgen_files.py \
            batch_${start_index}/intermediate_designs_inverse_folded/fold_out_design_npz \
            ${design_name} \
            ${start_index}
    fi
    """
}
