process BOLTZGEN_INVERSE_FOLDING {
    tag "batch_${start_index}"

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:4b1e659_filename-index-offset'

    //publishDir path: "${params.outdir}/boltzgen/batch_${start_index}", pattern: '**', mode: 'copy'
    publishDir path: "${params.outdir}/boltzgen/batches/inverse_folding", pattern: '**', mode: 'copy'

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
    ${projectDir}/bin/stage_boltzgen_inputs.py ${config_basename} input_files --config-dir .
    
    # Run boltzgen inverse folding step
    # HF_HOME is set to /models/boltzgen in container with pre-cached weights
    boltzgen run ${config_basename} \
        --output batch_${start_index}/ \
        --protocol ${protocol} \
        --steps inverse_folding \
        --start_index ${start_index} \
        --devices ${devices} \
        --num_workers ${num_workers} \
        --cache /models/boltzgen \
        ${task.ext.args ?: ''}
    """
}
