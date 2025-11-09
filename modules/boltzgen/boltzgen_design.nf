process BOLTZGEN_DESIGN {
    tag "batch_${start_index}"

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:4b1e659_filename-index-offset'

    //publishDir path: "${params.outdir}/boltzgen/batch_${start_index}", pattern: '**', mode: 'copy'
    publishDir path: "${params.outdir}/boltzgen/batches/design", pattern: '**', mode: 'copy'

    input:
    path config_yaml
    path input_files
    val design_name
    val protocol
    val n_designs
    val start_index
    val devices
    val num_workers

    output:
    path ("batch_${start_index}"), type: 'dir', emit: batch_dir

    script:
    def config_basename = config_yaml.name
    """
    # Create batch output directory
    mkdir -p batch_${start_index}
    
    # Stage config.yaml (skip if already exists with same name)
    if [ ! -f ${config_basename} ]; then
        cp ${config_yaml} ${config_basename}
    fi
    
    # Stage input files at correct relative paths
    ${projectDir}/bin/stage_boltzgen_inputs.py ${config_basename} input_files --config-dir .
    
    # Run boltzgen design step
    # HF_HOME is set to /models/boltzgen in container with pre-cached weights
    boltzgen run ${config_basename} \
        --output batch_${start_index}/ \
        --protocol ${protocol} \
        --steps design \
        --num_designs ${n_designs} \
        --devices ${devices} \
        --num_workers ${num_workers} \
        --cache /models/boltzgen \
        ${task.ext.args ?: ''}
    
    # Rename files to add start_index offset
    ${projectDir}/bin/rename_boltzgen_files.py \
        batch_${start_index}/intermediate_designs \
        ${design_name} \
        ${start_index} \
        --num-designs ${n_designs}
    
    # Ensure directory exists and is non-empty for Nextflow output detection
    touch batch_${start_index}/.nextflow_complete
    """
}
