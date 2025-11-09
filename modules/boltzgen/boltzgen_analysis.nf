process BOLTZGEN_ANALYSIS {
    tag "analysis"

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:71cf788'

    publishDir path: "${params.outdir}/boltzgen", pattern: '**', mode: 'copy'

    input:
    path merged_dir
    path config_yaml
    path input_files
    val design_name
    val protocol

    output:
    path 'merged', type: 'dir', emit: merged_dir

    script:
    def config_basename = config_yaml.name
    """
    # Copy merged directory structure
    cp -r ${merged_dir}/* merged/ || true
    
    # Stage config.yaml (skip if already exists with same name)
    if [ ! -f ${config_basename} ]; then
        cp ${config_yaml} ${config_basename}
    fi
    
    # Stage input files at correct relative paths
    ${projectDir}/bin/stage_boltzgen_inputs.py ${config_basename} input_files --config-dir .
    
    # Run boltzgen analysis step
    # HF_HOME is set to /models/boltzgen in container with pre-cached weights
    boltzgen run ${config_basename} \
        --output merged/ \
        --protocol ${protocol} \
        --steps analysis \
        --cache /models/boltzgen \
        ${task.ext.args ?: ''}
    """
}
