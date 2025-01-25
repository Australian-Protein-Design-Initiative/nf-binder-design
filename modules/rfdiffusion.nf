process RFDIFFUSION {
    publishDir path: "${params.outdir}/rfdiffusion", pattern: "pdbs/*.pdb", mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: "traj/*.pdb", mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: "configs/*/*.yaml", mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: "logs/*/*.log", mode: 'copy'

    input:
    path rfd_config_name
    path input_pdb
    path rfd_model_path
    val contigs
    val hotspot_res
    val batch_size
    val design_startnum

    output:
    path "pdbs/*.pdb", emit: pdbs
    path "traj/*.pdb", emit: trajs
    path "configs/*/*.yaml", emit: configs
    path "logs/*/*.log", emit: logs
    
    script:
    """
    if [[ ${params.require_gpu} == "true" ]]; then
       if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
           echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi
    fi

    nvidia-smi

    RUN_INF="python3.9 /app/RFdiffusion/scripts/run_inference.py"
    
    mkdir -p schedules

    # TODO: if rfd_config is a path to a .yml or .yaml file, use --config-path
    #       instead of --config-name

    \${RUN_INF} \
        --config-name=${rfd_config_name} \
        inference.output_prefix=outputs/${params.design_name} \
        inference.input_pdb=${input_pdb} \
        contigmap.contigs='${contigs}' \
        ppi.hotspot_res='${hotspot_res}' \
        inference.num_designs=${batch_size} \
        inference.design_startnum=${design_startnum} \
        denoiser.noise_scale_ca=0 \
        denoiser.noise_scale_frame=0 \
        inference.ckpt_override_path=${rfd_model_path} \
        inference.schedule_directory_path=schedules \
        ${params.rfd_extra_args}
    
    # inference.design_startnum=0
    # inference.model_directory_path=models 
    # inference.schedule_directory_path=schedules

    mkdir -p pdbs
    mv outputs/*.pdb pdbs/
    mv outputs/traj .
    
    # Move configs and logs for nicer per-batch output
    mkdir -p configs/${design_startnum}
    mv outputs/*/*/.hydra/*.yaml configs/${design_startnum}/    
    mkdir -p logs/${design_startnum}
    mv outputs/*/*/*.log logs/${design_startnum}/
    """
} 