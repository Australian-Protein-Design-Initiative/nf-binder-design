process RFDIFFUSION_PARTIAL {
    container 'ghcr.io/australian-protein-design-initiative/containers/rfdiffusion:pytorch2407'

    publishDir path: "${params.outdir}/rfdiffusion", pattern: 'pdbs/partial/*/*.pdb', mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: 'traj/partial/*/*.pdb{,.gz}', mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: 'configs/partial/*/*/*.yaml', mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: 'logs/partial/*/*/*.log', mode: 'copy'

    input:
    path rfd_config_name
    path input_pdb
    val rfd_model_path
    val contigs
    val hotspot_res
    val batch_size
    val design_startnum
    val unique_id
    val partial_T

    output:
    path 'pdbs/partial/*/*.pdb', emit: pdbs
    path 'traj/partial/*/*.pdb{,.gz}', emit: trajs
    path 'configs/partial/*/*/*.yaml', emit: configs
    path 'logs/partial/*/*/*.log', emit: logs
    path 'gpu_stats.csv', optional: true, emit: gpu_stats

    script:
    // we replace all '.' characters with '_' in the output name to avoid filename collisions
    // downstream, since af2_initial_guess splits on the first '.'
    def outname = "${input_pdb.baseName}_partial_T${partial_T}".replaceAll('\\.', '_')

    def rfd_model_path_arg = rfd_model_path ? "inference.ckpt_override_path=${rfd_model_path}" : ''
    def hotspot_res_arg = hotspot_res ? "ppi.hotspot_res='${hotspot_res}'" : ''

    """
    if [[ ${params.require_gpu} == "true" ]]; then
       if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
           echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi
        nvidia-smi
    fi

    # Find least-used GPU (by active processes and VRAM) and set CUDA_VISIBLE_DEVICES
    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    # Start GPU monitoring in background if enabled
    if [[ "${params.enable_gpu_stats}" == "true" ]]; then
        PARENT_DIR=\$(basename \$(dirname \$(pwd)))
        CURRENT_DIR=\$(basename \$(pwd))
        TASK_HASH="\${PARENT_DIR}/\${CURRENT_DIR}"
        TASK_HASH="\${TASK_HASH:0:9}"
        ${baseDir}/bin/monitor-gpu.py \
            --process-name "RFDIFFUSION_PARTIAL" \
            --task-hash "\${TASK_HASH}" \
            --task-index "${task.index}" \
            --interval ${params.gpu_stats_interval} \
            --output gpu_stats.csv &
        GPU_MONITOR_PID=\${!:-}
        if [[ -n "\${GPU_MONITOR_PID}" ]]; then
            trap "kill \${GPU_MONITOR_PID} 2>/dev/null || true" EXIT
        fi
    fi

    RUN_INF="python /app/RFdiffusion/scripts/run_inference.py"

    mkdir -p schedules

    # TODO: if rfd_config is a path to a .yml or .yaml file, use --config-path
    #       instead of --config-name

    # --config-name=${rfd_config_name}

    \${RUN_INF} \
        diffuser.partial_T=${partial_T} \
        inference.output_prefix=outputs/${outname}_${unique_id} \
        inference.input_pdb=${input_pdb} \
        contigmap.contigs='${contigs}' \
        ${hotspot_res_arg} \
        denoiser.noise_scale_ca=${params.rfd_noise_scale} \
        denoiser.noise_scale_frame=${params.rfd_noise_scale} \
        inference.num_designs=${batch_size} \
        inference.design_startnum=${design_startnum} \
        inference.schedule_directory_path=\$(pwd)/schedules \
        ${rfd_model_path_arg} \
        ${params.rfd_extra_args}

    # inference.model_directory_path=models

    mkdir -p pdbs/partial/${input_pdb.baseName} traj/partial/${input_pdb.baseName}
    mv outputs/*.pdb pdbs/partial/${input_pdb.baseName}/
    mv outputs/traj/*.pdb traj/partial/${input_pdb.baseName}/

    if [[ ${params.rfd_compress_trajectories} == "true" ]]; then
        gzip traj/partial/${input_pdb.baseName}/*.pdb || true
    fi

    # Move configs and logs for nicer per-batch output
    mkdir -p configs/partial/${input_pdb.baseName}/${design_startnum}
    mv outputs/*/*/.hydra/*.yaml configs/partial/${input_pdb.baseName}/${design_startnum}/
    mkdir -p logs/partial/${input_pdb.baseName}/${design_startnum}/
    mv outputs/*/*/*.log logs/partial/${input_pdb.baseName}/${design_startnum}/
    """
}
