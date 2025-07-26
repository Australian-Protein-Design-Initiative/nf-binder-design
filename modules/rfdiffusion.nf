process RFDIFFUSION {
    container 'ghcr.io/australian-protein-design-initiative/containers/rfdiffusion:pytorch2407'

    publishDir path: "${params.outdir}/rfdiffusion", pattern: 'pdbs/*.pdb', mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: 'traj/*.pdb{,.gz}', mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: 'configs/*/*.yaml', mode: 'copy'
    publishDir path: "${params.outdir}/rfdiffusion", pattern: 'logs/*/*.log', mode: 'copy'

    input:
    val rfd_config_name
    path input_pdb
    val rfd_model_path
    val contigs
    val hotspot_res
    val batch_size
    val design_startnum
    val unique_id

    output:
    path 'pdbs/*.pdb', emit: pdbs
    path 'traj/*.pdb{,.gz}', emit: trajs
    path 'configs/*/*.yaml', emit: configs
    path 'logs/*/*.log', emit: logs

    script:
    def rfd_model_path_arg = rfd_model_path ? "inference.ckpt_override_path=${rfd_model_path}" : ''
    def hotspot_res_arg = hotspot_res ? "ppi.hotspot_res='${hotspot_res}'" : ''
    def config_name = rfd_config_name ? "--config-name=${rfd_config_name}" : ''
    """
    if [[ ${params.require_gpu} == "true" ]]; then
       if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
           echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi

        nvidia-smi
    fi

    # Find least-used GPU (by active processes and VRAM) and set CUDA_VISIBLE_DEVICES
    # This is a bit of a hack, but nextflow (and hyperqueue) don't have good support for
    # allocating to specific GPUs on multi-GPU nodes, so this is the solution for now.
    # There are likely still race conditions where two processes could get allocated to the
    # same GPU, but in practise it seems to generally work well enough for runs on single
    # nodes with multiple GPUs.
    # We DON'T need this for SLURM, it should correctly allocate GPU resources itself.
    if [[ -n "${params.gpu_devices}" ]]; then
        ps aux | grep "run_inference.py"
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    RUN_INF="python /app/RFdiffusion/scripts/run_inference.py"

    mkdir -p schedules

    # TODO: if rfd_config is a path to a .yml or .yaml file, use --config-path
    #       instead of --config-name

    \${RUN_INF} \
        ${config_name} \
        inference.output_prefix=outputs/${params.design_name}_${unique_id} \
        inference.input_pdb=${input_pdb} \
        contigmap.contigs='${contigs}' \
        ${hotspot_res_arg} \
        denoiser.noise_scale_ca=${params.rfd_noise_scale} \
        denoiser.noise_scale_frame=${params.rfd_noise_scale} \
        inference.num_designs=${batch_size} \
        inference.design_startnum=${design_startnum} \
        inference.schedule_directory_path=schedules \
        ${rfd_model_path_arg} \
        ${params.rfd_extra_args}

    # inference.model_directory_path=models

    mkdir -p pdbs
    mv outputs/*.pdb pdbs/
    mv outputs/traj .

    if [[ ${params.rfd_compress_trajectories} == "true" ]]; then
        gzip traj/*.pdb || true
    fi

    # Move configs and logs for nicer per-batch output
    mkdir -p configs/${design_startnum}
    mv outputs/*/*/.hydra/*.yaml configs/${design_startnum}/
    mkdir -p logs/${design_startnum}
    mv outputs/*/*/*.log logs/${design_startnum}/
    """
}
