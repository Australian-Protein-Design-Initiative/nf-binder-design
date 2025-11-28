process BOLTZ {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1'
    publishDir "${params.outdir}/${step_name}", mode: 'copy'

    input:
    tuple val(meta), path(yaml_file), path(target_msa), path(binder_msa)
    path templates
    val step_name

    output:
    path ("boltz_results_${yaml_file.baseName}"), emit: results
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/${yaml_file.baseName}_model_0.pdb"), emit: pdb
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/confidence_${yaml_file.baseName}_model_0.json"), emit: confidence_json
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/*_ipsae.tsv"), emit: ipsae_tsv
    tuple val(meta), path("boltz_results_${yaml_file.baseName}/predictions/${yaml_file.baseName}/*_ipsae_byres.tsv"), emit: ipsae_byres_tsv

    script:
    def use_msa_server_flag = params.use_msa_server ? '--use_msa_server' : ''
    def args = task.ext.args ?: ''
    """
    # Find least-used GPU (by active processes and VRAM) and set CUDA_VISIBLE_DEVICES
    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    # Boltz model weights are stored in our container
    export BOLTZ_CACHE=/app/boltz/cache

    # Create various tmp/cache directories that are expected to be in \$HOME by default
    export NUMBA_CACHE_DIR="\$(pwd)/.numba_cache"
    mkdir -p \$NUMBA_CACHE_DIR
    export XDG_CONFIG_HOME="\$(pwd)/.config"
    mkdir -p \$XDG_CONFIG_HOME
    export TRITON_CACHE_DIR="\$(pwd)/.triton_cache"
    mkdir -p \$TRITON_CACHE_DIR

    # Prevent Python from using ~/.local/lib/ packages mounted inside the container
    export PYTHONNOUSERSITE=1

    # We could autodetect if we have a GPU, but lets leave this up to task.ext.args
    # instead of using --accelerator \${ACCELERATOR} \
    # if nvidia-smi >/dev/null 2>&1; then
    #    ACCELERATOR=gpu
    # else
    #     ACCELERATOR=cpu
    # fi

    boltz predict \
        ${args} \
        ${use_msa_server_flag} \
        --preprocessing-threads ${task.cpus} \
        --num_workers ${task.cpus} \
        --output_format pdb \
        ${yaml_file}

    ${projectDir}/bin/ipsae.py \
        boltz_results_*/predictions/*/pae_*.npz \
        boltz_results_*/predictions/*/*.pdb \
        10 10
    """
}
