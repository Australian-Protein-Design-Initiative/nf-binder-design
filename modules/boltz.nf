process BOLTZ {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1'
    publishDir "${params.outdir}/${protocol}/${type}", mode: 'copy'

    input:
    tuple val(meta), path(yaml_file), path(target_msa), path(binder_msa)
    path templates
    tuple val(protocol), val(type)

    output:
    path ("boltz_results_${meta.id}"), emit: results
    tuple val(meta), path("boltz_results_${meta.id}/predictions/${meta.id}/*.cif"), emit: predicted_structure
    tuple val(meta), path("boltz_results_${meta.id}/predictions/${meta.id}/confidence_${meta.id}_model_0.json"), emit: confidence_json

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
    export BOLTZ_CACHE=\${BOLTZ_CACHE:-/app/boltz/cache}

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
        --cache \$BOLTZ_CACHE \
        ${yaml_file}
        #--out_dir "results"
        #--output_format pdb
    """
}
