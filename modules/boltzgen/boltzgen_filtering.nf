process BOLTZGEN_FILTERING {
    tag "filtering"

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:0.2.0'

    // Only publish the final filtering outputs to avoid collisions when re-filtering an existing merged directory
    publishDir path: "${params.outdir}/boltzgen", pattern: 'filtered/final_ranked_designs', mode: 'copy', overwrite: true

    input:
    path merged_dir
    path config_yaml
    path input_files
    val design_name
    val protocol
    val budget
    val filtering_arg_list

    output:
    path 'filtered/final_ranked_designs', type: 'dir', emit: final_ranked_designs_dir
    path 'filtered/final_ranked_designs/config/filtering.yaml', type: 'file', emit: filtering_config_yaml
    path 'gpu_stats.csv', optional: true, emit: gpu_stats

    script:
    def config_basename = config_yaml.name
    """
    set -euo pipefail

    # Create various tmp/cache directories that are expected to be in \$HOME by default
    export NUMBA_CACHE_DIR="\$(pwd)/.numba_cache"
    mkdir -p \$NUMBA_CACHE_DIR
    export XDG_CONFIG_HOME="\$(pwd)/.config"
    mkdir -p \$XDG_CONFIG_HOME
    export TRITON_CACHE_DIR="\$(pwd)/.triton_cache"
    mkdir -p \$TRITON_CACHE_DIR
    export CUEQ_TRITON_CACHE_DIR="\$(pwd)/.cuequivariance_triton_cache"
    mkdir -p \$CUEQ_TRITON_CACHE_DIR

    nvidia-smi

    # Start GPU monitoring in background if enabled
    if [[ "${params.enable_gpu_stats}" == "true" ]]; then
        PARENT_DIR=\$(basename \$(dirname \$(pwd)))
        CURRENT_DIR=\$(basename \$(pwd))
        TASK_HASH="\${PARENT_DIR}/\${CURRENT_DIR}"
        TASK_HASH="\${TASK_HASH:0:9}"
        ${baseDir}/bin/monitor_gpu.py \
            --process-name "BOLTZGEN_FILTERING" \
            --task-hash "\${TASK_HASH}" \
            --task-index "${task.index}" \
            --interval ${params.gpu_stats_interval} \
            --output gpu_stats.csv &
        GPU_MONITOR_PID=\${!:-}
        if [[ -n "\${GPU_MONITOR_PID}" ]]; then
            trap "kill \${GPU_MONITOR_PID} 2>/dev/null || true" EXIT
        fi
    fi

    # Auto-detect MIG GPU
    BOLTZGEN_USE_KERNELS_FLAG=""
    MIG_UUID=\$(nvidia-smi -L 2>/dev/null | sed -n 's/.*(UUID: \\(MIG-[^)]*\\)).*/\\1/p' | head -n 1 || true)
    if [ -n "\$MIG_UUID" ]; then
    cat << 'EOF' > /tmp/mig_patch.py
import pynvml
import sys

try:
    # BoltzGen v0.1.4 has a bug that causes it to crash when MIG is used.
    # We override the function that crashes on MIG slices to return this constant.
    # (A100 has 6912 CUDA cores - in theory this number should work okay for other 
    #  GPUs irrespective of core count)
    pynvml.nvmlDeviceGetNumGpuCores = lambda handle: 6912
    print(">>> MIG PATCH: Successfully mocked nvmlDeviceGetNumGpuCores", file=sys.stderr)
except Exception as e:
    print(f">>> MIG PATCH: Failed to mock pynvml: {e}", file=sys.stderr)
EOF

        export PYTHONSTARTUP=/tmp/mig_patch.py
        export BOLTZGEN_USE_KERNELS_FLAG="--use_kernels false"
    fi

    # Copy merged directory structure
    mkdir -p filtered
    cp -r ${merged_dir}/* filtered/ || true
    
    # Stage config.yaml (skip if already exists with same name)
    if [ ! -f ${config_basename} ]; then
        cp ${config_yaml} ${config_basename}
    fi
    
    # Stage input files at correct relative paths
    ${projectDir}/bin/boltzgen/stage_boltzgen_inputs.py ${config_basename} input_files --config-dir .
    
    # Run boltzgen filtering step
    # HF_HOME is set to /models/boltzgen in container with pre-cached weights
    boltzgen run ${config_basename} \
        --output filtered/ \
        --protocol ${protocol} \
        --steps filtering \
        --budget ${budget} \
        --cache /models/boltzgen \
        \${BOLTZGEN_USE_KERNELS_FLAG} \
        ${filtering_arg_list.join(' ')} \
        ${task.ext.args ?: ''}
    
    # Move filtering.yaml to final_ranked_designs/config/ so when running filtering 
    # multiple times we can keep the config with that output
    mkdir -p filtered/final_ranked_designs/config
    if [ -f filtered/config/filtering.yaml ]; then
        cp filtered/config/filtering.yaml filtered/final_ranked_designs/config/filtering.yaml
    fi
    """
}
