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
    nvidia-smi

    # Auto-detect MIG GPU
    MIG_UUID=\$(nvidia-smi -L 2>/dev/null | sed -n 's/.*(UUID: \\(MIG-[^)]*\\)).*/\\1/p' | head -n 1)
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
    cp -r ${merged_dir}/* merged/ || true
    
    # Stage config.yaml (skip if already exists with same name)
    if [ ! -f ${config_basename} ]; then
        cp ${config_yaml} ${config_basename}
    fi
    
    # Stage input files at correct relative paths
    ${projectDir}/bin/boltzgen/stage_boltzgen_inputs.py ${config_basename} input_files --config-dir .
    
    # Run boltzgen analysis step
    # HF_HOME is set to /models/boltzgen in container with pre-cached weights
    boltzgen run ${config_basename} \
        --output merged/ \
        --protocol ${protocol} \
        --steps analysis \
        --cache /models/boltzgen \
        \${BOLTZGEN_USE_KERNELS_FLAG} \
        ${task.ext.args ?: ''}
    """
}
