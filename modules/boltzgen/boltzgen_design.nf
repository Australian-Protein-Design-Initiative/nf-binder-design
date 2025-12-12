process BOLTZGEN_DESIGN {
    tag "batch_${start_index}"

    container 'ghcr.io/australian-protein-design-initiative/containers/boltzgen:0.2.0'

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
    def devices_arg = devices ? "--devices ${devices}" : ""
    def num_workers_arg = num_workers ? "--num_workers ${num_workers}" : ""
    """
    set -euo pipefail
    
    nvidia-smi

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

    # Create batch output directory
    mkdir -p batch_${start_index}
    
    # Stage config.yaml (skip if already exists with same name)
    if [ ! -f ${config_basename} ]; then
        cp ${config_yaml} ${config_basename}
    fi
    
    # Stage input files at correct relative paths
    ${projectDir}/bin/boltzgen/stage_boltzgen_inputs.py ${config_basename} input_files --config-dir .

    # Run boltzgen design step
    # HF_HOME is set to /models/boltzgen in container with pre-cached weights
    boltzgen run ${config_basename} \
        --output batch_${start_index}/ \
        --protocol ${protocol} \
        --steps design \
        --num_designs ${n_designs} \
        ${devices_arg} \
        ${num_workers_arg} \
        --cache /models/boltzgen \
        \${BOLTZGEN_USE_KERNELS_FLAG} \
        ${task.ext.args ?: ''}

    if grep -q "WARNING: ran out of memory, skipping batch" .command.log; then
        echo "ERROR: Task failed due to out of memory condition" >&2
        exit 1
    fi
    
    # Rename files to add start_index offset
    ${projectDir}/bin/boltzgen/rename_boltzgen_files.py \
        batch_${start_index}/intermediate_designs \
        ${design_name} \
        ${start_index} \
        --num_designs ${n_designs}
    
    # Ensure directory exists and is non-empty for Nextflow output detection
    touch batch_${start_index}/.nextflow_complete
    """
}