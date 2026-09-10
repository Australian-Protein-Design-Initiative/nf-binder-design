process BOLTZ {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1-2'
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
    # Claim a GPU for this task's lifetime, then record which card we got
    # (bin/gpu_lock.sh). The claim is required, and fails the task if it cannot
    # be made. The recording is diagnostic, and must never fail the task -- the
    # `|| true` also suspends `set -e` for the whole function body, so nothing
    # inside it can abort the script either.
    if [[ -n "${params.gpu_devices}" ]]; then
        source ${projectDir}/bin/gpu_lock.sh
        nfbd_acquire_gpu "${params.gpu_devices}" "${params.gpu_lock_dir ?: workDir.toString() + '/.gpu_locks'}" ${task.ext.gpu_slots ?: params.gpu_slots_per_device} ${params.gpu_lock_timeout} || exit 1
    else
        source ${projectDir}/bin/gpu_lock.sh || true
    fi
    nfbd_record_gpu_trace "${params.gpu_trace_dir ?: workDir.toString() + '/.gpu_trace'}" "${task.process}" || true

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

    BOLTZ_PREDICT_LOG=.boltz_predict_console.log
    rm -f "\$BOLTZ_PREDICT_LOG"
    set +e
    boltz predict \
        ${args} \
        ${use_msa_server_flag} \
        --preprocessing-threads ${task.cpus} \
        --num_workers ${task.cpus} \
        --output_format pdb \
        ${yaml_file} 2>&1 | tee "\$BOLTZ_PREDICT_LOG"
    boltz_rc=\${PIPESTATUS[0]}
    set -e
    if grep -qF 'ran out of memory, skipping batch' "\$BOLTZ_PREDICT_LOG"; then
        echo 'BOLTZ: Boltz logged GPU OOM (batch skipped); failing.' >&2
        exit 1
    fi
    if [[ "\$boltz_rc" -ne 0 ]]; then
        exit "\$boltz_rc"
    fi

    ${projectDir}/bin/ipsae.py \\
        --update-summary boltz_results_*/predictions/*/confidence_*_model_0.json \\
        --format boltz \\
        boltz_results_*/predictions/*/pae_*.npz \\
        boltz_results_*/predictions/*/*.pdb \\
        10 10
    """
}
