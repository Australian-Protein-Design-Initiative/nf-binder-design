process BINDCRAFT2 {
    tag "${batch_id}"

    container 'ghcr.io/australian-protein-design-initiative/containers/bindcraft2:v1.0.1'

    // The campaign project_folder is named after batch_id, which is unique
    // across inputs and batch indices, so batches never collide when the
    // gather step stages them all into one parent.
    publishDir path: "${params.outdir}/bc2/batches", pattern: { "${batch_id}" }, mode: 'copy'
    publishDir path: { "${params.outdir}/bc2/batches/${batch_id}" }, pattern: 'bc2.log', mode: 'copy'
    publishDir path: { "${params.outdir}/bc2/batches/${batch_id}" }, pattern: 'bc2_campaign.staged.json', mode: 'copy'
    //
    // Do not add a publishDir here to flatten 3_Ranked/*.cif. Declaring the
    // batch directory as an output shadows the files nested inside it for
    // publishing, so such a pattern matches nothing and fails silently -- the
    // nested output channels below still work, only publishing does not.
    // BC2_CAMPAIGN_OUTPUT flattens them instead, which is also the correct
    // point to do it, since campaign_output re-ranks each campaign.

    input:
    tuple val(batch_id), val(bc2_args), val(stage_args), path(config_json), path(input_files)
    path metadata_json

    output:
    tuple val(batch_id), path(batch_id, type: 'dir'),             emit: batch_dir
    path "${batch_id}/3_Ranked/!_Ranked.csv",    optional: true,  emit: ranked_csv
    path "${batch_id}/2_Refolded/!_Refolded.csv", optional: true, emit: refolded_csv
    path "${batch_id}/1_Trajectories/!_Trajectories.csv", optional: true, emit: trajectory_csv
    path "${batch_id}/campaign_metadata.json",   optional: true,  emit: campaign_metadata
    path 'bc2_campaign.staged.json',                              emit: staged_config
    path 'bc2.log',                                               emit: log

    script:
    // metadata_json is optional; the workflow passes a placeholder file when unset.
    def metadata_arg = (metadata_json && metadata_json.name != 'NO_BC2_METADATA') ? "--metadata ${metadata_json}" : ''
    def af2_params_export = params.bc2_af2_params
        ? "export BINDCRAFT_AF2_PARAMS=\"${params.bc2_af2_params}\""
        : '# AlphaFold parameters are baked into the image at /opt/bindcraft/bindcraft/weights/alphafold'
    """
    set -euo pipefail

    if [[ ${params.require_gpu} == "true" ]]; then
       if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
           echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi

        nvidia-smi
    fi

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

    # BC2 writes caches under \$HOME by default, which is read-only or shared in
    # a container. Keep them in the task directory. This does not affect weight
    # lookup: the baked parameters are found inside the package, an earlier tier
    # than the XDG cache (bindcraft/model_weights.py).
    export XDG_CACHE_HOME="\$(pwd)/.cache"
    export JAX_COMPILATION_CACHE_DIR="\$(pwd)/.jax_cache"
    export MPLCONFIGDIR="\$(pwd)/.matplotlib"
    mkdir -p "\$XDG_CACHE_HOME" "\$JAX_COMPILATION_CACHE_DIR" "\$MPLCONFIGDIR"
    ${af2_params_export}

    # Rewrite targets[].target_path and binder_scaffold to basenames so they
    # resolve against the campaign JSON's own directory, which is this work
    # directory where Nextflow staged them.
    #
    # The output name must not collide with the staged input's. Nextflow stages
    # inputs as symlinks, so writing to the input's own name would follow the
    # link and overwrite the user's source config (or, in params mode, the
    # upstream task's cached output). BC2 then resolves target paths with
    # Path(...).resolve(), which follows symlinks, so a symlinked campaign file
    # would also resolve relative to the original's directory rather than this
    # one, and the staged inputs would not be found.
    #
    # Called through python3 rather than the shebang: this container ships only
    # /usr/bin/python3, with no `python` on PATH for `/usr/bin/env python` to find.
    # stage_args carries --hotspot-subsample when it is in play, and is empty
    # otherwise. Subsampling is done here rather than at config generation time
    # so it reaches config mode too, where the pipeline never writes the JSON.
    python3 ${projectDir}/bin/bindcraft2/stage_bc2_config.py stage ${config_json} \\
        ${stage_args} -o bc2_campaign.staged.json

    # pipefail is set above so a bindcraft crash still fails the task despite
    # `| tee` (the v1 module had exactly this bug, fixed in 0.3.1).
    #
    # BindCraft2 then intermittently dies during C++/JAX/CUDA interpreter
    # teardown, after the campaign has already completed and flushed its
    # outputs: exit 139 (SIGSEGV) or 134 (SIGABRT, "pure virtual method
    # called" / "corrupted double-linked list"). `set -e` would discard a
    # finished batch. Tolerate only that case. Any other status, or the same
    # status before `campaign done:` or without the summary files, still fails.
    rc=0
    bindcraft design bc2_campaign.staged.json \\
        ${bc2_args} \\
        ${metadata_arg} \\
        ${params.bc2_extra_args ?: ''} \\
        ${task.ext.args ?: ''} 2>&1 | tee bc2.log || rc=\$?

    if [[ "\$rc" -ne 0 ]]; then
        if [[ "\$rc" -eq 134 || "\$rc" -eq 139 ]] &&
            grep -qF 'campaign done:' bc2.log &&
            [[ -s "${batch_id}/summary.csv" && -s "${batch_id}/campaign_metadata.json" ]]; then
            echo "WARNING: tolerated a post-completion JAX/CUDA teardown crash (bindcraft design exited \$rc after printing 'campaign done:'; ${batch_id}/summary.csv and campaign_metadata.json are present and non-empty). Continuing." >&2
        else
            exit "\$rc"
        fi
    fi

    # A campaign that exits 0 having written nothing would otherwise surface
    # much later as an empty gather.
    if [[ ! -d "${batch_id}" ]]; then
        echo "ERROR: bindcraft produced no ${batch_id} project folder" >&2
        exit 1
    fi
    """
}
