process BINDCRAFT {
    container 'ghcr.io/australian-protein-design-initiative/containers/bindcraft:0366085-nv-cuda120'

    publishDir path: "${params.outdir}/bindcraft/batches/${batch_id}", pattern: 'results/**', mode: 'copy'
    publishDir path: "${params.outdir}/bindcraft/batches/${batch_id}", pattern: '*.{pdb,pdb.gz,json}', mode: 'copy'
    publishDir path: "${params.outdir}/bindcraft/batches/${batch_id}", pattern: 'bindcraft.log', mode: 'copy'
    publishDir(
        path: "${params.outdir}/bindcraft/accepted",
        pattern: 'results/Accepted/*.{pdb,pdb.gz}',
        mode: 'copy',
        saveAs: { filename -> java.nio.file.Paths.get(filename).getFileName().toString() }
    )

    input:
    tuple val(batch_id), path(input_pdb), path(settings_json)
    val advanced_settings_preset
    val filters_preset
    val compress_html
    val compress_pdb

    output:
    path 'batches/*', type: 'dir', followLinks: true, optional: true,       emit: batch_dir
    path 'results/Accepted/*.{pdb,pdb.gz}', optional: true,                 emit: accepted_pdbs
    path 'results/Rejected/*.{pdb,pdb.gz}', optional: true,                 emit: rejected_pdbs
    path 'results/Trajectory/Relaxed/*.{pdb,pdb.gz}', optional: true,       emit: relaxed_pdbs
    path 'results/Trajectory/LowConfidence/*.{pdb,pdb.gz}', optional: true, emit: low_confidence_pdbs
    path 'results/Trajectory/Clashing/*.{pdb,pdb.gz}', optional: true,      emit: clashing_pdbs
    path 'results/final_design_stats.csv', optional: true,         emit: final_stats_csv
    path 'results/trajectory_stats.csv', optional: true,           emit: trajectory_stats_csv
    path 'results/mpnn_design_stats.csv', optional: true,          emit: mpnn_design_stats_csv
    path 'results/failure_csv.csv', optional: true,                emit: failure_csv
    path 'bindcraft.log', optional: true,                          emit: bindcraft_log
    path '*.json', followLinks: true, includeInputs: true, optional: true, emit: settings_files
    path '*.pdb', followLinks: true, includeInputs: true, optional: true, emit: input_pdb
    path 'results/**', emit: all_results

    script:
    def advanced_settings_filename = advanced_settings_preset ? "/app/BindCraft/settings_advanced/${advanced_settings_preset}.json" : '/app/BindCraft/settings_advanced/default_4stage_multimer.json'
    def modified_advanced_settings_filename = "./${file(advanced_settings_filename).getName()}"
    def filters_filename = filters_preset ? "/app/BindCraft/settings_filters/${filters_preset}.json" : '/app/BindCraft/settings_filters/default_filters.json'
    def modified_filters_filename = "./${file(filters_filename).getName()}"
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

    ##
    # We modify the advanced settings to set the `max_trajectories` to the batch size
    # This way Nextflow can run a defined number of trajectories per task
    ##
    #######################################################################
    /opt/conda/envs/BindCraft/bin/python -c '
import json
import os

with open("${advanced_settings_filename}", "r") as f:
    settings = json.load(f)

settings["max_trajectories"] = ${params.bindcraft_batch_size}

with open("${modified_advanced_settings_filename}", "w") as f:
    json.dump(settings, f, indent=4)
'
    #######################################################################

    # Keep a copy of the filter settings - in the future we may modify these with
    # overrides in a similar way to the advanced settings
    cp ${filters_filename} ${modified_filters_filename}

    # So that matplotlib doesn't complain
    mkdir -p ./.matplotlib
    export MPLCONFIGDIR=./.matplotlib

    # Ensure that BindCraft finds the correct version of ffmpeg
    export PATH=/opt/conda/envs/BindCraft/bin/:\$PATH

    /opt/conda/envs/BindCraft/bin/python /app/BindCraft/bindcraft.py \
        --settings ${settings_json} \
        --filters /app/BindCraft/settings_filters/default_filters.json \
        --advanced ${modified_advanced_settings_filename} \
        --filters ${modified_filters_filename} \
        ${task.ext.args ?: ''} 2>&1 | tee bindcraft.log

    # Tag per-batch stats CSVs with the input structure filename (incl. extension)
    /opt/conda/envs/BindCraft/bin/python ${baseDir}/bin/bindcraft/add_bindcraft_target_column.py \
        --target "${input_pdb.name}" \
        --results-dir results

    if [[ ${compress_html} == "true" ]]; then
        find ./results -type f -name '*.html' -exec gzip -9 {} +
    fi

    if [[ ${compress_pdb} == "true" ]]; then
        find ./results -type f -name '*.pdb' ! -path './results/Accepted/*.pdb' -exec gzip -9 {} +
    fi

    # Prepare a batch-scoped results directory for downstream reporting staging
    # TODO: Can we use a symlink instead of copying here ?
    mkdir -p "batches/${batch_id}"
    cp -a ./results "batches/${batch_id}/"
    """
}
