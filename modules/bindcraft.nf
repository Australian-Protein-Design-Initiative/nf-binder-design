process BINDCRAFT {
    container 'ghcr.io/australian-protein-design-initiative/containers/bindcraft:05702c4_nv-cuda12'

    publishDir path: "${params.outdir}/bindcraft/batches/${batch_id}", mode: 'copy'
    publishDir path: "${params.outdir}/bindcraft/accepted", pattern: 'results/Accepted/*.pdb', mode: 'copy'

    input:
    path input_pdb
    path settings_json
    val advanced_settings_preset
    val batch_id

    output:
    path "${params.outdir}/Accepted/*.pdb", emit: accepted_pdbs, optional: true
    path "${params.outdir}/Rejected/*.pdb", emit: rejected_pdbs, optional: true
    path "${params.outdir}/Trajectory/Relaxed/*.pdb", emit: relaxed_pdbs, optional: true
    path "${params.outdir}/Trajectory/LowConfidence/*.pdb", emit: low_confidence_pdbs, optional: true
    path "${params.outdir}/Trajectory/Clashing/*.pdb", emit: clashing_pdbs, optional: true
    path "${params.outdir}/final_design_stats.csv", emit: final_stats_csv, optional: true
    path "${params.outdir}/trajectory_stats.csv", emit: trajectory_stats_csv, optional: true
    path "${params.outdir}/mpnn_design_stats.csv", emit: mpnn_design_stats_csv, optional: true
    path "${params.outdir}/failure_csv.csv", emit: failure_csv, optional: true
    path 'results/**', emit: all_results

    script:
    def advanced_settings_filename = advanced_settings_preset ? "/app/BindCraft/settings_advanced/${advanced_settings_preset}.json" : '/app/BindCraft/settings_advanced/default_4stage_multimer.json'
    def modified_advanced_settings_filename = "./${file(advanced_settings_filename).getName()}"
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

    ##
    # We modify the advanced settings to set the `max_trajectories` to the batch size
    # This way Nextflow can run a defined number of trajectories per task
    # We can monitor the accepted_pdbs channel count if we want to stop after
    # a fixed number of Accepted designs or stop if the accept rate is too low
    # (as per normal BindCraft behaviour)
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

    # So that matplotlib doesn't complain
    mkdir -p ./.matplotlib
    export MPLCONFIGDIR=./.matplotlib

    # Ensure that BindCraft finds the correct version of ffmpeg
    export PATH=/opt/conda/envs/BindCraft/bin/:\$PATH

    /opt/conda/envs/BindCraft/bin/python /app/BindCraft/bindcraft.py \
        --settings ${settings_json} \
        --filters /app/BindCraft/settings_filters/default_filters.json \
        --advanced ${modified_advanced_settings_filename} \
        ${task.ext.args ?: ''}
    """
}
