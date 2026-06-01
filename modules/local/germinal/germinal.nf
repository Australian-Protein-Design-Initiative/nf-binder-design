process GERMINAL {
    tag "${batch_id}"

    // container 'ghcr.io/australian-protein-design-initiative/containers/germinal:018c75f_post-pr67-review'
    container 'ghcr.io/australian-protein-design-initiative/containers/germinal:20260611-104bbdd7'

    publishDir(
        path: "${params.outdir}/germinal/batches/${batch_id}",
        pattern: 'out/results/**',
        mode: 'copy',
        saveAs: { filename ->
            filename.toString().replace('\\', '/')
                .replaceFirst(/^(.*\/)?out\\/results\\//, '')
        }
    )
    publishDir path: "${params.outdir}/germinal/batches/${batch_id}", pattern: 'germinal.log', mode: 'copy'
    publishDir(
        path: "${params.outdir}/germinal/accepted/structures",
        pattern: 'out/results/**/accepted/structures/*.pdb',
        mode: 'copy',
        saveAs: { filename -> java.nio.file.Paths.get(filename).getFileName().toString() }
    )
    publishDir(
        path: "${params.outdir}/germinal/config",
        pattern: 'out/results/**/final_config.yaml',
        mode: 'copy',
        saveAs: { filename -> 'final_config.yaml' },
        enabled: { task.tag == '0' },
    )

    input:
    tuple val(batch_id), val(batch_n_traj)
    path config
    path pdb_dir, stageAs: 'input_pdbs'
    val experiment_name
    val max_passing_designs
    val max_hallucinated_trajectories

    output:
    path "${batch_id}", type: 'dir', followLinks: true, optional: true, emit: batch_dir
    path 'out/results/**/all_trajectories.csv', optional: true, emit: all_trajectories_csv
    path 'out/results/**/failure_counts.csv', optional: true, emit: failure_counts_csv
    path 'out/results/**/accepted/designs.csv', optional: true, emit: accepted_designs_csv
    path 'out/results/**/trajectories/designs.csv', optional: true, emit: trajectories_designs_csv
    path 'out/results/**/redesign_candidates/designs.csv', optional: true, emit: redesign_designs_csv
    path 'out/results/**/accepted/structures/*.pdb', optional: true, emit: accepted_pdbs
    path 'germinal.log', optional: true, emit: germinal_log
    path 'out/results/**', emit: all_results

    script:
    def config_name = config.name.replaceAll(/\.(yaml|yml)$/, '')
    """
    set -euo pipefail

    if [[ ${params.require_gpu} == "true" ]]; then
        if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
            echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi
        nvidia-smi
    fi

    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.germinal_gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    export XLA_PYTHON_CLIENT_PREALLOCATE=false
    export XLA_CLIENT_MEM_FRACTION=0.5
    export PATH="/opt/conda/envs/germinal/bin:/opt/conda/bin:/opt/conda/envs/protenix/bin:\$PATH"
    export MPLCONFIGDIR=./.matplotlib
    mkdir -p ./.matplotlib

    cp -rL /workspace/configs ./configs
    cp -L ${config} "./configs/${config.name}"

    mkdir -p pdbs
    cp -a input_pdbs/. pdbs/
    for scaffold in nb scfv; do
        if [[ ! -f "pdbs/\${scaffold}.pdb" && -f "/workspace/pdbs/\${scaffold}.pdb" ]]; then
            cp "/workspace/pdbs/\${scaffold}.pdb" "pdbs/\${scaffold}.pdb"
        fi
    done

    # Hydra resolves relative --config-path from /workspace (run_germinal.py location), not \$PWD
    CONFIG_DIR="\$(pwd)/configs"
    PDB_DIR="\$(pwd)/pdbs"
    PROJECT_DIR="\$(pwd)/out"

    python -u /workspace/run_germinal.py \\
        --config-path="\${CONFIG_DIR}" \\
        --config-name=${config_name} \\
        project_dir="\${PROJECT_DIR}" \\
        results_dir=results \\
        experiment_name=${experiment_name} \\
        pdb_dir="\${PDB_DIR}" \\
        dssp_path=/workspace/params/dssp \\
        dalphaball_path=/workspace/params/DAlphaBall.gcc \\
        af_params_dir=/workspace/params \\
        max_trajectories=${batch_n_traj} \\
        max_hallucinated_trajectories=${max_hallucinated_trajectories} \\
        max_passing_designs=${max_passing_designs} \\
        ${task.ext.args ?: ''} 2>&1 | tee germinal.log

    mkdir -p "${batch_id}"
    cp -a "\${PROJECT_DIR}/results/${experiment_name}" "${batch_id}/"
    """
}
