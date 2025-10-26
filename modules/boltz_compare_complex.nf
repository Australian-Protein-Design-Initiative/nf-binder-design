process BOLTZ_COMPARE_COMPLEX {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1-2'
    publishDir "${params.outdir}/boltz_refold/predict/complex", pattern: 'boltz_results_*', mode: 'copy'

    input:
    tuple val(meta), path(pdb)
    val binder_chain
    val target_chain
    val create_target_msa
    val use_msa_server
    path target_msa
    path binder_msa
    path templates
    path refold_target_fasta

    output:
    path ("boltz_results_${meta.id}"), emit: results
    tuple val(meta), path("boltz_results_${meta.id}/predictions/${meta.id}/${meta.id}_model_0.pdb"), emit: pdb
    tuple val(meta), path("rmsd_target_aligned_binder_${meta.id}.tsv"), emit: rmsd_target_aligned
    tuple val(meta), path("rmsd_complex_${meta.id}.tsv"), emit: rmsd_complex
    tuple val(meta), path("confidence_${meta.id}.tsv"), emit: confidence_tsv

    script:
    def design_id = meta.id
    def target_id = "${design_id}_target"
    def binder_id = "${design_id}_binder"
    def yaml_file = "${design_id}.yml"

    def use_msa_server_flag = use_msa_server ? '--use_msa_server' : ''
    def templates_flag = templates ? "--templates '${templates}'" : ''

    // Determine MSA flags - swap these because we're swapping target/binder roles
    def target_msa_flag = '--binder_msa empty'
    // Binder MSA (always empty)

    def binder_msa_flag = ''
    // Target MSA (may be created or from server)
    if (create_target_msa && use_msa_server) {
        binder_msa_flag = ''
    }
    else if (create_target_msa && !use_msa_server) {
        binder_msa_flag = ''
    }
    else {
        binder_msa_flag = '--target_msa empty'
    }

    // Determine target sequence source - swap target/binder roles so output has binder=A, target=B
    def target_source_args = ''
    if (refold_target_fasta) {
        target_source_args = "--binder_from_fasta '${refold_target_fasta}' --binder_chains '${target_chain}'"
    }
    else {
        target_source_args = "--binder_from_pdb '${pdb}' --binder_chains '${target_chain}'"
    }

    """
    set -euo pipefail

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

    # Step 1: Create Boltz YAML from PDB
    # Swapped: extract binder (chain ${binder_chain}) as target, target (chain ${target_chain}) as binder
    # This results in output structure with binder=chain A, target=chain B
    /usr/bin/python3 ${projectDir}/bin/create_boltz_yaml.py \\
        --target_id '${binder_id}' \\
        --binder_id '${target_id}' \\
        --target_from_pdb '${pdb}' \\
        --target_chains '${binder_chain}' \\
        ${target_source_args} \\
        ${binder_msa_flag} \\
        ${target_msa_flag} \\
        --output_yaml '${yaml_file}' \\
        ${use_msa_server_flag} \\
        ${templates_flag}

    # Step 2: Run Boltz prediction
    boltz predict \\
        --preprocessing-threads ${task.cpus} \\
        --num_workers ${task.cpus} \\
        --output_format pdb \\
        ${yaml_file}

    # Step 3: Run RMSD calculations
    # Create directories for rmsd4all
    mkdir -p fixed/
    mkdir -p mobile/

    # Symlink PDBs into directories
    ln -s "\$(readlink -f ${pdb})" "fixed/\$(basename ${pdb})"
    ln -s "\$(readlink -f boltz_results_${meta.id}/predictions/${meta.id}/${meta.id}_model_0.pdb)" "mobile/\$(basename boltz_results_${meta.id}/predictions/${meta.id}/${meta.id}_model_0.pdb)"

    # Run target-aligned binder RMSD: superimpose on target (B), score binder (A)
    /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${target_chain} \\
        --mobile-superimpose-chains ${target_chain} \\
        --score-chains ${binder_chain} \\
        --mobile-score-chains ${binder_chain} \\
        fixed/ mobile/ > rmsd_target_aligned_binder_${meta.id}.tsv

    # Run complex RMSD: superimpose and score both chains (A,B)
    /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${binder_chain},${target_chain} \\
        --mobile-superimpose-chains ${binder_chain},${target_chain} \\
        --score-chains ${binder_chain},${target_chain} \\
        --mobile-score-chains ${binder_chain},${target_chain} \\
        fixed/ mobile/ > rmsd_complex_${meta.id}.tsv

    # Step 4: Parse confidence JSON
    /usr/bin/python3 ${projectDir}/bin/parse_boltz_confidence.py \\
        --json "boltz_results_${meta.id}/predictions/${meta.id}/confidence_${meta.id}_model_0.json" \\
        --id "${meta.id}" \\
        --target "${target_id}" \\
        --binder "${binder_id}" > confidence_${meta.id}.tsv
    """
}
