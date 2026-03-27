process BOLTZ_COMPARE_COMPLEX {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1-2'
    publishDir "${params.outdir}/boltz_refold/predict/complex", pattern: 'boltz_results_*', mode: 'copy'
    publishDir "${params.outdir}/boltz_refold/rmsd/aligned_rmsd_target_aligned_binder", pattern: 'aligned_rmsd_target_aligned_binder/*.pdb', mode: 'copy', saveAs: { file(it).name }
    publishDir "${params.outdir}/boltz_refold/rmsd/aligned_rmsd_complex", pattern: 'aligned_rmsd_complex/*.pdb', mode: 'copy', saveAs: { file(it).name }

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
    path ("aligned_rmsd_target_aligned_binder/*.pdb"), emit: aligned_pdbs_target_aligned_binder, optional: true
    path ("aligned_rmsd_complex/*.pdb"), emit: aligned_pdbs_complex, optional: true
    tuple val(meta), path("confidence_${meta.id}.tsv"), emit: confidence_tsv
    tuple val(meta), path("boltz_results_${meta.id}/predictions/${meta.id}/*_ipsae.tsv"), emit: ipsae_tsv
    tuple val(meta), path("boltz_results_${meta.id}/predictions/${meta.id}/*_ipsae_byres.tsv"), emit: ipsae_byres_tsv

    script:
    def design_id = meta.id
    def target_id = "${design_id}_target"
    def binder_id = "${design_id}_binder"
    def yaml_file = "${design_id}.yml"

    def use_msa_server_flag = use_msa_server ? '--use_msa_server' : ''
    def templates_flag = templates ? "--templates '${templates}'" : ''

    def output_transformed_flag_target = params.output_rmsd_aligned ? "--output-transformed aligned_rmsd_target_aligned_binder/" : ''
    def output_transformed_flag_complex = params.output_rmsd_aligned ? "--output-transformed aligned_rmsd_complex/" : ''

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

    // Target sequence for YAML "binder" slot (Boltz output chain B); see create_boltz_yaml.py protein ids A/B
    def target_source_args = ''
    if (refold_target_fasta.name != 'empty') {
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

    # Create Boltz YAML from PDB: binder sequence goes to YAML target slot (Boltz chain A), target sequence to YAML binder slot (Boltz chain B)
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

    # Run Boltz prediction; fail if log shows GPU OOM batch skip (boltz may still exit 0)
    BOLTZ_PREDICT_LOG=.boltz_predict_console.log
    rm -f "\$BOLTZ_PREDICT_LOG"
    set +e
    boltz predict \\
        --preprocessing-threads ${task.cpus} \\
        --num_workers ${task.cpus} \\
        --output_format pdb \\
        ${use_msa_server_flag} \\
        ${yaml_file} 2>&1 | tee "\$BOLTZ_PREDICT_LOG"
    boltz_rc=\${PIPESTATUS[0]}
    set -e
    if grep -qF 'ran out of memory, skipping batch' "\$BOLTZ_PREDICT_LOG"; then
        echo 'BOLTZ_COMPARE_COMPLEX: Boltz logged GPU OOM (batch skipped); failing.' >&2
        exit 1
    fi
    if [[ "\$boltz_rc" -ne 0 ]]; then
        exit "\$boltz_rc"
    fi

    # Run RMSD calculations
    # Create directories for rmsd4all
    mkdir -p fixed/
    mkdir -p mobile/

    # Symlink PDBs into directories
    ln -s "\$(readlink -f ${pdb})" "fixed/\$(basename ${pdb})"
    ln -s "\$(readlink -f boltz_results_${meta.id}/predictions/${meta.id}/${meta.id}_model_0.pdb)" "mobile/\$(basename boltz_results_${meta.id}/predictions/${meta.id}/${meta.id}_model_0.pdb)"

    # Boltz complex PDB chain IDs are always A,B from create_boltz_yaml.py: YAML "target" slot is
    # chain A (here: input binder sequence), YAML "binder" slot is chain B (input target sequence).
    # Input design structure keeps original chain IDs (${target_chain}=target, ${binder_chain}=binder).
    # Mobile must use A=binder and B=target, not the input chain letters.
    # Target-aligned binder RMSD: superimpose on target, score binder
    /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${target_chain} \\
        --mobile-superimpose-chains B \\
        --score-chains ${binder_chain} \\
        --mobile-score-chains A \\
        ${output_transformed_flag_target} \\
        fixed/ mobile/ > rmsd_target_aligned_binder_${meta.id}.tsv

    # Complex RMSD: Boltz mobile is A=binder, B=target. Fixed side keeps input chain IDs; when those
    # match Boltz (e.g. rfd AF2ig A=binder B=target), order aligns. Otherwise sequence-based alignment
    # still pairs homologous chains; target-aligned row above is the robust pose metric.
    /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${binder_chain},${target_chain} \\
        --mobile-superimpose-chains A,B \\
        --score-chains ${binder_chain},${target_chain} \\
        --mobile-score-chains A,B \\
        ${output_transformed_flag_complex} \\
        fixed/ mobile/ > rmsd_complex_${meta.id}.tsv

    /usr/bin/python3 ${projectDir}/bin/ipsae.py \\
        --update-summary "boltz_results_${meta.id}/predictions/${meta.id}/confidence_${meta.id}_model_0.json" \\
        --binder-chain A \\
        --target-chain B \\
        --format boltz \\
        boltz_results_${meta.id}/predictions/${meta.id}/pae*.npz \\
        boltz_results_${meta.id}/predictions/${meta.id}/*.pdb \\
        10 10

    # Parse confidence JSON
    /usr/bin/python3 ${projectDir}/bin/parse_boltz_confidence.py \\
        --json "boltz_results_${meta.id}/predictions/${meta.id}/confidence_${meta.id}_model_0.json" \\
        --id "${meta.id}" \\
        --target "${target_id}" \\
        --binder "${binder_id}" > confidence_${meta.id}.tsv
    """
}
