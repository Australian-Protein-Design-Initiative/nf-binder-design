process BOLTZ_COMPARE_BINDER_MONOMER {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/boltz:v2.2.1-2'
    publishDir "${params.outdir}/boltz_refold/predict/binder_monomer", pattern: 'boltz_results_*', mode: 'copy'

    input:
    tuple val(meta), path(af2ig_pdb), path(boltz_complex_pdb)
    val binder_chain

    output:
    path ("boltz_results_${meta.id}_monomer"), emit: results
    tuple val(meta), path("boltz_results_${meta.id}_monomer/predictions/${meta.id}_monomer/${meta.id}_monomer_model_0.pdb"), emit: pdb
    tuple val(meta), path("rmsd_monomer_vs_af2ig_${meta.id}.tsv"), emit: rmsd_monomer_vs_af2ig
    tuple val(meta), path("rmsd_monomer_vs_complex_${meta.id}.tsv"), emit: rmsd_monomer_vs_complex
    tuple val(meta), path("confidence_monomer_${meta.id}.tsv"), emit: confidence_tsv

    script:
    def design_id = meta.id
    def binder_id = "${design_id}_binder_monomer"
    def yaml_file = "${design_id}_monomer.yml"

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

    # Step 1: Create Boltz YAML for binder monomer only
    # (no target specified = monomer mode)
    /usr/bin/python3 ${projectDir}/bin/create_boltz_yaml.py \\
        --binder_id '${binder_id}' \\
        --binder_from_pdb '${af2ig_pdb}' \\
        --binder_chains '${binder_chain}' \\
        --binder_msa empty \\
        --output_yaml '${yaml_file}'

    # Step 2: Run Boltz prediction for monomer
    boltz predict \\
        --preprocessing-threads ${task.cpus} \\
        --num_workers ${task.cpus} \\
        --output_format pdb \\
        ${yaml_file}

    # Step 3: Run RMSD calculations
    # Create directories for rmsd4all
    mkdir -p fixed/
    mkdir -p mobile/

    # Get the monomer PDB path
    MONOMER_PDB="boltz_results_${meta.id}_monomer/predictions/${meta.id}_monomer/${meta.id}_monomer_model_0.pdb"

    # Run RMSD: monomer vs AF2IG binder (chain A)
    # Create fresh directories for this comparison
    rm -rf fixed/ mobile/
    mkdir -p fixed/
    mkdir -p mobile/
    
    ln -s "\$(readlink -f ${af2ig_pdb})" "fixed/\$(basename ${af2ig_pdb})"
    ln -s "\$(readlink -f \$MONOMER_PDB)" "mobile/\$(basename \$MONOMER_PDB)"

    /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${binder_chain} \\
        --mobile-superimpose-chains ${binder_chain} \\
        --score-chains ${binder_chain} \\
        --mobile-score-chains ${binder_chain} \\
        fixed/ mobile/ > rmsd_monomer_vs_af2ig_${meta.id}.tsv

    # Run RMSD: monomer vs Boltz complex binder (chain A)
    # Create fresh directories for this comparison
    rm -rf fixed/ mobile/
    mkdir -p fixed/
    mkdir -p mobile/
    
    ln -s "\$(readlink -f ${boltz_complex_pdb})" "fixed/\$(basename ${boltz_complex_pdb})"
    ln -s "\$(readlink -f \$MONOMER_PDB)" "mobile/\$(basename \$MONOMER_PDB)"

    /usr/bin/python3 ${projectDir}/bin/rmsd4all.py \\
        --tm-score \\
        --superimpose-chains ${binder_chain} \\
        --mobile-superimpose-chains ${binder_chain} \\
        --score-chains ${binder_chain} \\
        --mobile-score-chains ${binder_chain} \\
        fixed/ mobile/ > rmsd_monomer_vs_complex_${meta.id}.tsv

    # Step 4: Parse confidence JSON
    /usr/bin/python3 ${projectDir}/bin/parse_boltz_confidence.py \\
        --json "boltz_results_${meta.id}_monomer/predictions/${meta.id}_monomer/confidence_${meta.id}_monomer_model_0.json" \\
        --id "${meta.id}_monomer" \\
        --binder "${binder_id}" > confidence_monomer_${meta.id}.tsv
    """
}
