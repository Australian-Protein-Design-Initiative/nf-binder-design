process ALPHAFOLD2 {
    tag "${meta.id}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif'

    // Predictions land under "out/${meta.id}/" (the "out/" wrapper avoids
    // colliding with the precomputed-MSA dir that stages in as "${meta.id}").
    // The output is declared as a recursive "out/${meta.id}/**" glob so each
    // nested file is its own publish item - that lets saveAs filter per-file
    // (saveAs/pattern only see top-level output items, never files nested inside
    // a directory item). We strip the task-local "out/" prefix, and drop msas/ +
    // features.pkl, which are published once by ALPHAFOLD2_JACKHMMER_MSA under
    // af2/msas/ (the "out/" dir carries a copy only because the precomputed MSAs
    // are staged into it).
    publishDir(
        path: "${params.outdir}/af2/predictions",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString().replaceFirst(/^out\//, '')
            (rel.startsWith("${meta.id}/msas/") || rel == "${meta.id}/features.pkl") ? null : rel
        }
    )

    input:
    tuple val(meta), path(fasta), path(msa_dir)

    output:
    tuple val(meta), path("out/${meta.id}/**"), emit: predictions

    script:
    // Same DB flags as ALPHAFOLD2_JACKHMMER_MSA - see the comment there for why
    // they're required (and duplicated) on both stages.
    def d = params.alphafold2_db_path
    def db_flags = [
        "--data_dir=${d}",
        "--uniref90_database_path=${d}/uniref90/uniref90.fasta",
        "--mgnify_database_path=${d}/mgnify/mgy_clusters_2022_05.fa",
        "--bfd_database_path=${d}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        "--uniref30_database_path=${d}/uniref30/UniRef30_2021_03",
        "--pdb70_database_path=${d}/pdb70/pdb70",
        "--template_mmcif_dir=${d}/pdb_mmcif/mmcif_files",
        "--obsolete_pdbs_path=${d}/pdb_mmcif/obsolete.dat",
        "--max_template_date=${params.alphafold2_max_template_date}",
        "--db_preset=${params.alphafold2_db_preset}",
        "--model_preset=${params.alphafold2_model_preset}",
    ].join(' ')

    // Natively-supported CLI knobs only (confirmed via --helpfull; see plans/
    // on-this-new-branch-groovy-moonbeam.md step 0). --num_recycle and
    // model-subset selection have no native flags in this build, so they are
    // deliberately not wired here.
    def random_seed_flag = params.alphafold2_random_seed ? "--random_seed=${params.alphafold2_random_seed}" : ''
    def models_to_relax_flag = "--models_to_relax=${params.alphafold2_models_to_relax}"
    def num_multimer_predictions_flag = params.alphafold2_model_preset == 'multimer'
        ? "--num_multimer_predictions_per_model=${params.alphafold2_num_predictions_per_model}"
        : ''
    """
    if [[ ${params.require_gpu} == "true" ]]; then
        if [[ \$(nvidia-smi -L) =~ "No devices found" ]]; then
            echo "No GPU detected! Failing fast rather than going slow (since --require_gpu=true)"
            exit 1
        fi

        nvidia-smi
    fi

    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    mkdir -p out
    # -L dereferences the symlink-staged MSA dir (stageInMode='symlink') into real,
    # writable files - AlphaFold rewrites features.pkl and writes prediction
    # outputs alongside it, which fails/corrupts -resume caching against a symlink.
    cp -rL "${msa_dir}" "out/${meta.id}"

    python /app/alphafold/run_alphafold.py \
        --fasta_paths=${fasta} \
        --output_dir=\$PWD/out \
        --use_precomputed_msas=true \
        --generate_msas_only=false \
        --use_gpu_relax=true \
        ${random_seed_flag} \
        ${models_to_relax_flag} \
        ${num_multimer_predictions_flag} \
        ${db_flags}
    """
}
