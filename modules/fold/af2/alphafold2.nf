process ALPHAFOLD2 {
    tag "${meta.id} run${meta.af2_run}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/alphafold_cuda12_upstream-c77e5d2_custom-57618c5.sif'

    // Predictions land under "out/${meta.id}/" (the "out/" wrapper avoids
    // colliding with the precomputed-MSA dir that stages in as "${meta.id}").
    // The output is declared as a recursive "out/${meta.id}/**" glob so each
    // nested file is its own publish item - that lets saveAs filter per-file
    // (saveAs/pattern only see top-level output items, never files nested inside
    // a directory item). We strip the task-local "out/" prefix, and drop msas/ +
    // features.pkl, which are published once by ALPHAFOLD2_JACKHMMER_MSA /
    // COLABFOLD_A3M_TO_AF2_MSAS (features.pkl under fold/af2/msas/; raw msas/
    // under fold/msa/<method>/). The "out/" dir carries a copy only because the
    // precomputed MSAs are staged into it.
    // When --n_predictions fans AF2 out into multiple seeded runs (meta.af2_namespaced),
    // each run's output is nested under run_<n>/ so they don't collide.
    // In --af2_keep_models=best mode only the top-ranked model and the ranking
    // table are kept per run; 'all' keeps everything (minus the MSAs/features).
    publishDir(
        path: "${params.outdir}/fold/af2/predictions",
        mode: 'copy',
        saveAs: { filename ->
            def rel = filename.toString().replaceFirst(/^out\//, '')
            if (rel.startsWith("${meta.id}/msas/") || rel == "${meta.id}/features.pkl") {
                return null
            }
            if (meta.af2_keep_models == 'best'
                && rel != "${meta.id}/ranked_0.pdb"
                && rel != "${meta.id}/ranked_0.cif"
                && rel != "${meta.id}/ranking_debug.json") {
                return null
            }
            if (meta.af2_namespaced) {
                def tail = rel.substring(meta.id.toString().length() + 1) // strip "${meta.id}/"
                def runtag = meta.af2_seed != null ? "run_${meta.af2_run}_seed_${meta.af2_seed}" : "run_${meta.af2_run}"
                return "${meta.id}/${runtag}/${tail}"
            }
            return rel
        }
    )
    // Flat <outdir>/fold/predictions/: relaxed cif when relaxing, unrelaxed
    // (or ranked_0.cif for keep=best) when --af2_no_relax.
    publishDir(
        path: "${params.outdir}/fold/predictions",
        mode: 'copy',
        saveAs: { filename ->
            def bn = filename.toString().replaceFirst(/^.*\//, '')
            def keep_best = meta.af2_keep_models == 'best'
            def no_relax = params.af2_no_relax instanceof Boolean \
                ? params.af2_no_relax \
                : params.af2_no_relax.toString().toBoolean()
            def match = false
            if (no_relax) {
                match = keep_best ? (bn == 'ranked_0.cif') : (bn ==~ /unrelaxed_model_.*\.cif/)
            } else {
                match = bn ==~ /relaxed_model_.*\.cif/
            }
            if (!match) { return null }
            return "af2_${meta.id}_run${meta.af2_run}_${bn}"
        }
    )

    input:
    tuple val(meta), path(fasta), path(msa_dir)

    output:
    tuple val(meta), path("out/${meta.id}/**"), emit: predictions

    script:
    // Same DB flags as ALPHAFOLD2_JACKHMMER_MSA - see the comment there for why
    // they're required (and duplicated) on both stages.
    def d = params.af2_db_path
    def db_flags = [
        "--data_dir=${d}",
        "--uniref90_database_path=${d}/uniref90/uniref90.fasta",
        "--mgnify_database_path=${d}/mgnify/mgy_clusters_2022_05.fa",
        "--bfd_database_path=${d}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        "--uniref30_database_path=${d}/uniref30/UniRef30_2021_03",
        "--pdb70_database_path=${d}/pdb70/pdb70",
        "--template_mmcif_dir=${d}/pdb_mmcif/mmcif_files",
        "--obsolete_pdbs_path=${d}/pdb_mmcif/obsolete.dat",
        "--max_template_date=${params.af2_max_template_date}",
        "--db_preset=${params.af2_db_preset}",
        "--model_preset=${params.af2_model_preset}",
    ].join(' ')

    // Natively-supported CLI knobs only (confirmed via --helpfull; see plans/
    // on-this-new-branch-groovy-moonbeam.md step 0). --num_recycle and
    // model-subset selection have no native flags in this build, so they are
    // deliberately not wired here.
    // Seed only pinned when --af2_random_seed was set (the subworkflow fans it
    // out distinctly per run). Unset -> no flag, so AF2 draws its own random
    // seed per run and the command stays stable for -resume.
    def random_seed_flag = meta.af2_seed != null ? "--random_seed=${meta.af2_seed}" : ''
    // --models_to_relax follows --af2_keep_models (best|all) unless --af2_no_relax.
    def no_relax = params.af2_no_relax instanceof Boolean \
        ? params.af2_no_relax \
        : params.af2_no_relax.toString().toBoolean()
    def models_to_relax = no_relax ? 'none' : params.af2_keep_models
    def models_to_relax_flag = "--models_to_relax=${models_to_relax}"
    def use_gpu_relax_flag = "--use_gpu_relax=${no_relax ? 'false' : 'true'}"
    def num_multimer_predictions_flag = params.af2_model_preset == 'multimer'
        ? "--num_multimer_predictions_per_model=${params.af2_num_predictions_per_model}"
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
        ${use_gpu_relax_flag} \
        ${random_seed_flag} \
        ${models_to_relax_flag} \
        ${num_multimer_predictions_flag} \
        ${db_flags}
    """
}
