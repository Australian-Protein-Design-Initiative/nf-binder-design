process ALPHAFOLD2 {
    tag "${meta.af2_tool ?: 'af2'} ${meta.id} run${meta.af2_run}${meta.msa_depth_tag ? " msa${meta.msa_depth_tag}" : ''}"

    container 'https://bioinformatics.erc.monash.edu/home/andrewperry/containers/ghcr.io-australian-protein-design-initiative-containers-alphafold2-2.3.2-custom.img'

    // Predictions land under "out/${meta.id}/" (the "out/" wrapper avoids
    // colliding with the precomputed-MSA dir that stages in as "${meta.id}").
    // The output is declared as a recursive "out/${meta.id}/**" glob so each
    // nested file is its own publish item - that lets saveAs filter per-file
    // (saveAs/pattern only see top-level output items, never files nested inside
    // a directory item). We strip the task-local "out/" prefix, and drop msas/ +
    // features.pkl, which are published once by ALPHAFOLD2_JACKHMMER_MSA /
    // COLABFOLD_A3M_TO_AF2_MSAS / FOLD_ASSEMBLE_AF2_MULTIMER_MSAS (features.pkl under fold/af2/msas/; raw msas/
    // under fold/msa/<method>/). The "out/" dir carries a copy only because the
    // precomputed MSAs are staged into it.
    // When --n_predictions / --msa_subsample fans AF2 out (meta.af2_namespaced),
    // each run's output is nested under run_<n>[_msa_...]/ so they don't collide.
    // In --af2_keep_models=best mode only the top-ranked model and the ranking
    // table are kept per run; 'all' keeps everything (minus the MSAs/features).
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/${meta.af2_tool ?: 'af2'}/predictions",
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
                if (meta.msa_depth_tag) {
                    runtag = "${runtag}_msa_${meta.msa_depth_tag}"
                }
                return "${meta.id}/${runtag}/${tail}"
            }
            return rel
        }
    )
    // Flat <outdir>/fold/predictions/: relaxed cif when relaxing, unrelaxed
    // (or ranked_0.cif for keep=best) when --af2_no_relax.
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/predictions",
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
            return "${FoldNaming.af2Prefix(meta)}${bn}"
        }
    )
    // Sequence IDs used in each MSA depth job (when --msa_subsample is on).
    publishDir(
        path: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/msa_ids",
        mode: 'copy',
        pattern: '*_ids.txt'
    )

    input:
    tuple val(meta), path(fasta), path(msa_dir), path(a3m)

    output:
    // Recursive glob so each nested file is its own publish item (see the
    // publishDir comment). This also feeds ENGENS and FOLD_SCORE_AF2, which
    // both get every model's files (the glob stages flat into their work dir,
    // so FOLD_SCORE_AF2 reads pkl/pae/pdb with --run-dir .).
    // NB: do NOT also declare the parent dir (path("out/${meta.id}")) as an
    // output - a whole-directory output is published as a single unit, which
    // suppresses the per-file publishing both saveAs blocks below rely on (the
    // flat fold/predictions/ gather then silently publishes nothing).
    tuple val(meta), path("out/${meta.id}/**"), emit: predictions
    path("*_ids.txt"), emit: msa_ids, optional: true

    script:
    // Same DB flags as ALPHAFOLD2_JACKHMMER_MSA - see the comment there for why
    // they're required (and duplicated) on both stages, and why multimer swaps
    // pdb70 for pdb_seqres + uniprot and sets model_preset=multimer.
    def d = params.af2_db_path
    def data_dir = params.af2_data_dir ?: d
    def tool = meta.af2_tool ?: 'af2'
    // 'chainbreak' folds a complex with the MONOMER weights: the chains are
    // concatenated and separated only by an af2_chain_break_offset jump in
    // residue_index (see af2_monomer_features_from_msas.py). So a multi-chain input
    // is deliberately NOT multimer here, and it takes the monomer DB flags (pdb70 /
    // hhsearch) rather than the multimer ones (pdb_seqres + uniprot / hmmsearch).
    def chainbreak = meta.af2_complex_mode == 'chainbreak'
    def is_multimer = !chainbreak && ((meta.n_chains ?: 1) > 1 || params.af2_model_preset == 'multimer')
    def model_preset = is_multimer \
        ? 'multimer' \
        : (chainbreak ? params.af2_monomer_model_preset : params.af2_model_preset)
    def db_flags_list = [
        "--data_dir=${data_dir}",
        "--uniref90_database_path=${d}/uniref90/uniref90.fasta",
        "--mgnify_database_path=${d}/${params.af2_mgnify_subpath}",
        "--bfd_database_path=${d}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
        "--uniref30_database_path=${d}/${params.af2_uniref30_subpath}",
    ]
    if (is_multimer) {
        db_flags_list += [
            "--uniprot_database_path=${d}/${params.af2_uniprot_subpath}",
            "--pdb_seqres_database_path=${d}/${params.af2_pdb_seqres_subpath}",
        ]
    }
    else {
        // hhsearch.py:61 globs "<path>_*" - the monomer-only existence check that the
        // multimer path (hmmsearch) never reaches.
        db_flags_list += ["--pdb70_database_path=${d}/${params.af2_pdb70_subpath}"]
    }
    db_flags_list += [
        "--template_mmcif_dir=${d}/pdb_mmcif/mmcif_files",
        "--obsolete_pdbs_path=${d}/pdb_mmcif/obsolete.dat",
        "--max_template_date=${params.af2_max_template_date}",
        "--db_preset=${params.af2_db_preset}",
        "--model_preset=${model_preset}",
    ]
    def db_flags = db_flags_list.join(' ')

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
    def num_multimer_predictions_flag = is_multimer
        ? "--num_multimer_predictions_per_model=${params.af2_num_predictions_per_model}"
        : ''
    // Shallow --msa_subsample jobs: rebuild features.pkl from subsampled a3m
    // (empty templates). Full / non-subsample: reuse staged features.pkl.
    def do_subsample = meta.msa_max_seq != null
    def write_msa_ids = meta.msa_depth_tag != null
    def msa_ids_file = write_msa_ids \
        ? "af2_${meta.id}_run${meta.af2_run}_msa${meta.msa_depth_tag}_ids.txt" \
        : ''

    // The three chain-break stanzas below are assembled in Groovy rather than emitted
    // as `if [[ "${chainbreak}" == "true" ]]` shell blocks, so that when af2_mono is
    // NOT in play the rendered command is byte-for-byte what it was before this engine
    // existed. Nextflow hashes the rendered script, so any always-present shell guard -
    // even a dead one - would invalidate the -resume cache of every prior af2 run.
    // Adding an opt-in engine must not cost existing users their caches.
    def build_multimer_features = chainbreak ? '' : """        # fold_pulldown assemble used to omit features.pkl (expecting AF2 to
        # rebuild it). This container's predict_structure() loads the pickle
        # and never re-reads msas/, so build it here if the MSA stage did not.
        if [[ ! -f "out/${meta.id}/features.pkl" ]]; then
            python ${projectDir}/bin/fold/af2_multimer_features_from_msas.py \\
                --fasta ${fasta} \\
                --msas-dir "out/${meta.id}"
        fi
"""
    // Anything staged here is multimer-format (paired + block diagonal); the monomer
    // models cannot read it, so chain-break mode always rebuilds features.pkl.
    def chainbreak_features = !chainbreak ? '' : """
    rm -f "out/${meta.id}/features.pkl"
    python ${projectDir}/bin/fold/af2_monomer_features_from_msas.py \\
        --fasta ${fasta} \\
        --msas-dir "out/${meta.id}" \\
        --chain-break-offset ${params.af2_chain_break_offset} \\
        --layout-out chain_layout.json
"""
    // The monomer models emit one chain with our +offset numbering still in it. Split
    // it back into real chains before anything downstream (FOLD_SCORE_AF2 -> ipsae, and
    // every consumer of fold/predictions/) tries to select by chain.
    def split_chainbreak = !chainbreak ? '' : """

    python ${projectDir}/bin/fold/af2_split_chainbreak.py \\
        --layout chain_layout.json \\
        --also-mmcif \\
        "out/${meta.id}"/*.pdb"""
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

    mkdir -p out
    if [[ "${do_subsample}" == "true" ]]; then
        python3 ${projectDir}/bin/fold/subsample_a3m.py \
            --a3m "${a3m}" \
            --max-seq ${meta.msa_max_seq} \
            --max-extra-seq ${meta.msa_max_extra_seq} \
            --seed ${meta.msa_subsample_seed} \
            -o subsampled.a3m \
            --ids-output "${msa_ids_file}"
        mkdir -p "out/${meta.id}/msas"
        cp subsampled.a3m "out/${meta.id}/msas/subsampled.a3m"
        python ${projectDir}/bin/fold/colabfold_a3m_to_af2_msas.py \
            --fasta ${fasta} \
            --a3m subsampled.a3m \
            --output-dir "out/${meta.id}"
    else
        # -L dereferences the symlink-staged MSA dir (stageInMode='symlink') into real,
        # writable files - AlphaFold rewrites features.pkl and writes prediction
        # outputs alongside it, which fails/corrupts -resume caching against a symlink.
        cp -rL "${msa_dir}" "out/${meta.id}"
        if [[ "${write_msa_ids}" == "true" ]]; then
            python3 ${projectDir}/bin/fold/subsample_a3m.py \
                --a3m "${a3m}" \
                --ids-only \
                --ids-output "${msa_ids_file}"
        fi
${build_multimer_features}    fi
${chainbreak_features}
    python /app/alphafold/run_alphafold.py \
        --fasta_paths=${fasta} \
        --output_dir=\$PWD/out \
        --use_precomputed_msas=true \
        --generate_msas_only=false \
        ${use_gpu_relax_flag} \
        ${random_seed_flag} \
        ${models_to_relax_flag} \
        ${num_multimer_predictions_flag} \
        ${db_flags}${split_chainbreak}
    """
}
