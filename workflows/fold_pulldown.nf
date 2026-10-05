#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Fold pulldown: co-fold every binder against every target with one or more of
AF2 / AF2-monomer (af2_mono) / Boltz-2 / RF3 / Protenix / AlphaFold3 / OpenFold3,
then summarise interface scores.

Usage via main.nf:
  nextflow run main.nf --method fold_pulldown \
      --targets targets.fasta --binders binders.fasta \
      --methods boltz,rf3 --create_target_msa true
*/

params.method = 'fold_pulldown'
params.targets = false
params.binders = false
params.help = false
params.outdir = 'results'
params.methods = 'boltz'
params.msa_method = 'jackhmmer_af2'
params.fold_publish_dir = 'fold_pulldown'

params.create_target_msa = false
params.create_binder_msa = false
params.n_predictions = false

// --- Ranking (see bin/fold_pulldown_summarise.py) ---
params.consensus_metric = 'ipsae'
params.z_stat = 'max'
params.z_scope = 'target'
params.min_pool = 10

// --- AF2 (predict uses meta.n_chains=2 -> multimer; jackhmmer on targets stays monomer) ---
params.af2_db_path = '/mnt/datasets/alphafold/alphafold_20211129'
params.af2_model_preset = 'multimer'
params.af2_db_preset = 'full_dbs'
params.af2_max_template_date = '2024-01-01'
params.af2_random_seed = false
params.af2_num_predictions_per_model = 1
params.af2_uniref30_subpath = 'uniclust30/uniclust30_2018_08/uniclust30_2018_08'
params.af2_uniprot_subpath = 'uniprot/uniprot.fasta'
params.af2_pdb_seqres_subpath = 'pdb_seqres/pdb_seqres.txt'
params.af2_mgnify_subpath = 'mgnify/mgy_clusters_2018_12.fa'
// pdb70 is reached only by the monomer presets (--methods af2_mono): multimer
// template search uses hmmsearch over pdb_seqres instead of hhsearch over pdb70.
params.af2_pdb70_subpath = 'pdb70/pdb70'
// AF2 model parameters. The alphafold2:2.3.2-custom container bundles them at
// /models/alphafold2, exposed as /app/alphafold/params by symlink, and
// alphafold/model/data.py resolves `<data_dir>/params/params_<model>.npz` - so
// the in-container default needs no host params dir. Point this at a host
// AlphaFold download to override. (Falls back to af2_db_path if set false.)
params.af2_data_dir = '/app/alphafold'
params.af2_keep_models = 'best'
params.af2_no_relax = false
// AF2's result_model_*.pkl carry the full model output (distogram, MSA and
// structure-module tensors) at ~90 MB each -- 860 MB per prediction, and by far
// the largest thing the run writes. ptm/iptm/ranking_confidence appear ONLY in
// there and nowhere in the published JSON, so they are published by default;
// set false when the merged scores TSV is enough and disk is the binding
// constraint. The pickles still exist in the Nextflow work dir either way.
params.af2_publish_pkl = true
// --- AF2 monomer chain-break mode (--methods af2_mono) ---
// Fold a complex with the monomer weights: chains concatenated, separated only by a
// jump in residue_index. AF2 clips relative positions at 32, so any offset above that
// reads as "not covalently connected"; 200 is the dl_binder_design convention.
params.af2_monomer_model_preset = 'monomer_ptm'
params.af2_chain_break_offset = 200

params.colabfold_msa_publish_name = 'result'

// --- Boltz ---
params.use_msa_server = false
params.templates = false
params.boltz_recycling = false
params.boltz_batch_size = false
params.boltz_sampling_steps = false
params.boltz_seed = false

// --- RF3 ---
params.rf3_ckpt_path = '/models/foundry/rf3_foundry_01_24_latest_remapped.ckpt'
params.rf3_num_steps = 50
params.rf3_n_recycles = 10
params.rf3_batch_size = false
params.rf3_early_stopping_plddt_threshold = 0.5
params.rf3_seed = false

// --- Protenix ---
params.protenix_seeds = false
params.protenix_cycle = 10
params.protenix_step = 200
params.protenix_batch_size = false
params.protenix_model_name = 'protenix_base_default_v1.0.0'
params.protenix_use_msa = true
params.protenix_need_atom_confidence = true

// --- MSA subsample (not supported for pulldown multimers) ---
params.msa_subsample = false
params.msa_subsample_include_full = true

// EnGens off by default (would emit one report per complex = N x M)
params.skip_engens = true
params.engens_dimred = 'umap'
params.engens_clustering = 'hdbscan'
params.engens_min_structures = 3
params.engens_max_clusters = 10
params.engens_gmm_ic = 'aic'
params.engens_seed = false
params.engens_superpose_method = 'blosum62'
params.engens_featurizers = 'default,3di'

params.use_remote_server = false
params.uniref30 = false
params.colabfold_envdb = false

params.gpu_devices = ''
params.gpu_slots_per_device = 1
params.gpu_lock_timeout = 14400

include { FOLD_PULLDOWN_MSA } from '../subworkflows/local/fold_pulldown_msa'
include { FOLD_PREDICT } from '../subworkflows/local/fold_predict'
include { FOLD_PULLDOWN_MERGE_SCORES } from '../modules/local/fold/common/fold_pulldown_merge_scores'
include { FOLD_PULLDOWN_REPORTING } from '../modules/local/common/fold_pulldown_reporting'

workflow FOLD_PULLDOWN {

    main:

    if (params.help || params.targets == false || params.binders == false) {
        log.info(
            """
        ==================================================================
        FOLD PULLDOWN PIPELINE
        ==================================================================

        Co-fold every binder against every target with one or more structure
        predictors (AF2, Boltz-2, RF3, Protenix, AlphaFold3, OpenFold3), then summarise interface
        scores (iptm, ipsae) across replicates and models.

        Required arguments:
            --targets             FASTA file of target sequences
            --binders             FASTA file of binder sequences

        Optional arguments:
            --outdir              Output directory [default: ${params.outdir}]
            --methods             Comma-separated af2,af2_mono,boltz,rf3,protenix,af3,openfold3,esmfold2,esmfold2_fast [default: ${params.methods}]
                                   af2      = AlphaFold2 multimer (every pair here is a 2-chain complex,
                                              so af2 always uses AF2's native multimer weights/pipeline).
                                   af2_mono = AF2 MONOMER weights on a concatenated complex, chains
                                              separated only by a residue_index jump. Shares af2's
                                              MSAs; only features.pkl differs. Without an initial
                                              guess the monomer models often fail to dock at all and
                                              only ranking separates the good pose, so pair it with
                                              --af2_keep_models best. Not an independent engine:
                                              it shares weights lineage with af2.
            --af2_chain_break_offset  residue_index jump per chain break, must exceed AF2's
                                   relative-position clip of 32 [default: ${params.af2_chain_break_offset}]
            --msa_method          jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]
            --create_target_msa   Build MSA for each target [default: ${params.create_target_msa}]
            --create_binder_msa   Build MSA for each binder [default: ${params.create_binder_msa}]
            --n_predictions       Structures per complex per method [default: unset -> engine defaults]
            --use_msa_server      Boltz fetches its own MSA [default: ${params.use_msa_server}]
            --templates           Templates directory with .cif files (reaches Boltz too) [default: ${params.templates}]
            --skip_engens         Skip EnGens clustering [default: ${params.skip_engens}]

            Every method-specific flag from --method fold --help (--af2_*, --boltz_*,
            --rf3_*, --protenix_*, --af3_*, --openfold3_*, --esmfold2_*) applies here too, e.g.
            --af2_keep_models (which of AF2's 5 models/run to keep toward
            --n_predictions) [default: ${params.af2_keep_models}].

            Ranking (summary table):
            --consensus_metric    ipsae|iptm; metric averaged into consensus_z [default: ${params.consensus_metric}]
            --z_stat              max|mean; per-complex statistic over samples that is
                                   standardised [default: ${params.z_stat}]
            --z_scope             target|global; standardise within (target, tool) or over all
                                   targets together [default: ${params.z_scope}]
            --min_pool            flag z-score pools smaller than this as z_pool_small, since
                                   over k complexes |z| cannot exceed (k-1)/sqrt(k)
                                   [default: ${params.min_pool}]

            AF2 (--methods includes af2 or af2_mono) needs the 2021 DB snapshot with uniprot/:
            --af2_db_path         [default: ${params.af2_db_path}]
            --af2_uniref30_subpath / --af2_mgnify_subpath / --af2_uniprot_subpath /
            --af2_pdb_seqres_subpath / --af2_pdb70_subpath   Per-DB subpaths under
                                   --af2_db_path, as for --method fold (see --method fold --help)
            --af2_data_dir        Directory containing params/ (model weights); bundled in the
                                   container by default [default: ${params.af2_data_dir}]

            AlphaFold3 (--methods includes af3; weights are NOT bundled - see models/download_af3_weights.sh):
            --af3_model_dir       Directory holding exactly one af3.bin.zst / af3.bin
                                   [default: ${params.af3_model_dir}]
            --af3_batch_size / --af3_seeds / --af3_num_recycles / --af3_flash_attention /
            --af3_jax_cache_dir   As for --method fold (see --method fold --help)

            OpenFold3 (--methods includes openfold3; weights are bundled in the container):
            --openfold3_batch_size / --openfold3_seeds /
            --openfold3_kernel_cache_dir   As for --method fold (see --method fold --help)

            ESMFold2 (--methods includes esmfold2 and/or esmfold2_fast; weights are bundled in the containers):
            --esmfold2_weights_dir / --esmfold2_batch_size /
            --esmfold2_seeds / --esmfold2_single_sequence / --esmfold2_num_loops /
            --esmfold2_num_sampling_steps / --esmfold2_msa_max_depth /
            --esmfold2_kernel_backend      As for --method fold (see --method fold --help)

            ColabFold MSA (--msa_method mmseqs2_colabfold):
            --use_remote_server   Query ColabFold API [default: ${params.use_remote_server}]
            --uniref30 / --colabfold_envdb   Local ColabFold DB paths

            --gpu_devices           GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
            --gpu_slots_per_device  Concurrent tasks allowed per GPU [default: ${params.gpu_slots_per_device}]
            --gpu_lock_timeout      Seconds a task waits for a free GPU [default: ${params.gpu_lock_timeout}]

        Example:
            nextflow run main.nf --method fold_pulldown \\
                --targets targets.fasta --binders binders.fasta \\
                --methods boltz,rf3,protenix --create_target_msa true \\
                --msa_method jackhmmer_af2 -profile slurm,m3

        """.stripIndent()
        )
        exit(1)
    }

    def methods = FoldValidation.parseMethods(params.methods)
    // AF2 in fold_pulldown assembles the target MSA itself (jackhmmer dir or
    // ColabFold/mmseqs2 a3m via FOLD_ASSEMBLE_AF2_MULTIMER_MSAS) and always
    // folds the binder chain query-only, so neither FoldValidation's generic
    // "AF2 multimer requires jackhmmer_af2" error nor its ColabFold "will run
    // UNPAIRED" warning apply here - opts.pulldown: true skips both.
    def (errors, warnings) = FoldValidation.validate(params, methods, [
        hasMultimer: true,
        af2DbPath: params.af2_db_path,
        pulldown: true,
    ])

    warnings.each { log.warn("fold_pulldown: ${it}") }
    if (errors) {
        error("fold_pulldown: ${errors.join('\nfold_pulldown: ')}")
    }

    if (MsaSubsample.isEnabled(params.msa_subsample)) {
        error("fold_pulldown: --msa_subsample is not supported (multimer pairs only).")
    }

    def target_ids = FoldIds.extractFastaIds(file(params.targets))
    def binder_ids = FoldIds.extractFastaIds(file(params.binders))
    def id_errors = FoldIds.validatePulldownIds(target_ids, binder_ids)
    if (id_errors) {
        error("fold_pulldown: ${id_errors.join('\nfold_pulldown: ')}")
    }

    ch_targets_meta = Channel.fromPath(params.targets)
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }

    ch_binders_meta = Channel.fromPath(params.binders)
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }

    ch_targets_fasta_paths = Channel.fromPath(params.targets).splitFasta(file: true)
    ch_binders_fasta_paths = Channel.fromPath(params.binders).splitFasta(file: true)

    ch_targets = ch_targets_meta.merge(ch_targets_fasta_paths)
    ch_binders = ch_binders_meta.merge(ch_binders_fasta_paths)

    FOLD_PULLDOWN_MSA(ch_targets, ch_binders, methods, params.msa_method)

    ch_pairs_tsv = FOLD_PULLDOWN_MSA.out.pairs_rows.collectFile(
        name: 'pairs.tsv',
        storeDir: "${params.outdir}/${params.fold_publish_dir ?: 'fold_pulldown'}",
        newLine: false,
        seed: "id\ttarget\tbinder\n",
    )

    FOLD_PREDICT(
        FOLD_PULLDOWN_MSA.out.af2_in,
        FOLD_PULLDOWN_MSA.out.for_boltz,
        FOLD_PULLDOWN_MSA.out.for_rf3,
        FOLD_PULLDOWN_MSA.out.for_protenix,
        FOLD_PULLDOWN_MSA.out.for_af3,
        FOLD_PULLDOWN_MSA.out.for_openfold3,
        FOLD_PULLDOWN_MSA.out.for_esmfold2,
        FOLD_PULLDOWN_MSA.out.for_esmfold2_fast,
        methods,
    )

    FOLD_PULLDOWN_MERGE_SCORES(FOLD_PREDICT.out.scores, ch_pairs_tsv)

    FOLD_PULLDOWN_REPORTING(
        file("${projectDir}/assets/fold_pulldown_reporting.qmd"),
        FOLD_PULLDOWN_MERGE_SCORES.out.scores,
        FOLD_PULLDOWN_MERGE_SCORES.out.summary,
    )

    emit:
    scores = FOLD_PULLDOWN_MERGE_SCORES.out.scores
    summary = FOLD_PULLDOWN_MERGE_SCORES.out.summary
}
