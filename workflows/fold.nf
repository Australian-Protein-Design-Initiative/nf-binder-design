#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Multi-method structure folding: predicts structures for FASTA inputs with any
combination of --methods af2,boltz,rf3,protenix, sharing a single MSA-generation
stage (FOLD_MSA) with a selectable --msa_method.

Usage via main.nf:
  nextflow run main.nf --method fold --input 'input/*.fasta' --outdir results \
      --methods af2,boltz,rf3,protenix --msa_method jackhmmer_af2 -profile slurm,m3
*/

params.method = 'fold'
params.input = false
params.help = false
params.outdir = 'results'
params.methods = 'af2'
params.msa_method = 'jackhmmer_af2'
params.fold_publish_dir = 'fold'

// Total predicted structures per input, per method. Left UNSET by default so
// each method falls back to its own per-fold default.
params.n_predictions = false

// --- AF2 ---
params.af2_db_path = '/mnt/datasets/alphafold/alphafold_20240229'
params.af2_model_preset = 'monomer_ptm'
params.af2_db_preset = 'full_dbs'
params.af2_max_template_date = '2024-01-01'
params.af2_random_seed = false
params.af2_num_predictions_per_model = 1
params.af2_uniref30_subpath = 'uniref30/UniRef30_2021_03'
params.af2_uniprot_subpath = 'uniprot/uniprot.fasta'
params.af2_pdb_seqres_subpath = 'pdb_seqres/pdb_seqres.txt'
params.af2_mgnify_subpath = 'mgnify/mgy_clusters_2022_05.fa'
// pdb70 is reached only by the monomer presets (--methods af2_mono): multimer
// template search uses hmmsearch over pdb_seqres instead of hhsearch over pdb70.
params.af2_pdb70_subpath = 'pdb70/pdb70'
// See fold_pulldown.nf: alphafold2:2.3.2-custom bundles the model parameters
// and exposes them at /app/alphafold/params, so no host params dir is needed.
params.af2_data_dir = '/app/alphafold'
// --- AF2 monomer chain-break mode (--methods af2_mono) ---
// Fold a complex with the monomer weights: chains concatenated, separated only by a
// jump in residue_index. AF2 clips relative positions at 32, so any offset above that
// reads as "not covalently connected"; 200 is the dl_binder_design convention.
params.af2_monomer_model_preset = 'monomer_ptm'
params.af2_chain_break_offset = 200
params.af2_keep_models = 'best'
params.af2_no_relax = false

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

// --- MSA subsample ---
params.msa_subsample = false
params.msa_subsample_include_full = true

// --- EnGens ---
params.skip_engens = false
params.engens_dimred = 'umap'
params.engens_clustering = 'hdbscan'
params.engens_min_structures = 3
params.engens_max_clusters = 10
params.engens_gmm_ic = 'aic'
params.engens_seed = false
params.engens_superpose_method = 'blosum62'
params.engens_featurizers = 'default,3di'

// --- ColabFold MSA ---
params.use_remote_server = false
params.uniref30 = false
params.colabfold_envdb = false

params.gpu_devices = ''
params.gpu_slots_per_device = 1
params.gpu_lock_timeout = 14400

include { FOLD as FOLD_CORE } from '../subworkflows/local/fold'

workflow FOLD {

    main:

    if (params.help || params.input == false) {
        log.info(
            """
        ==================================================================
        FOLD WORKFLOW (multi-method structure prediction)
        ==================================================================

        Predict structures for one or more FASTA files with any combination of
        AlphaFold2, Boltz-2, RosettaFold3 and Protenix, sharing one
        MSA-generation stage.

        Multimer: a multi-record FASTA folds as a protein complex (one record
        = one chain -> chain IDs A, B, C, ...; homo-oligomers = repeated
        records; up to 26 chains). Header-derived MSA pairing needs
        --msa_method jackhmmer_af2; ColabFold multimer should use
        --use_msa_server. AF2 multimer needs the 2021 DB snapshot
        (see --af2_db_path).

        Required arguments:
            --input                            Single FASTA file, glob, or directory of FASTA files.
                                                Each file is one prediction unit (multi-record = one
                                                complex, folded as chains A, B, C, ...).

        Optional arguments:
            --outdir                           Output directory [default: ${params.outdir}]
            --methods                          Comma-separated list of af2,af2_mono,boltz,rf3,protenix [default: ${params.methods}]
                                                af2      = AlphaFold2-multimer.
                                                af2_mono = AF2 MONOMER weights on a concatenated complex,
                                                           chains separated only by a residue_index jump.
                                                           Shares af2's MSAs; only features.pkl differs.
                                                           See --af2_chain_break_offset. Without an initial
                                                           guess the monomer models often fail to dock at
                                                           all, and only ranking separates the good pose,
                                                           so pair it with --af2_keep_models best.
            --msa_method                       jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]
            --n_predictions                    Total structures per input, per method. Unset (default) => each
                                                engine uses its own default: Boltz/RF3/Protenix emit 5 each,
                                                AF2 does one run keeping per --af2_keep_models. Set N to pin
                                                every diffusion engine to N.
                                                [default: unset]

            AlphaFold2 (--methods includes af2):
            --af2_db_path                       AlphaFold2 database directory [default: ${params.af2_db_path}]
            --af2_model_preset                  monomer|monomer_ptm|monomer_casp14|multimer [default: ${params.af2_model_preset}]
            --af2_db_preset                     full_dbs|reduced_dbs [default: ${params.af2_db_preset}]
            --af2_max_template_date             Maximum template release date [default: ${params.af2_max_template_date}]
            --af2_random_seed                   Fix the data pipeline's random seed [default: unset]
            --af2_keep_models                   Which of AF2's 5 models/run to keep toward --n_predictions
                                                (also which to relax, unless --af2_no_relax):
                                                'all'  = keep 5/run  -> ceil(N/5) runs;
                                                'best' = keep 1/run  -> N runs [default: ${params.af2_keep_models}]
            --af2_no_relax                      Skip Amber relaxation [default: ${params.af2_no_relax}]
            --af2_pdb70_subpath                 pdb70 prefix under --af2_db_path; monomer presets only
                                                [default: ${params.af2_pdb70_subpath}]

            AF2 monomer chain-break (--methods includes af2_mono):
            --af2_monomer_model_preset          monomer|monomer_ptm|monomer_casp14 [default: ${params.af2_monomer_model_preset}]
            --af2_chain_break_offset            residue_index jump at each chain break; must exceed AF2's
                                                relative-position clip of 32 [default: ${params.af2_chain_break_offset}]

            Boltz-2 (--methods includes boltz):
            --use_msa_server                   Use Boltz's own MMseqs2 MSA server [default: ${params.use_msa_server}]
            --templates                        Templates directory with .cif files [default: ${params.templates}]
            --boltz_recycling                  Boltz --recycling_steps override [default: boltz's own default]
            --boltz_batch_size                 Samples per Boltz job (--diffusion_samples)
            --boltz_sampling_steps              Boltz --sampling_steps override [default: boltz's own default]
            --boltz_seed                        Boltz --seed [default: unset]

            RosettaFold3 (--methods includes rf3):
            --rf3_ckpt_path                     RF3 checkpoint path [default: ${params.rf3_ckpt_path}]
            --rf3_num_steps                     [default: ${params.rf3_num_steps}]
            --rf3_n_recycles                    [default: ${params.rf3_n_recycles}]
            --rf3_batch_size                    Samples per RF3 job (diffusion_batch_size)
            --rf3_early_stopping_plddt_threshold [default: ${params.rf3_early_stopping_plddt_threshold}]
            --rf3_seed                          RF3 hydra seed= [default: unset]

            Protenix (--methods includes protenix):
            --protenix_seeds                    Single seed [default: unset]
            --protenix_cycle                    Pairformer cycles [default: ${params.protenix_cycle}]
            --protenix_step                     Diffusion steps [default: ${params.protenix_step}]
            --protenix_batch_size               Samples per Protenix job (--sample)
            --protenix_model_name               Checkpoint name [default: ${params.protenix_model_name}]
            --protenix_use_msa                  Feed shared a3m to Protenix [default: ${params.protenix_use_msa}]
            --protenix_need_atom_confidence     Write full-confidence JSON (PAE matrix) per sample [default: ${params.protenix_need_atom_confidence}]

            MSA subsample:
            --msa_subsample                     false (default), true (CF-random depths), or custom list
            --msa_subsample_include_full        Also keep one full-MSA job [default: ${params.msa_subsample_include_full}]

            ColabFold MSA (--msa_method mmseqs2_colabfold):
            --use_remote_server                 Query the ColabFold MMseqs2 API [default: ${params.use_remote_server}]
            --uniref30                          UniRef30 database path [default: ${params.uniref30}]
            --colabfold_envdb                   ColabFold environment database path [default: ${params.colabfold_envdb}]

            EnGens (runs by default after prediction):
            --skip_engens                       Skip EnGens clustering [default: ${params.skip_engens}]
            --engens_clustering                 hdbscan (default), gmm, km, or combinations
            --engens_dimred                     umap [default: ${params.engens_dimred}]
            --engens_min_structures             [default: ${params.engens_min_structures}]
            --engens_max_clusters               [default: ${params.engens_max_clusters}]
            --engens_gmm_ic                     aic|bic [default: ${params.engens_gmm_ic}]
            --engens_seed                       Optional RNG seed [default: unset]
            --engens_featurizers                default,3di,pb [default: ${params.engens_featurizers}]

            --gpu_devices                        GPU devices [default: ${params.gpu_devices}]
            --gpu_slots_per_device               Concurrent tasks allowed per GPU [default: ${params.gpu_slots_per_device}]
            --gpu_lock_timeout                    Seconds a task waits for a free GPU [default: ${params.gpu_lock_timeout}]

        Example:
            nextflow run main.nf --method fold --input 'input/*.fasta' --outdir results \\
                --methods af2,boltz,rf3,protenix --msa_method jackhmmer_af2 -profile slurm,m3

        """.stripIndent()
        )
        exit(1)
    }

    def methods = FoldValidation.parseMethods(params.methods)

    def p = params.input
    def resolved = file(p).isDirectory() ? file("${p}/*.{fasta,fa,faa}") : file(p)
    List input_paths = (resolved instanceof List) ? resolved : [resolved]

    if (!input_paths) {
        error("fold: no FASTA files found for --input '${p}'")
    }

    def MAX_CHAINS = 26
    def chain_counts = [:]
    input_paths.each { f ->
        def (n_chains, empty_records) = FoldValidation.countFastaChains(f)
        if (n_chains == 0) {
            error("fold: ${f} contains no FASTA records.")
        }
        if (n_chains > MAX_CHAINS) {
            error(
                "fold: ${f} has ${n_chains} FASTA records but at most ${MAX_CHAINS} chains " +
                "(A-Z) are supported. Reduce the chain count."
            )
        }
        if (empty_records > 0) {
            error("fold: ${f} has ${empty_records} FASTA record(s) with an empty sequence.")
        }
        chain_counts[f.toString()] = n_chains
    }

    def has_multimer = chain_counts.values().any { it > 1 }
    def (errors, warnings) = FoldValidation.validate(params, methods, [
        hasMultimer: has_multimer,
        af2DbPath: params.af2_db_path,
    ])
    warnings.each { log.warn("fold: ${it}") }
    if (errors) {
        error("fold: ${errors.join('\nfold: ')}")
    }

    ch_input = Channel.fromList(input_paths).map { f ->
        [[id: f.baseName, n_chains: chain_counts[f.toString()]], f]
    }

    FOLD_CORE(ch_input, methods, params.msa_method)

    emit:
    scores = FOLD_CORE.out.scores
    predictions = FOLD_CORE.out.predictions
}
