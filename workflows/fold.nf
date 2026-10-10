#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Multi-method structure folding: predicts structures for FASTA inputs with any
combination of --methods af2,af2_mono,boltz,rf3,protenix,af3,openfold3,esmfold2,esmfold2_fast, sharing
a single MSA-generation stage (FOLD_MSA) with a selectable --msa_method.

Usage via main.nf:
  nextflow run main.nf --method fold --input 'input/*.fasta' --outdir results \
      --methods af2,boltz,rf3,protenix,af3,openfold3,esmfold2 --msa_method jackhmmer_af2 -profile slurm,m3
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
// AF2's result_model_*.pkl carry the full model output (distogram, MSA and
// structure-module tensors) at ~90 MB each -- 860 MB per prediction, and by far
// the largest thing the run writes. ptm/iptm/ranking_confidence appear ONLY in
// there and nowhere in the published JSON, so they are published by default;
// set false when the merged scores TSV is enough and disk is the binding
// constraint. The pickles still exist in the Nextflow work dir either way.
params.af2_publish_pkl = true

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
        AlphaFold2, Boltz-2, RosettaFold3, Protenix, AlphaFold3 and OpenFold3,
        sharing one MSA-generation stage.

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
            --methods                          Comma-separated list of af2,af2_mono,boltz,rf3,protenix,af3,openfold3,esmfold2,esmfold2_fast [default: ${params.methods}]
                                                af2      = AlphaFold2. Monomer inputs use --af2_model_preset
                                                           (monomer_ptm by default); multi-chain inputs use
                                                           AF2's native multimer weights/pipeline instead.
                                                af2_mono = AF2 MONOMER weights on a concatenated complex,
                                                           chains separated only by a residue_index jump.
                                                           Shares af2's MSAs; only features.pkl differs.
                                                           See --af2_chain_break_offset. Without an initial
                                                           guess the monomer models often fail to dock at
                                                           all, and only ranking separates the good pose,
                                                           so pair it with --af2_keep_models best.
            --msa_method                       jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]
            --n_predictions                    Total structures per input, per method. Unset (default) => each
                                                engine uses its own default: Boltz/RF3/Protenix/AF3/OpenFold3 emit 5 each,
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
            --af2_publish_pkl                   Publish AF2's result_model_*.pkl (~90 MB each; the only
                                                source of ptm/iptm/ranking_confidence)
                                                [default: ${params.af2_publish_pkl}]
            --af2_num_predictions_per_model      AF2's --num_multimer_predictions_per_model; multimer only
                                                [default: ${params.af2_num_predictions_per_model}]
            --af2_data_dir                       Directory containing params/ (model weights); bundled in
                                                the container by default, so this rarely needs overriding
                                                [default: ${params.af2_data_dir}]
            Per-DB subpaths under --af2_db_path (override individually for a non-default
            snapshot, e.g. the 2021 multimer snapshot):
            --af2_uniref30_subpath                [default: ${params.af2_uniref30_subpath}]
            --af2_mgnify_subpath                  [default: ${params.af2_mgnify_subpath}]
            --af2_uniprot_subpath                 multimer only [default: ${params.af2_uniprot_subpath}]
            --af2_pdb_seqres_subpath               multimer only [default: ${params.af2_pdb_seqres_subpath}]
            --af2_pdb70_subpath                  monomer presets only [default: ${params.af2_pdb70_subpath}]

            AF2 monomer chain-break (--methods includes af2_mono):
            --af2_monomer_model_preset          monomer|monomer_ptm|monomer_casp14 [default: ${params.af2_monomer_model_preset}]
            --af2_chain_break_offset            residue_index jump at each chain break; must exceed AF2's
                                                relative-position clip of 32 [default: ${params.af2_chain_break_offset}]

            Templates (every engine except esmfold2 / esmfold2_fast):
            --templates                        Template structures (.pdb/.cif; dir or glob), matched to chains
                                               by sequence alignment; not used by esmfold2 [default: ${params.templates}]
            --template_min_identity            Min identity over aligned residues [default: ${params.template_min_identity}]
            --template_min_coverage            Min fraction of the chain covered [default: ${params.template_min_coverage}]
            --template_min_aligned             Min aligned residues (or the full chain, if shorter) [default: ${params.template_min_aligned}]
            --template_max_per_chain           Templates kept per chain [default: ${params.template_max_per_chain}]
            --boltz_template_force             Hold templated chains near the template (Boltz force) [default: false]
            --boltz_template_threshold         Boltz force threshold in Angstrom [default: ${params.boltz_template_threshold}]

            Boltz-2 (--methods includes boltz):
            --use_msa_server                   Use Boltz's own MMseqs2 MSA server [default: ${params.use_msa_server}]
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
            --protenix_seeds                    Base seed; batch i uses seed+i [default: unset]
            --protenix_cycle                    Pairformer cycles [default: ${params.protenix_cycle}]
            --protenix_step                     Diffusion steps [default: ${params.protenix_step}]
            --protenix_batch_size               Samples per Protenix job (--sample)
            --protenix_model_name               Checkpoint name [default: ${params.protenix_model_name}]
            --protenix_use_msa                  Feed shared a3m to Protenix [default: ${params.protenix_use_msa}]
            --protenix_need_atom_confidence     Write full-confidence JSON (PAE matrix) per sample [default: ${params.protenix_need_atom_confidence}]

            AlphaFold3 (--methods includes af3; weights are NOT bundled - see models/download_af3_weights.sh):
            --af3_model_dir                     Directory holding exactly one af3.bin.zst / af3.bin
                                                [default: ${params.af3_model_dir}]
            --af3_batch_size                    Samples per AF3 job (--num_diffusion_samples)
            --af3_seeds                         Base model seed; batch i uses seed+i [default: 1].
                                                A COMMA LIST ('1,2,3') instead puts every seed in one
                                                job's modelSeeds, as Germinal does; --n_predictions is
                                                then ignored and the run makes
                                                n_seeds x --af3_batch_size structures per complex.
            --af3_num_recycles                  [default: ${params.af3_num_recycles}]
            --af3_flash_attention               auto|triton|cudnn|xla; auto picks xla (plus the XLA
                                                workaround) on pre-Ampere GPUs [default: ${params.af3_flash_attention}]
            --af3_jax_cache_dir                 Persistent JAX compilation cache dir [default: unset]

            AlphaFold3 input and run mode (-profile af3_germinal_parity sets the
            combination Germinal uses: pairing off, template search on, data pipeline on):
            --af3_paired_msa                    false emits "pairedMsa": "" per chain, i.e. no
                                                cross-chain pairing [default: ${params.af3_paired_msa}]
            --af3_templates                     inline|none|search. inline embeds templates matched by
                                                --templates; search omits the key so AF3's data pipeline
                                                searches pdb_seqres/mmcif_files on every chain
                                                [default: ${params.af3_templates}]
            --af3_run_data_pipeline             Run AF3's data pipeline rather than inference only.
                                                MSAs still come from the JSON; this is what enables
                                                AF3's template search [default: ${params.af3_run_data_pipeline}]
            --af3_db_dir                        AlphaFold3 public databases (~630 GB), required when the
                                                data pipeline runs. AF3 validates all nine default paths
                                                before reading the input, so the full set must be present
                                                even though only the template databases are used.
                                                On M3: /mnt/datasets/alphafold3/3.0.0 [default: unset]

            OpenFold3 (--methods includes openfold3; weights are bundled in the container):
            --openfold3_batch_size              Samples per OpenFold3 job (--num-diffusion-samples)
            --openfold3_seeds                   Base model seed; batch i uses seed+i [default: 42]
            --openfold3_kernel_cache_dir        Persistent Triton kernel cache dir [default: unset]

            ESMFold2 (--methods includes esmfold2 and/or esmfold2_fast; weights are bundled in the containers).
            esmfold2 is MSA-conditioned (biohub/ESMFold2); esmfold2_fast (biohub/ESMFold2-Fast) always folds
            from sequence alone. The options below apply to both unless noted.
            --esmfold2_weights_dir              External HF cache dir (HF_HOME) instead of the in-image weights [default: unset]
            --esmfold2_batch_size               Samples per ESMFold2 job (num_diffusion_samples)
            --esmfold2_seeds                    Base model seed; batch i uses seed+i [default: 42]
            --esmfold2_single_sequence          Fold esmfold2 from sequence alone, no MSAs [default: ${params.esmfold2_single_sequence}]
            --esmfold2_num_loops                Trunk loops [default: esm's own, 20]
            --esmfold2_num_sampling_steps       Diffusion steps [default: esm's own, 200]
            --esmfold2_msa_max_depth            esmfold2 MSA rows kept per loop [default: esm's own, 1024]
            --esmfold2_kernel_backend           fused (default), cuequivariance or none

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
            --engens_superpose_method            Superposition scheme for geometric featurizers [default: ${params.engens_superpose_method}]

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
    // A glob makes file() return a List, which has no isDirectory().
    def given = file(p)
    def resolved = (!(given instanceof List) && given.isDirectory()) ? file("${p}/*.{fasta,fa,faa}") : given
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

    def id_errors = FoldIds.validateFoldIds(input_paths.collect { it.baseName })
    if (id_errors) {
        error("fold: ${id_errors.join('\nfold: ')}")
    }

    ch_input = Channel.fromList(input_paths).map { f ->
        [[id: f.baseName, n_chains: chain_counts[f.toString()]], f]
    }

    FOLD_CORE(ch_input, methods, params.msa_method)

    emit:
    scores = FOLD_CORE.out.scores
    predictions = FOLD_CORE.out.predictions
}
