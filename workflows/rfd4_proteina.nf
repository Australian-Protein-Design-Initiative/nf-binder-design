#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
RFD4-Proteina binder design workflow

Usage via main.nf:
  nextflow run main.nf --method rfd4_proteina --rfd4_config binder.yaml
  nextflow run main.nf --method rfd4_proteina --input_pdb target.cif --contigs "50-120, /0, A17-131" --hotspot_res "A56,A115,A123"

RFD4-Proteina co-designs sequence and structure, so there is no inverse-folding step between
generation and scoring: each design comes out of the model as a complex that already has a
binder sequence.
*/

// Default parameters
params.input_pdb = false
params.outdir = 'results'

params.design_name = 'rfd4'
params.contigs = ''
params.hotspot_res = false

// RFD4-Proteina params
params.rfd4_config = false
params.rfd4_n_designs = 1
params.rfd4_batch_size = 1
params.rfd4_designs_per_task = 1
// 0 (no separate dataloader processes) is the default on purpose: with workers, rfd4-generate
// has been seen to hang before the first sampling step, with the GPU idle and the main thread
// blocked on the dataloader. Input featurisation is cheap relative to sampling, so there is
// little to gain from raising it.
params.rfd4_num_workers = 0
params.rfd4_seed = 0
params.rfd4_nsteps = false
params.rfd4_recycle_strategy = false
params.rfd4_inference_experiment = false
// Which metrics run inline, as a comma-separated tag list (eg 'ppi_refolding,hotspot_conditioning').
// Unset means none, and the metrics group is dropped entirely so no refolding oracle is built:
// RFD4-Proteina's own refolding needs an AlphaFold 3 / OpenFold3 environment that is not in the
// container. Use --full_refold_with boltz for refold scoring instead.
params.rfd4_metrics_tags = false
// false runs the trunk eagerly (compile_lengths=null), skipping the torch.compile warm-up
params.rfd4_compile_lengths = true
params.rfd4_partial_trajectory_fraction = false
params.rfd4_inverse_folding = false
params.rfd4_extra_args = ''
params.rfd4_filters = false

// Chain ids in the generated complex. Derived from --contigs in params mode; set these when
// using a condition spec whose target chain is not A.
params.rfd4_target_chain = false
params.rfd4_binder_chain = false

// --- Refold scoring (REFOLD) ---
// Nextflow does not evaluate top-level params statements in an included script, so this
// block lives here rather than in subworkflows/local/refold.nf. Copy it when adopting
// that subworkflow in another workflow.
//
// One engine list, whichever engines it names: each is folded, compared against the
// design and interface-scored the same way. The fold stack's own parameters
// (--n_predictions, --msa_method, --create_target_msa, --use_msa_server, --templates and
// each engine's flags) govern how the refold folds.
params.full_refold_with = ''  // comma-separated engines, eg 'boltz,esmfold2'
params.full_refold_max = false

// The refold's outputs and the fold stack's own outputs share a directory: the fold
// stack reads fold_publish_dir directly in its modules, so the two are kept equal.
params.refold_publish_dir = 'rfd4_proteina/refold'
params.fold_publish_dir = 'rfd4_proteina/refold'

params.refold_target_id = 'target'
params.refold_metrics = 'iptm_max,iptm_mean,ipsae_max,ipsae_mean,plddt_max,plddt_mean,ptm_max,pae_mean,n'
// Several samples per design (--n_predictions > 1) reduce to one row per design. The
// confidence columns are aggregates over the samples (iptm_max and friends). The
// structural columns - RMSD against the design, and the interface metrics - are not:
// they describe one set of coordinates, so they are all read from the single
// best-ranked sample, by --consensus_metric. A statistic per metric would give a row
// that no predicted structure matches, and the interface metrics have no shared
// direction to take a best in anyway. min / mean / max / first still override it.
// refold_sample_scores.tsv carries every sample's own values.
params.refold_structure_stat = 'selected'
// Fold the binder on its own as well, for the "does it hold its fold without the target"
// comparison. Doubles the folding work, so it can be turned off.
params.refold_binder_monomer = true
params.refold_monomer_chain = 'A'
params.refold_monomer_single_sequence = true
params.refold_interface_scoring = true
// Keep the superimposed structures from the RMSD comparisons, for looking at.
params.output_rmsd_aligned = false
params.refold_fasta_score_columns = 'consensus_z'

params.uniref30 = false
params.colabfold_envdb = false

// --- deprecated: --refold_oracles was a second parameter for the same thing ---
params.refold_oracles = ''
params.refold_oracles_max = false
params.refold_oracle_metrics = false
params.refold_oracle_target_id = false
// Refolding takes its MSA and template choices from the fold stack now, so these have no
// equivalent of their own. They are still read, to fail with a message that says so.
params.full_refold_use_msa_server = false
params.full_refold_create_target_msa = false
params.full_refold_target_templates = false
params.full_refold_target_fasta = false
params.full_refold_alignment = false

// --- Shared fold-stack parameters (the fold / fold_pulldown entrypoints declare
// their own copies; everything else these engines read lives in nextflow.config) ---
params.msa_method = 'jackhmmer_af2'
params.create_target_msa = false
params.create_binder_msa = false
// 5 structures per design per engine. Each is a separate job of one sample with its
// own seed, not one job of five samples, because the engines that batch samples share
// a single trunk / recycling pass across them: five samples from one seed vary only in
// the diffusion noise, while five seeds re-run the whole prediction. The per-engine
// batch sizes below are what forces the one-sample-per-job split (see
// subworkflows/local/*_fold.nf foldPredictionBatches). The disadvantage is five model
// loads instead of one, which for Boltz is roughly 45 s of fixed cost per sample.
params.n_predictions = 5
params.use_msa_server = false
params.templates = false
params.msa_subsample = false
params.msa_subsample_include_full = true

// One sample per job and a base seed, which together turn --n_predictions into
// N seeds x 1 sample: the fold subworkflows split n_predictions into jobs of
// batch_size samples and give job i seed base+i. These carry the full_refold_
// prefix because the plain names are declared in nextflow.config for the fold /
// fold_pulldown entrypoints, where a config declaration beats a workflow script's
// default; the fold stack reads the prefixed name wherever it is declared (see
// lib/FoldParams.groovy). Protenix's copies are this workflow's own declarations
// further down, so they take the refold's defaults directly.
params.full_refold_boltz_batch_size = 1
params.full_refold_boltz_seed = 42
params.full_refold_esmfold2_batch_size = 1
params.full_refold_esmfold2_seeds = 42
params.full_refold_chai1_batch_size = 1
params.full_refold_chai1_seeds = 42
params.full_refold_openfold3_batch_size = 1
params.full_refold_openfold3_seeds = 42
params.full_refold_af3_batch_size = 1

// Ranking, as in fold_pulldown (see bin/fold_pulldown_summarise.py)
params.consensus_metric = 'ipsae'
params.z_stat = 'max'
params.z_scope = 'target'
params.min_pool = 10

// --- AF2 ---
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
params.af2_pdb70_subpath = 'pdb70/pdb70'
params.af2_data_dir = '/app/alphafold'
params.af2_keep_models = 'all'
params.af2_no_relax = false
params.af2_publish_pkl = true
params.af2_monomer_model_preset = 'monomer_ptm'
params.af2_chain_break_offset = 200

// --- RF3 ---
// Prefixed, unlike the other engines': rfd3 declares --rf3_batch_size, --rf3_num_steps
// and --rf3_n_recycles for its own design-time RF3 step with different meanings, and a
// workflow's parameters are one flat namespace. The fold stack reads the prefixed names
// wherever they are declared (lib/FoldParams.groovy), so plain --rf3_* has no effect
// here.
params.full_refold_rf3_ckpt_path = '/models/foundry/rf3_foundry_01_24_latest_remapped.ckpt'
params.full_refold_rf3_num_steps = 50
params.full_refold_rf3_n_recycles = 10
params.full_refold_rf3_batch_size = 1
params.full_refold_rf3_early_stopping_plddt_threshold = 0.5
params.full_refold_rf3_seed = 42

// --- Protenix ---
params.protenix_seeds = 42
params.protenix_cycle = 10
params.protenix_step = 200
params.protenix_batch_size = 1
params.protenix_model_name = 'protenix_base_default_v1.0.0'
params.protenix_use_msa = true
params.protenix_need_atom_confidence = true

// --- EnGens off: one report per design x oracle is rarely what a design run wants ---
params.skip_engens = true
params.engens_dimred = 'umap'
params.engens_clustering = 'hdbscan'
params.engens_min_structures = 3
params.engens_max_clusters = 10
params.engens_gmm_ic = 'aic'
params.engens_seed = false
params.engens_superpose_method = 'blosum62'
params.engens_featurizers = 'default,3di'

// FoldSeek params
params.do_foldseek = false
params.foldseek_search_max = false  // max designs for FoldSeek (defaults to full_refold_max)
params.foldseek_database = 'CATH50'
params.foldseek_databases_path = false
params.foldseek_use_webserver = false
params.foldseek_mode = '3diaa'
params.foldseek_maxaccept = 1
params.foldseek_gzip_output = false
params.foldseek_include_html_output = false
params.foldseek_cath_names_path = false

params.require_gpu = true
params.gpu_devices = ''
params.gpu_slots_per_device = 1
params.gpu_lock_timeout = 14400

include { UNIQUE_ID } from '../modules/local/common/unique_id'
include { PARSE_RFD4_CONFIG } from '../modules/local/rfd4_proteina/parse_rfd4_config'
include { GENERATE_RFD4_CONFIG } from '../modules/local/rfd4_proteina/generate_rfd4_config'
include { PREPARE_RFD4_INPUT } from '../modules/local/rfd4_proteina/prepare_rfd4_input'
include { RFD4_GENERATE } from '../modules/local/rfd4_proteina/rfd4_generate'
include { FILTER_DESIGNS as RFD4_FILTER_DESIGNS } from '../modules/local/rfd4_proteina/filter_designs'
include { COMBINE_RFD4_SCORES } from '../modules/local/rfd4_proteina/combine_rfd4_scores'
include { resolveRfd4Chains } from '../modules/local/rfd4_proteina/rfd4_utils'
// Needs the refold / fold-stack params declared above.
include { REFOLD } from '../subworkflows/local/refold'
include { FOLDSEEK_SEARCH } from '../subworkflows/local/foldseek_search'
include { FOLDSEEK_PREPARE_QUERIES } from '../modules/local/foldseek/foldseek_prepare_queries'

workflow RFD4_PROTEINA {

    main:

    if (params.rfd4_config && params.input_pdb) {
        throw new Exception('--rfd4_config and --input_pdb are mutually exclusive. Use --rfd4_config (input paths are read from the spec) or --input_pdb with --contigs.')
    }

    // Refolding is one path with one engine list now, so the parameters that existed
    // only because it was two paths are gone. Each is still read, so a run that sets one
    // stops with the replacement named rather than quietly ignoring it.
    [
        full_refold_use_msa_server: '--use_msa_server',
        full_refold_create_target_msa: '--create_target_msa',
        full_refold_target_templates: '--templates',
        full_refold_alignment: '--create_target_msa with --msa_method mmseqs2_colabfold',
        full_refold_target_fasta: 'a target sequence in the design structures themselves',
    ].each { name, replacement ->
        if (params[name]) {
            error(
                "rfd4_proteina: --${name} has been replaced. Refolding takes its MSA and " +
                "template choices from the fold stack now; use ${replacement}."
            )
        }
    }

    // Show help message
    if (params.input_pdb == false && params.rfd4_config == false) {
        log.info(
            """
        ==================================================================
        PROTEIN BINDER DESIGN PIPELINE - RFD4-Proteina
        ==================================================================

        Required (one of):
            --rfd4_config         Path to an RFD4-Proteina condition spec (YAML or JSON)
                                  Input structure paths are read from the spec.
            --input_pdb           Input CIF/PDB file (params mode; a spec is generated for you)

        Config mode (--rfd4_config):
            The spec is the baseline interface and gives full access to RFD4-Proteina's
            conditioning. See https://rfd4-proteina.org/latest/reference/yaml-json.html

        Params mode (--input_pdb, generates a spec automatically):
            --contigs             Contig string, e.g. "50-120, /0, A17-131"
                                  A length or length range is the designed binder; a chain id
                                  with residue numbers is preserved from the input. [default: ${params.contigs}]
            --hotspot_res         Hotspot residues, eg "A56,A115,A123" [default: ${params.hotspot_res}]

        RFD4-Proteina options:
            --outdir              Output directory [default: ${params.outdir}]
            --design_name         Name of the design [default: ${params.design_name}]
            --rfd4_n_designs      Total number of designs [default: ${params.rfd4_n_designs}]
            --rfd4_designs_per_task  Designs generated per Nextflow task (num_replicates) [default: ${params.rfd4_designs_per_task}]
            --rfd4_batch_size     Designs processed in parallel on the GPU within a task [default: ${params.rfd4_batch_size}]
            --rfd4_num_workers    Dataloader worker processes. Leave at 0; rfd4-generate can hang
                                  before the first sampling step when workers are used. [default: ${params.rfd4_num_workers}]
            --rfd4_seed           Base random seed; each task uses this plus its index [default: ${params.rfd4_seed}]
            --rfd4_nsteps         generation.args.nsteps, the number of integration steps [default: model default, 60]
            --rfd4_recycle_strategy  Trunk recycle schedule, eg thinking_medium, high/midpeak [default: model default]
            --rfd4_inference_experiment  Experiment preset, eg quickstart [default: inference_on_user_inputs]
            --rfd4_metrics_tags   Comma-separated inline metric tags. Unset drops the metrics group;
                                  RFD4-Proteina's refolding metrics need an AlphaFold 3 / OpenFold3
                                  environment that is not in the container. [default: none]
            --rfd4_compile_lengths  Set false to run the trunk eagerly, skipping the torch.compile
                                  warm-up. Worth it for a smoke test, not for a real run. [default: ${params.rfd4_compile_lengths}]
            --rfd4_partial_trajectory_fraction  Partial diffusion: run only this final fraction of
                                  nsteps, refining the input rather than building from noise [default: disabled]
            --rfd4_inverse_folding  Hold the input backbone fixed and resample sequence only [default: ${params.rfd4_inverse_folding}]
            --rfd4_extra_args     Additional Hydra key=value overrides for rfd4-generate [default: ${params.rfd4_extra_args}]
            --rfd4_filters        Semicolon-separated filters for generated designs, eg "rg<25" [default: disabled]
            --rfd4_target_chain   Target chain id in the generated complex [default: derived from --contigs, else A]
            --rfd4_binder_chain   Binder chain id in the generated complex [default: derived from --contigs, else B]

        Refold scoring (--full_refold_with):
            Folds each design's binder sequence against its target with the engines you
            name, compares the result with the design, and scores the rebuilt interface.
            Every engine gives the same columns in combined_scores.tsv, prefixed with the
            engine name: confidence (iptm, ipSAE, plddt, ptm), RMSD against the design
            (target-aligned binder and whole complex), the binder folded on its own, and
            BindCraft interface metrics. A cross-engine consensus_z ranks designs when
            more than one engine is named.

            --full_refold_with     Comma-separated engines: af2, af2_mono, boltz, rf3,
                                   protenix, af3, openfold3, chai1, esmfold2,
                                   esmfold2_fast [default: ${params.full_refold_with}]
            --full_refold_max      Maximum designs to refold [default: ${params.full_refold_max}]
            --refold_metrics       Summary columns carried into combined_scores.tsv
                                   [default: ${params.refold_metrics}]
            --refold_structure_stat  How several samples per design reduce to one RMSD
                                   and one interface score: selected|min|first|mean|max
                                   [default: ${params.refold_structure_stat}, both read
                                   from the best-ranked sample]
            --refold_binder_monomer  Also fold each binder alone, for the monomer-vs-complex
                                   and monomer-vs-design comparisons [default: ${params.refold_binder_monomer}]
            --refold_monomer_single_sequence  Fold that monomer from sequence alone, with no
                                   MSA search [default: ${params.refold_monomer_single_sequence}]
            --refold_interface_scoring  BindCraft interface metrics on the refolded complex
                                   [default: ${params.refold_interface_scoring}]
            --output_rmsd_aligned  Keep the superimposed structures from the RMSD
                                   comparisons [default: ${params.output_rmsd_aligned}]
            --n_predictions        Structures per design per engine, run as that many
                                   seeds of one sample each [default: ${params.n_predictions}]
            --msa_method           jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]
            --create_target_msa    Build an MSA for the target [default: ${params.create_target_msa}]
            --use_msa_server       Let the engine fetch its own MSA [default: ${params.use_msa_server}]
            --consensus_metric / --z_stat / --z_scope / --min_pool
                                   Ranking, as in --method fold_pulldown
            --colabfold_envdb      ColabFold environment database (local MSA) [default: ${params.colabfold_envdb}]
            --uniref30             UniRef30 database (local MSA) [default: ${params.uniref30}]

            Every engine flag from --method fold --help (--af2_*, --boltz_*,
            --protenix_*, --af3_*, --openfold3_*, --chai1_*, --esmfold2_*) applies here
            too. RF3's are the exception: use --full_refold_rf3_* (ckpt_path, num_steps,
            n_recycles, batch_size, early_stopping_plddt_threshold, seed), since plain
            --rf3_* names are taken by the design-time RF3 step in other workflows.

            --refold_oracles is deprecated: it named the same thing as --full_refold_with,
            which now takes an engine list.

        FoldSeek structural search (optional, enabled with --do_foldseek):
            --do_foldseek                 Enable FoldSeek search on the same designs selected for full refolding
            --foldseek_search_max         Maximum designs to search with FoldSeek [default: same as --full_refold_max]
            --foldseek_database           Database to search [default: ${params.foldseek_database}]
            --foldseek_databases_path     Path to local databases directory (auto-downloaded if unset)
            --foldseek_use_webserver      Use FoldSeek web API instead of local search [default: ${params.foldseek_use_webserver}]
            --foldseek_mode               Search mode [default: ${params.foldseek_mode}]
            --foldseek_maxaccept          Max accepted alignments per query [default: ${params.foldseek_maxaccept}]
            --foldseek_gzip_output        Gzip TSV output [default: ${params.foldseek_gzip_output}]
            --foldseek_include_html_output  Include HTML report [default: ${params.foldseek_include_html_output}]
            --foldseek_cath_names_path    Path to local cath-names.txt (auto-downloaded if unset)

        Other options:
            --require_gpu           Fail tasks without a GPU [default: ${params.require_gpu}]
            --gpu_devices           GPU devices to use (comma-separated or 'all') [default: ${params.gpu_devices}]
            --gpu_slots_per_device  Concurrent tasks allowed per GPU [default: ${params.gpu_slots_per_device}]
            --gpu_lock_timeout      Seconds a task waits for a free GPU [default: ${params.gpu_lock_timeout}]

        """.stripIndent()
        )
        exit(1)
    }

    // One engine list. --refold_oracles named the same thing and is accepted for now,
    // folded into --full_refold_with rather than running a second path beside it.
    def refold_engines = FoldValidation.parseMethods(params.full_refold_with ?: '')
    if (params.refold_oracles) {
        log.warn(
            "rfd4_proteina: --refold_oracles is deprecated and will be removed; it is the " +
            "same thing as --full_refold_with, which now takes a list of engines. Adding " +
            "'${params.refold_oracles}' to it."
        )
        refold_engines = (refold_engines + FoldValidation.parseMethods(params.refold_oracles)).unique()
    }
    // A local, because params are read directly by the processes and assigning to one
    // here would not reach them.
    def refold_max = params.full_refold_max
    if (params.refold_oracles_max) {
        log.warn('rfd4_proteina: --refold_oracles_max is deprecated; use --full_refold_max.')
        refold_max = params.full_refold_max ?: params.refold_oracles_max
    }
    [
        refold_oracle_metrics: '--refold_metrics',
        refold_oracle_target_id: '--refold_target_id',
    ].each { name, replacement ->
        if (params[name]) {
            error("rfd4_proteina: --${name} has been renamed; use ${replacement}.")
        }
    }

    if (refold_engines) {
        // Every design pair here is a 2-chain complex and the binder is folded query-only,
        // exactly as in --method fold_pulldown, so pulldown: true applies for the same reasons.
        def (refold_errors, refold_warnings) = FoldValidation.validate(params, refold_engines, [
            hasMultimer: true,
            af2DbPath: params.af2_db_path,
            pulldown: true,
        ])
        refold_warnings.each { log.warn("rfd4_proteina: --full_refold_with: ${it}") }
        if (refold_errors) {
            error("rfd4_proteina: --full_refold_with: ${refold_errors.join('\nrfd4_proteina: --full_refold_with: ')}")
        }
    }

    def designs_per_task = (params.rfd4_designs_per_task as int)
    if (designs_per_task < 1) {
        throw new Exception('--rfd4_designs_per_task must be >= 1')
    }
    if ((params.rfd4_batch_size as int) < 1) {
        throw new Exception('--rfd4_batch_size must be >= 1')
    }
    if ((params.rfd4_n_designs as int) < 1) {
        throw new Exception('--rfd4_n_designs must be >= 1')
    }

    def chain_pair = resolveRfd4Chains(params)
    def rfd4TargetChain = chain_pair[0]
    def rfd4BinderChain = chain_pair[1]

    UNIQUE_ID()
    // UNIQUE_ID and PARSE_RFD4_CONFIG take only value channels (or none), so their outputs
    // are value channels already and need no .first() to be reusable across tasks.
    def ch_unique_id = UNIQUE_ID.out.id_file.map { it.text.trim() }
    def ch_design_name = ch_unique_id.map { uid -> "${params.design_name}_${uid}" }

    // Initialised rather than bare-declared: an uninitialised `def` in a workflow body makes the
    // Nextflow 24.x parser lose track of `Channel` for the rest of the body.
    def ch_config_spec = Channel.empty()
    def ch_input_structures = Channel.empty()
    if (params.rfd4_config) {
        def config_file = file(params.rfd4_config)
        ch_config_spec = Channel.fromPath(params.rfd4_config).first()

        // A spec may be YAML, which the head node cannot parse, so resolve its input paths in a
        // task that has PyYAML (same approach as PARSE_BOLTZGEN_CONFIG).
        PARSE_RFD4_CONFIG(ch_config_spec, Channel.value(config_file.parent.toString()))
        ch_input_structures = PARSE_RFD4_CONFIG.out
            .map { stdout_text ->
                def paths = stdout_text.trim().split('\n').findAll { it.trim() }.unique()
                if (paths.isEmpty()) {
                    throw new Exception("No 'input' path found in --rfd4_config ${params.rfd4_config}")
                }
                paths.collect { p -> file(p.trim()) }
            }
    }
    else {
        if (!params.contigs) {
            throw new Exception('Either --rfd4_config or --contigs must be provided')
        }

        // RFD4-Proteina reads only .cif/.cif.gz/.bcif, so a PDB input is converted first
        PREPARE_RFD4_INPUT(Channel.fromPath(params.input_pdb).first())
        def ch_input_cif = PREPARE_RFD4_INPUT.out.cif

        GENERATE_RFD4_CONFIG(
            params.design_name,
            ch_input_cif,
            params.contigs,
            params.hotspot_res ?: false,
            params.rfd4_partial_trajectory_fraction ?: false,
        )
        ch_config_spec = GENERATE_RFD4_CONFIG.out.config_spec
        ch_input_structures = ch_input_cif.map { f -> [f] }
    }

    // Chunk the requested designs across tasks, giving each its own seed so -resume is
    // deterministic and no two tasks generate the same design.
    def n_tasks = Math.ceil((params.rfd4_n_designs as double) / designs_per_task).toInteger()
    def task_replicates = (0..<n_tasks).collect { index ->
        Math.min(designs_per_task, (params.rfd4_n_designs as int) - index * designs_per_task)
    }
    def ch_batches = Channel.of(0..<n_tasks).flatten()

    RFD4_GENERATE(
        ch_config_spec,
        ch_input_structures,
        ch_design_name,
        ch_batches,
        ch_batches.map { index -> task_replicates[index] },
        params.rfd4_batch_size,
        ch_batches.map { index -> (params.rfd4_seed as int) + index },
        params.rfd4_extra_args,
    )

    def ch_generated = RFD4_GENERATE.out.cifs.flatten()

    def ch_designs = Channel.empty()
    if (params.rfd4_filters) {
        RFD4_FILTER_DESIGNS(
            ch_generated,
            params.rfd4_filters,
            rfd4BinderChain,
            'filter',
        )
        ch_designs = RFD4_FILTER_DESIGNS.out.accepted
    }
    else {
        ch_designs = ch_generated
    }

    def ch_designs_with_meta = ch_designs.map { cif ->
        tuple([id: cif.name.replaceFirst(/\.(cif|cif\.gz|cif\.zst)$/, '')], cif)
    }

    // Refold scoring: fold each design's binder against its target with every engine in
    // --full_refold_with, compare the result with the design, and score the interface.
    def ch_refold_scores = Channel.empty()
    if (refold_engines) {
        REFOLD(
            ch_designs_with_meta,
            rfd4BinderChain,
            rfd4TargetChain,
            refold_max,
            refold_engines,
            params.msa_method,
        )
        ch_refold_scores = REFOLD.out.scores
    }

    if (params.do_foldseek) {
        def foldseek_max = (params.foldseek_search_max ?: refold_max)
        def ch_foldseek_input = foldseek_max
            ? ch_designs.take(foldseek_max as int)
            : ch_designs

        FOLDSEEK_PREPARE_QUERIES(ch_foldseek_input)
        FOLDSEEK_SEARCH(FOLDSEEK_PREPARE_QUERIES.out.design_chains, Channel.of([id: 'rfd4_foldseek']))
    }

    def ch_design_table = RFD4_GENERATE.out.design_table
        .collectFile(name: 'rfd4_designs.tsv', keepHeader: true, skip: 1, sort: false)

    // collect() yields a List; combine() flattens Lists into the parent tuple. Wrap as [list]
    // so one combine element stays a List for path(design_cifs, stageAs: 'cifs/*').
    // ifEmpty takes a value, not a channel: given Channel.of(...) it emits the channel object
    // itself, and staging then fails with "Not a valid path value type". Reached whenever the
    // filters reject every design.
    def ch_design_cifs = ch_designs
        .ifEmpty(file("${projectDir}/assets/dummy_files/empty"))
        .collect()
        .map { cifs -> [cifs as List] }

    def ch_combine_input = ch_design_table
        .combine(ch_refold_scores.ifEmpty(file("${projectDir}/assets/dummy_files/combine_placeholder_refold")))
        .combine(Channel.value(rfd4BinderChain))
        .combine(ch_design_cifs)
        .map { designs, refold_scores, binder_chain, cifs_wrapped ->
            def cifs = (cifs_wrapped instanceof List && cifs_wrapped.size() == 1 && cifs_wrapped[0] instanceof List)
                ? cifs_wrapped[0]
                : cifs_wrapped
            tuple(designs, refold_scores, binder_chain, cifs)
        }

    COMBINE_RFD4_SCORES(ch_combine_input)

    emit:
    rfd4_cifs = RFD4_GENERATE.out.cifs
    design_table = ch_design_table
    combined_scores = COMBINE_RFD4_SCORES.out.combined_scores
    binders_fasta = COMBINE_RFD4_SCORES.out.binders_fasta
    foldseek_tsv = params.do_foldseek ? FOLDSEEK_SEARCH.out.tsv : Channel.empty()
    foldseek_tsv_annotated = params.do_foldseek ? FOLDSEEK_SEARCH.out.tsv_annotated : Channel.empty()
    foldseek_html = params.do_foldseek ? FOLDSEEK_SEARCH.out.html : Channel.empty()
}
