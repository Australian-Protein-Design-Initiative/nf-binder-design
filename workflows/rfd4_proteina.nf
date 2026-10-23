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

// Boltz full refolding params
params.full_refold_with = '' // comma-separated list of methods, e.g., 'boltz'
params.full_refold_max = false
params.full_refold_filter_sort = 'pair_pae_min'
params.full_refold_use_msa_server = false
params.full_refold_create_target_msa = false
params.full_refold_target_templates = false
params.full_refold_target_fasta = false
params.full_refold_alignment = false  // external A3M file; bypasses MMseqs2 when set
params.uniref30 = false
params.colabfold_envdb = false

// --- Oracle refold scoring (REFOLD_ORACLES) ---
// Nextflow does not evaluate top-level params statements in an included script, so this
// block lives here rather than in subworkflows/local/refold_oracles.nf. Copy it when
// adopting that subworkflow in another workflow.
params.fold_publish_dir = 'rfd4_proteina/refold_oracles'

params.refold_oracles = ''
params.refold_oracles_max = false
params.refold_oracle_metrics = 'iptm_max,iptm_mean,ipsae_max,ipsae_mean,n'
params.refold_oracle_target_id = 'target'

// --- Shared fold-stack parameters (the fold / fold_pulldown entrypoints declare
// their own copies; everything else these engines read lives in nextflow.config) ---
params.msa_method = 'jackhmmer_af2'
params.create_target_msa = false
params.create_binder_msa = false
params.n_predictions = 1
params.use_msa_server = false
params.templates = false
params.msa_subsample = false
params.msa_subsample_include_full = true

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
include { BOLTZ_REFOLD_SCORING_RFD4 } from '../subworkflows/local/boltz_refold_scoring_rfd4'
// Needs the oracle/fold-stack params declared above.
include { REFOLD_ORACLES } from '../subworkflows/local/refold_oracles'
include { FOLDSEEK_SEARCH } from '../subworkflows/local/foldseek_search'
include { FOLDSEEK_PREPARE_QUERIES } from '../modules/local/foldseek/foldseek_prepare_queries'

workflow RFD4_PROTEINA {

    main:

    if (params.rfd4_config && params.input_pdb) {
        throw new Exception('--rfd4_config and --input_pdb are mutually exclusive. Use --rfd4_config (input paths are read from the spec) or --input_pdb with --contigs.')
    }

    if (params.full_refold_create_target_msa && params.full_refold_alignment) {
        throw new Exception('--full_refold_create_target_msa and --full_refold_alignment are mutually exclusive. Use one or the other.')
    }
    if (params.full_refold_create_target_msa && !params.full_refold_alignment && !params.full_refold_use_msa_server && (!params.colabfold_envdb || !params.uniref30)) {
        throw new Exception('--full_refold_create_target_msa with local MSA (no --full_refold_alignment) requires --colabfold_envdb and --uniref30.')
    }
    if (params.full_refold_use_msa_server && !params.full_refold_create_target_msa && !params.full_refold_alignment) {
        throw new Exception('--full_refold_use_msa_server only has an effect when either --full_refold_create_target_msa is true or --full_refold_alignment is set.')
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

        Full refolding options (Boltz-2):
            --full_refold_with         Comma-separated list of refolding methods (valid options: 'boltz') [default: ${params.full_refold_with}]
            --full_refold_max          Maximum designs to refold [default: ${params.full_refold_max}]
            --boltz_refold_batch_size  Designs per Boltz refold task; >1 amortises Boltz model loading over the batch [default: ${params.boltz_refold_batch_size}]
            --full_refold_create_target_msa  Create target MSA for refolding (runs MMseqs2) [default: ${params.full_refold_create_target_msa}]
            --full_refold_alignment    External A3M file for refold target; bypasses MMseqs2 when set [default: ${params.full_refold_alignment}]
            --full_refold_use_msa_server     Use ColabFold MMseqs2 API for refold target MSA (skip local MMseqs) [default: ${params.full_refold_use_msa_server}]
            --full_refold_target_fasta       Full-length target FASTA for Boltz refold (headers match design basenames) [default: ${params.full_refold_target_fasta}]
            --full_refold_target_templates   Directory of full-length target template PDB/CIF for Boltz [default: ${params.full_refold_target_templates}]
            --colabfold_envdb          ColabFold environment database (local MSA when not using MSA server) [default: ${params.colabfold_envdb}]
            --uniref30                 UniRef30 database (local MSA when not using MSA server) [default: ${params.uniref30}]

        Oracle refold scoring (prototype; independent of --full_refold_with, both may be used):
            Co-folds each design's binder sequence against the target with any set of folding
            oracles, reusing the --method fold_pulldown machinery, and adds tool-prefixed
            iptm / ipSAE columns plus a cross-tool consensus_z to combined_scores.tsv. It does
            not give the RMSD-to-design or binder-monomer metrics --full_refold_with boltz does,
            since it works from the binder sequence rather than the design structure.

            --refold_oracles           Comma-separated oracles: af2, af2_mono, boltz, rf3,
                                       protenix, af3, openfold3, chai1, esmfold2, esmfold2_fast
                                       [default: disabled]
            --refold_oracles_max       Maximum designs to score [default: --full_refold_max]
            --refold_oracle_metrics    Summary columns carried into combined_scores.tsv
                                       [default: ${params.refold_oracle_metrics}]
            --n_predictions            Structures per design per oracle [default: ${params.n_predictions}]
            --msa_method               jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]
            --create_target_msa        Build an MSA for the target [default: ${params.create_target_msa}]
            --use_msa_server           Let the engine fetch its own MSA [default: ${params.use_msa_server}]
            --consensus_metric / --z_stat / --z_scope / --min_pool
                                       Ranking, as in --method fold_pulldown

            Every engine flag from --method fold --help (--af2_*, --boltz_*, --rf3_*,
            --protenix_*, --af3_*, --openfold3_*, --chai1_*, --esmfold2_*) applies here too.

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

    def refold_methods = params.full_refold_with ? params.full_refold_with.toString().split(',').collect { it.trim() } : []
    def refold_oracles = params.refold_oracles
        ? FoldValidation.parseMethods(params.refold_oracles)
        : []
    if (refold_oracles) {
        // Every design pair here is a 2-chain complex and the binder is folded query-only,
        // exactly as in --method fold_pulldown, so pulldown: true applies for the same reasons.
        def (oracle_errors, oracle_warnings) = FoldValidation.validate(params, refold_oracles, [
            hasMultimer: true,
            af2DbPath: params.af2_db_path,
            pulldown: true,
        ])
        oracle_warnings.each { log.warn("rfd4_proteina: --refold_oracles: ${it}") }
        if (oracle_errors) {
            error("rfd4_proteina: --refold_oracles: ${oracle_errors.join('\nrfd4_proteina: --refold_oracles: ')}")
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

    def ch_boltz_complex = Channel.empty()
    def ch_boltz_monomer = Channel.empty()
    def ch_boltz_rmsd_target_aligned_binder = Channel.empty()
    def ch_boltz_rmsd_monomer_vs_complex = Channel.empty()

    if (refold_methods.contains('boltz')) {
        BOLTZ_REFOLD_SCORING_RFD4(
            ch_designs_with_meta,
            rfd4BinderChain,
            rfd4TargetChain,
            params.full_refold_max,
            params.full_refold_create_target_msa,
            params.full_refold_use_msa_server,
            params.full_refold_alignment,
            params.full_refold_target_fasta,
            params.full_refold_target_templates,
            params.colabfold_envdb,
            params.uniref30,
            params.outdir,
        )
        ch_boltz_complex = BOLTZ_REFOLD_SCORING_RFD4.out.boltz_scores_complex
        ch_boltz_monomer = BOLTZ_REFOLD_SCORING_RFD4.out.boltz_scores_monomer
        ch_boltz_rmsd_target_aligned_binder = BOLTZ_REFOLD_SCORING_RFD4.out.rmsd_target_aligned_binder
        ch_boltz_rmsd_monomer_vs_complex = BOLTZ_REFOLD_SCORING_RFD4.out.rmsd_monomer_vs_complex
    }

    // Oracle-agnostic refold scoring: co-fold each design's binder sequence against the
    // target with any set of folding oracles. Independent of --full_refold_with, which
    // drives the Boltz-specific RMSD/ipSAE/BindCraft path above; a run can use both.
    def ch_oracle_scores = Channel.empty()
    if (refold_oracles) {
        REFOLD_ORACLES(
            ch_designs_with_meta,
            rfd4BinderChain,
            rfd4TargetChain,
            params.refold_oracles_max ?: params.full_refold_max,
            refold_oracles,
            params.msa_method,
        )
        ch_oracle_scores = REFOLD_ORACLES.out.scores
    }

    if (params.do_foldseek) {
        def foldseek_max = (params.foldseek_search_max ?: params.full_refold_max)
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
    def ch_design_cifs = ch_designs
        .ifEmpty(Channel.of(file("${projectDir}/assets/dummy_files/empty")))
        .collect()
        .map { cifs -> [cifs as List] }

    def ch_combine_input = ch_design_table
        .combine(ch_boltz_complex.ifEmpty(file("${projectDir}/assets/dummy_files/combine_placeholder_boltz_complex")))
        .combine(ch_boltz_monomer.ifEmpty(file("${projectDir}/assets/dummy_files/combine_placeholder_boltz_monomer")))
        .combine(ch_boltz_rmsd_target_aligned_binder.ifEmpty(file("${projectDir}/assets/dummy_files/combine_placeholder_boltz_rmsd_target_aligned_binder")))
        .combine(ch_boltz_rmsd_monomer_vs_complex.ifEmpty(file("${projectDir}/assets/dummy_files/combine_placeholder_boltz_monomer_vs_complex")))
        .combine(ch_oracle_scores.ifEmpty(file("${projectDir}/assets/dummy_files/combine_placeholder_refold_oracles")))
        .combine(Channel.value(rfd4BinderChain))
        .combine(ch_design_cifs)
        .map { designs, boltz_cx, boltz_mono, boltz_rmsd_tab, boltz_rmsd_mono, oracle_scores, binder_chain, cifs_wrapped ->
            def cifs = (cifs_wrapped instanceof List && cifs_wrapped.size() == 1 && cifs_wrapped[0] instanceof List)
                ? cifs_wrapped[0]
                : cifs_wrapped
            tuple(designs, boltz_cx, boltz_mono, boltz_rmsd_tab, boltz_rmsd_mono, oracle_scores, binder_chain, cifs)
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
