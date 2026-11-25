#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
RFDiffusion3-based binder design workflow

Usage via main.nf:
  nextflow run main.nf --method rfd3 --rfd3_config config.json
  nextflow run main.nf --method rfd3 --input_pdb target.pdb --contigs "A17-131,/0,50-120" --hotspot_res "A56,A115,A123"

*/

// Default parameters
params.input_pdb = false
params.outdir = 'results'

params.design_name = 'rfd3'
params.contigs = ''
params.hotspot_res = false

// RFDiffusion3 params
params.rfd3_config = false
params.rfd3_n_designs = 1
params.rfd3_batch_size = 1
params.rfd3_hotspot_subsample = 1.0
// BindCraft-compatible alias; when unset (null), --rfd3_hotspot_subsample is used
params.hotspot_subsample = null
params.rfd3_step_scale = 3
params.rfd3_gamma_0 = 0.2
params.rfd3_is_non_loopy = null  // null = omit from config; true/false = add to config (params mode only; with --rfd3_config set in JSON)
params.rfd3_extra_args = ''

params.rfd3_filters = false

// RosettaFold3 params
// Default path is the one that exists in the "...-weights" container
//params.rf3_ckpt_path = '/weights/rf3_foundry_01_24_latest_remapped.ckpt'
params.rf3_ckpt_path = '/models/foundry/rf3_foundry_01_24_latest_remapped.ckpt'

// RosettaFold3 MSA/template params
params.rf3_create_target_msa = false
params.rf3_use_msa_server = false
params.rf3_target_fasta = false
params.rf3_target_msa = false  // external A3M file; bypasses MMseqs2 when set
// Removed: --rf3_alignment is --rf3_target_msa. Still declared so that setting it stops
// the run with the replacement named, rather than drawing Nextflow's generic
// unknown-parameter warning while the run carries on without the MSA it was given.
params.rf3_alignment = false
params.rf3_target_template = false
params.rf3_template_selection = false

// MPNN params - new (mpnn_*)
params.mpnn_model_type = 'protein_mpnn'
params.mpnn_legacy_weights = true
params.mpnn_designed_chains = 'auto'
// All MPNN batch/temperature/noise/omit params below default to `false` as an
// explicit "not set by user" sentinel; resolveParam() (rfd3_utils.nf) applies
// the true default (see buildMpnnArgs fallback args) only when neither the
// modern nor legacy name was set. Do not change these to real default values
// or the legacy pmpnn_* names will be silently ignored again.
params.mpnn_batch_size = false
params.mpnn_temperature = false  // true default: 0.1
params.mpnn_structure_noise = false  // true default: 0
params.mpnn_omit = false  // true default: 'CX'
params.mpnn_checkpoint_path = false
params.mpnn_preset = false
params.mpnn_weights_noise = false

// Legacy MPNN params (pmpnn_*) - for backward compatibility with rfd workflow
params.pmpnn_seqs_per_struct = false  // mpnn_batch_size; true default: 1
params.pmpnn_temperature = false      // mpnn_temperature; true default: 0.1
params.pmpnn_augment_eps = false      // mpnn_structure_noise; true default: 0
params.pmpnn_omit_aas = false         // mpnn_omit (with 1-letter -> 3-letter conversion); true default: 'CX'
params.pmpnn_weights = false          // mpnn_checkpoint_path

// RosettaFold3 params
params.rf3_early_stopping_plddt_threshold = 0.5  // exits early if mean pLDDT < 0.5 after the first recycle
params.rf3_num_steps = 50                        // default for rf3 cli is 200, however 50 is faster with no difference in quality
params.rf3_n_recycles = 10                       // default for rf3 cli is 10
params.rf3_diffusion_batch_size = 5              // default for rf3 cli is 5
params.rf3_batch_size = 1                        // MPNN designs per rf3 fold task (distinct from --rfd3_batch_size for diffusion)

// Refold scoring. Its knobs - --refold_with, --refold_filter, --refold_sort_by,
// --refold_max and what it reports - are shared by every design workflow and live in
// conf/refold.config; nextflow.config derives --refold_publish_dir and --fold_publish_dir
// from --method.

params.require_gpu = true
params.gpu_devices = ''
params.gpu_slots_per_device = 1
params.gpu_lock_timeout = 14400

include { UNIQUE_ID } from '../modules/local/common/unique_id'
include { RFDIFFUSION3 } from '../modules/local/rfd3/rfdiffusion3'
include { GENERATE_RFD3_CONFIG } from '../modules/local/rfd3/generate_rfd3_config'
include { MPNN } from '../modules/local/rfd3/mpnn'
include { PDB_TO_FASTA } from '../modules/local/common/pdb_to_fasta'
include { MMSEQS_COLABFOLDSEARCH } from '../modules/local/common/mmseqs_colabfoldsearch'
include { GENERATE_RF3_INPUT_JSON } from '../modules/local/rfd3/generate_rf3_input'
include { PREPARE_RF3_TEMPLATE; RENAME_RF3_TEMPLATE_CHAINS } from '../modules/local/rfd3/prepare_rf3_template'
include { ROSETTAFOLD3 } from '../modules/local/rfd3/rosettafold3'
include { RFD3_RMSD } from '../modules/local/rfd3/rfd3_rmsd'
include { COMBINE_RFD3_SCORES } from '../modules/local/rfd3/combine_rfd3_scores'
include { FILTER_DESIGNS as RFD3_FILTER_DESIGNS } from '../modules/local/rfd3/filter_designs'
include {
    buildMpnnArgs;
    canonicalizeMpnnWeightsNoiseParam;
    normaliseContigToV3;
    mpnnDesignedChainsFirst;
    extractRfd3InputPaths;
    resolveRfd3TargetBinderChains;
    validateRfd3MpnnPresetParams;
} from '../modules/local/rfd3/rfd3_utils'
include { REFOLD } from '../subworkflows/local/refold'
include { SELECT_FOR_REFOLD as SELECT_FOR_FOLDSEEK } from '../subworkflows/local/refold_select'
include {
    checkRefoldParamRenames
    refoldEngines
    refoldHelpText
    refoldSelection
} from '../modules/local/common/refold_select_utils'
include { FOLDSEEK_SEARCH } from '../subworkflows/local/foldseek_search'
include { FOLDSEEK_PREPARE_QUERIES } from '../modules/local/foldseek/foldseek_prepare_queries'
include { RFD3_STAGE_FOLDSEEK_CIF } from '../modules/local/rfd3/stage_foldseek_cif'
include { foldValidationValidate } from '../modules/local/fold/common/fold_validation'

workflow RFD3 {

    main:

    if (params.rfd3_config && params.input_pdb) {
        throw new Exception('--rfd3_config and --input_pdb are mutually exclusive. Use --rfd3_config (input path is read from the config) or --input_pdb with --contigs.')
    }

    if (params.rfd3_config && params.rfd3_is_non_loopy != null) {
        throw new Exception('--rfd3_is_non_loopy cannot be used with --rfd3_config. Put "is_non_loopy": true or "is_non_loopy": false in your JSON config file instead.')
    }

    if (params.rf3_alignment) {
        throw new Exception('--rf3_alignment has been renamed to --rf3_target_msa, for consistency with every other MSA parameter in the pipeline.')
    }

    if (params.rf3_create_target_msa && params.rf3_target_msa) {
        throw new Exception('--rf3_create_target_msa and --rf3_target_msa are mutually exclusive. Use one or the other.')
    }
    if (params.rf3_create_target_msa && !params.rf3_target_msa && !params.rf3_use_msa_server && (!params.colabfold_envdb || !params.uniref30)) {
        throw new Exception('--rf3_create_target_msa with local MSA (no --rf3_target_msa) requires --colabfold_envdb and --uniref30.')
    }
    if (params.rf3_use_msa_server && (!params.rf3_create_target_msa || params.rf3_target_msa)) {
        throw new Exception('--rf3_use_msa_server only has an effect when --rf3_create_target_msa is true and --rf3_target_msa is not set.')
    }
    if (params.rf3_target_template && params.rf3_target_fasta) {
        throw new Exception('--rf3_target_template and --rf3_target_fasta are mutually exclusive. The target sequence comes from the template file when --rf3_target_template is set.')
    }

    def mpnn_weights_noise = canonicalizeMpnnWeightsNoiseParam(params)
    validateRfd3MpnnPresetParams(params, mpnn_weights_noise)

    def hotspot_sub_frac = (params.hotspot_subsample != null)
        ? (params.hotspot_subsample as Double)
        : (params.rfd3_hotspot_subsample as Double)
    if (hotspot_sub_frac < 0.0 || hotspot_sub_frac > 1.0) {
        throw new Exception('--rfd3_hotspot_subsample / --hotspot_subsample must be between 0.0 and 1.0')
    }
    if (hotspot_sub_frac < 1.0 && !params.rfd3_config && (!params.hotspot_res || params.hotspot_res == false)) {
        log.warn('Hotspot subsampling < 1.0 has no effect in params mode without --hotspot_res')
    }
    if (hotspot_sub_frac < 1.0 && params.rfd3_config) {
        log.warn('Hotspot subsampling < 1.0 only affects design entries in the config that define select_hotspots')
    }

    // false/null means unset (fold workflows use false as sentinel); default to 1
    def rf3_batch_int = (params.rf3_batch_size == null || params.rf3_batch_size == false) ? 1 : (params.rf3_batch_size as int)
    if (rf3_batch_int < 1) {
        throw new Exception('--rf3_batch_size must be >= 1')
    }

    // Refolding takes its MSA and template choices from the fold stack now, so the
    // refold's own copies of them are gone; setting one stops the run with the
    // replacement named.
    checkRefoldParamRenames(params, 'rfd3')

    // Show help message
    if (params.input_pdb == false && params.rfd3_config == false) {
        log.info(
            """
        ==================================================================
        PROTEIN BINDER DESIGN PIPELINE - RFDiffusion3
        ==================================================================

        Required (one of):
            --rfd3_config         Path to RFDiffusion3 JSON/YAML config file
                                  (input path is read from the config, e.g. "designname": {"input": "path/to/target.pdb"})
            --input_pdb           Input PDB/CIF file (params mode only)

        Config mode (--rfd3_config):
            Input file path is taken from the config; do not use --input_pdb.

        Params mode (--input_pdb, generates config automatically):
            --contigs             Contig string for rfd3 (v3 format: "A17-131,/0,50-120")
                                  v1 format is also accepted and auto-translated [default: ${params.contigs}]
            --hotspot_res         Hotspot residues, eg "A56,A115,A123" [default: ${params.hotspot_res}]

        RFDiffusion3 options:
            --outdir              Output directory [default: ${params.outdir}]
            --design_name         Name of the design [default: ${params.design_name}]
            --rfd3_n_designs      Total number of designs [default: ${params.rfd3_n_designs}]
            --rfd3_batch_size     Designs per RFDiffusion3 batch (diffusion_batch_size) [default: ${params.rfd3_batch_size}]
            --rf3_batch_size      MPNN designs per RosettaFold3 rf3 fold task [default: ${params.rf3_batch_size}]
            --rfd3_step_scale     inference_sampler.step_scale [default: ${params.rfd3_step_scale}]
            --rfd3_gamma_0        inference_sampler.gamma_0 [default: ${params.rfd3_gamma_0}]
            --rfd3_hotspot_subsample  Fraction of hotspot residues to keep per RFDiffusion3 batch (ceil(N*fraction), min 1), aliased to --hotspot_subsample [default: ${params.rfd3_hotspot_subsample}]
            --rfd3_is_non_loopy  Set "is_non_loopy" in config: true/false to add, unset to omit (params mode only) [default: omit]
            --rfd3_extra_args     Additional CLI arguments for rfd3 [default: ${params.rfd3_extra_args}]
            --rfd3_filters        Semicolon-separated filters for RFD3 backbones (binder chain follows contig order; use --mpnn_designed_chains if needed), e.g. "rg<25" [default: disabled]

        MPNN options (new names / legacy names):
            --mpnn_model_type / (n/a)                  Model type [default: ${params.mpnn_model_type}]
            --mpnn_legacy_weights / (n/a)              Use legacy weights [default: ${params.mpnn_legacy_weights}]
            --mpnn_designed_chains / (n/a)             Chains for MPNN --designed_chains, or 'auto' (infer binder polymer from contig) [default: ${params.mpnn_designed_chains}]
            --mpnn_batch_size / --pmpnn_seqs_per_struct  Sequences per structure [default: 1]
            --mpnn_temperature / --pmpnn_temperature   Sampling temperature [default: 0.1]
            --mpnn_structure_noise / --pmpnn_augment_eps  Inference-time Gaussian noise on input coordinates (Å); not the training-noise weight tier [default: 0]
            --mpnn_omit / --pmpnn_omit_aas             Omit residue types (1-letter eg "CX") [default: CX]
            --mpnn_checkpoint_path / --pmpnn_weights   Custom weights path; when set, overrides --mpnn_preset / --mpnn_weights_noise [default: ${params.mpnn_checkpoint_path}]
            --mpnn_preset                              Weight family: vanilla, soluble, hyper, or false for default (ProteinMPNN v48_020) [default: ${params.mpnn_preset}]
            --mpnn_weights_noise                       Training-noise tier string: prefer '005','010','020','030' (see docs); use -params-file JSON strings or quoted CLI values to avoid coercion to integers [default: ${params.mpnn_weights_noise}]

        RF3 MSA/template options:
            --rf3_create_target_msa    Create target MSA for RF3 (once per target; runs MMseqs2) [default: ${params.rf3_create_target_msa}]
            --rf3_target_msa           External A3M file for target; bypasses MMseqs2 when set [default: ${params.rf3_target_msa}]
            --rf3_use_msa_server       Use external MSA server for RF3 when creating MSA (skip local MMseqs) [default: ${params.rf3_use_msa_server}]
            --rf3_target_fasta         Target FASTA file for RF3 MSA creation [default: ${params.rf3_target_fasta}]
            --rf3_target_template      Alternative template PDB/CIF for RF3 (not trimmed to contigs; mutually exclusive with --rf3_target_fasta) [default: uses --input_pdb / config input]
            --rf3_template_selection   RF3 template_selection (comma-separated AtomSelection tokens). With --rf3_target_template use chain IDs as in that file; the workflow maps them to rfd3TargetChain for RF3. Omit for whole target chain (auto-mapped). Default without --rf3_target_template: whole target chain (rfd3TargetChain).

        ${refoldHelpText(params)}

        FoldSeek structural search (optional, enabled with --do_foldseek_search):
            --do_foldseek_search          Enable FoldSeek search on the same designs selected for full refolding
            --foldseek_search_max         Maximum designs to search with FoldSeek [default: same as --refold_max]
            --foldseek_database           Database to search [default: CATH50]
                                          (CATH50, PDB, Alphafold/UniProt50, etc.)
            --foldseek_databases_path     Path to local databases directory (auto-downloaded if unset)
            --foldseek_use_webserver      Use FoldSeek web API instead of local search [default: false]
            --foldseek_mode               Search mode [default: 3diaa]
            --foldseek_maxaccept          Max accepted alignments per query [default: 1]
            --foldseek_gzip_output        Gzip TSV output [default: false]
            --foldseek_include_html_output  Include HTML report [default: false]
            --foldseek_cath_names_path    Path to local cath-names.txt (auto-downloaded if unset)
                                          CATH annotation is automatic when database name starts with "CATH"

        Other options:
            --require_gpu           Fail tasks without a GPU [default: ${params.require_gpu}]
            --gpu_devices           GPU devices to use (comma-separated or 'all') [default: ${params.gpu_devices}]
            --gpu_slots_per_device  Concurrent tasks allowed per GPU [default: ${params.gpu_slots_per_device}]
            --gpu_lock_timeout      Seconds a task waits for a free GPU [default: ${params.gpu_lock_timeout}]

        """.stripIndent()
        )
        exit(1)
    }
    // Refolding is the held-out second opinion after RF3 has already cofolded every
    // design, so --refold_with takes any engine the fold stack knows, not just Boltz.
    refold_engines = refoldEngines(params, 'rfd3', log)
    refold_selection = refoldSelection(params, 'rfd3', log)

    if (refold_engines) {
        // Every design pair is a 2-chain complex and the binder is folded query-only,
        // exactly as in --method fold_pulldown.
        refold_validation = foldValidationValidate(params, refold_engines, [
            hasMultimer: true,
            af2DbPath: params.af2_db_path,
            pulldown: true,
        ])
        refold_errors = refold_validation[0]
        refold_warnings = refold_validation[1]
        refold_warnings.each { log.warn("rfd3: --refold_with: ${it}") }
        if (refold_errors) {
            error("rfd3: --refold_with: ${refold_errors.join('\nrfd3: --refold_with: ')}")
        }
    }

    def ch_input_pdb
    def target_pdb_path
    if (params.rfd3_config) {
        def config_file = file(params.rfd3_config)
        def input_paths = extractRfd3InputPaths(config_file.toString(), config_file.parent.toString())
        if (input_paths.isEmpty()) {
            throw new Exception("No 'input' path found in --rfd3_config ${params.rfd3_config}")
        }
        target_pdb_path = file(input_paths[0])
        ch_input_pdb = Channel.fromPath(target_pdb_path).first()
    } else {
        target_pdb_path = file(params.input_pdb)
        ch_input_pdb = Channel.fromPath(params.input_pdb).first()
    }

    // Generate unique ID for this run
    UNIQUE_ID()
    def ch_unique_id = UNIQUE_ID.out.id_file.map { it.text.trim() }.first()
    def ch_design_name = ch_unique_id.map { uid -> "${params.design_name}_${uid}" }

    // Calculate number of batches
    def n_batches = Math.ceil(params.rfd3_n_designs / params.rfd3_batch_size).toInteger()
    if (n_batches < 1) { n_batches = 1 }

    def rf3ChainPair = resolveRfd3TargetBinderChains(params)
    def rfd3TargetChain = rf3ChainPair[0]
    def rfd3BinderChain = rf3ChainPair[1]
    def mpnnRaw = params.mpnn_designed_chains?.toString()?.trim()
    def effectiveMpnnDesignedChains = (mpnnRaw && !mpnnRaw.equalsIgnoreCase('auto')) ? mpnnRaw : rfd3BinderChain
    def mpnnBinderSeqChain = mpnnDesignedChainsFirst(effectiveMpnnDesignedChains)
    def mpnn_args = buildMpnnArgs(params, effectiveMpnnDesignedChains, mpnn_weights_noise)
    def rf3TemplateSelection = params.rf3_template_selection ? params.rf3_template_selection.toString() : rfd3TargetChain
    def rf3UserTemplateSelectionRaw = params.rf3_template_selection ? params.rf3_template_selection.toString() : ''
    def rf3UserTselB64 = java.util.Base64.getEncoder().encodeToString(rf3UserTemplateSelectionRaw.getBytes('UTF-8'))

    def ch_rfd3_startnum
    if (params.rfd3_config) {
        // Config mode: user provides their own JSON/YAML config
        ch_config = Channel.fromPath(params.rfd3_config).first()

        ch_rfd3_startnum = Channel.of(0..n_batches - 1)

        RFDIFFUSION3(
            ch_config,
            ch_input_pdb,
            ch_design_name,
            params.rfd3_batch_size,
            1,  // n_batches=1 per task, we parallelise via Nextflow
            ch_rfd3_startnum,
            params.rfd3_step_scale,
            params.rfd3_gamma_0,
            params.rfd3_extra_args,
            hotspot_sub_frac,
        )
    }
    else {
        // Params mode: generate config from --contigs, --hotspot_res, etc.
        if (!params.contigs) {
            throw new Exception('Either --rfd3_config or --contigs must be provided')
        }

        def is_non_loopy_mode = (params.rfd3_is_non_loopy == null)
            ? 'omit'
            : (params.rfd3_is_non_loopy ? 'non_loopy' : 'allow_loopy')

        GENERATE_RFD3_CONFIG(
            params.design_name,
            ch_input_pdb,
            normaliseContigToV3(params.contigs),
            params.hotspot_res ?: false,
            false,  // partial_t - not used in de novo design
            is_non_loopy_mode,
        )

        ch_rfd3_startnum = Channel.of(0..n_batches - 1)

        RFDIFFUSION3(
            GENERATE_RFD3_CONFIG.out.config_json,
            ch_input_pdb,
            ch_design_name,
            params.rfd3_batch_size,
            1,  // n_batches=1 per task, we parallelise via Nextflow
            ch_rfd3_startnum,
            params.rfd3_step_scale,
            params.rfd3_gamma_0,
            params.rfd3_extra_args,
            hotspot_sub_frac,
        )
    }

    def ch_rfd3_cifs = RFDIFFUSION3.out.cifs.flatten()

    def ch_backbones
    if (params.rfd3_filters) {
        RFD3_FILTER_DESIGNS(
            ch_rfd3_cifs,
            params.rfd3_filters,
            rfd3BinderChain,
            'filter',
        )
        ch_backbones = RFD3_FILTER_DESIGNS.out.accepted
    }
    else {
        ch_backbones = ch_rfd3_cifs
    }

    MPNN(
        ch_backbones,
        mpnn_args,
    )

    def ch_mpnn_with_meta = MPNN.out.cifs.flatten().map { c -> tuple([id: c.baseName], c) }

    // Prepare RF3 template before MSA so --rf3_target_template MSA uses the same structure as RF3 (single chain, target chain ID).
    ch_rf3_target_chain_val = Channel.value(rfd3TargetChain)
    def ch_rf3_template_val
    def ch_rf3_template_selection_for_json
    if (params.rf3_target_template) {
        RENAME_RF3_TEMPLATE_CHAINS(
            Channel.fromPath(params.rf3_target_template).first()
                .combine(ch_rf3_target_chain_val)
                .combine(Channel.value(rf3UserTselB64))
                .map { s, tc, b64 -> tuple(s, tc, b64) }
        )
        ch_rf3_template_val = RENAME_RF3_TEMPLATE_CHAINS.out.structure.first()
        ch_rf3_template_selection_for_json = RENAME_RF3_TEMPLATE_CHAINS.out.selection_mapped
            .map { p -> p.text.trim() }
            .first()
    } else {
        def ch_prepare_input
        if (params.rfd3_config) {
            ch_prepare_input = Channel.of(file(target_pdb_path)).first()
                .combine(Channel.fromPath(params.rfd3_config).first())
                .combine(ch_rf3_target_chain_val)
                .map { s, cfg, tc -> tuple(s, cfg, '', tc) }
        } else {
            ch_prepare_input = Channel.of(file(target_pdb_path)).first()
                .combine(Channel.of(file("${projectDir}/assets/dummy_files/empty_templates")))
                .combine(Channel.of(normaliseContigToV3(params.contigs)))
                .combine(ch_rf3_target_chain_val)
                .map { s, empty, contig, tc -> tuple(s, empty, contig, tc) }
        }
        PREPARE_RF3_TEMPLATE(ch_prepare_input)
        ch_rf3_template_val = PREPARE_RF3_TEMPLATE.out.pdb.first()
        ch_rf3_template_selection_for_json = Channel.value(rf3TemplateSelection)
    }

    // Target MSA: once per target (single target supported). External a3m file, or build once early, or omit.
    def ch_single_rf3_msa
    if (params.rf3_target_msa) {
        ch_single_rf3_msa = Channel.fromPath(params.rf3_target_msa).first()
    } else if (params.rf3_create_target_msa) {
        def ch_target_fasta_with_meta
        if (params.rf3_target_fasta) {
            def target_fasta_file = file(params.rf3_target_fasta)
            ch_target_fasta_with_meta = Channel.of(tuple([id: 'target'], target_fasta_file))
        } else {
            // Derive the target FASTA from the prepared RF3 template (ch_rf3_template_val),
            // which is trimmed to the contigs and has its target chain renamed to
            // rfd3TargetChain by PREPARE_RF3_TEMPLATE / RENAME_RF3_TEMPLATE_CHAINS. This keeps the
            // MSA query in sync with the structure RF3 actually folds, regardless of the input
            // PDB's original target chain ID (which is whatever the contig specifies, e.g. "C").
            PDB_TO_FASTA(ch_rf3_template_val, rfd3TargetChain)
            ch_target_fasta_with_meta = PDB_TO_FASTA.out.map { f -> tuple([id: 'target'], f) }
        }
        rf3_msa_db = params.rf3_use_msa_server ? file("${projectDir}/assets/dummy_files/empty") : params.colabfold_envdb
        rf3_msa_uniref = params.rf3_use_msa_server ? file("${projectDir}/assets/dummy_files/empty") : params.uniref30
        MMSEQS_COLABFOLDSEARCH(ch_target_fasta_with_meta, params.rf3_use_msa_server, rf3_msa_db, rf3_msa_uniref, 'rfd3/mmseqs2')
        ch_single_rf3_msa = MMSEQS_COLABFOLDSEARCH.out.a3m.map { m, a3m ->
            def files = (a3m instanceof List) ? a3m : [a3m]
            def primary = files.find { it.toString().contains('result') } ?: files[0]
            primary
        }
    } else {
        ch_single_rf3_msa = Channel.of(file("${projectDir}/assets/dummy_files/empty_target_msa"))
    }

    // Build RF3 input tuples: (meta, structure_cif, target_msa, template_structure, target_chain, binder_chain, binder_seq_chain, template_selection)
    def ch_rf3_msa_val = ch_single_rf3_msa.first()
    ch_rf3_input = ch_mpnn_with_meta
        .combine(ch_rf3_msa_val)
        .combine(ch_rf3_template_val)
        .combine(Channel.value(rfd3TargetChain))
        .combine(Channel.value(rfd3BinderChain))
        .combine(Channel.value(mpnnBinderSeqChain))
        .combine(ch_rf3_template_selection_for_json)
        .map { meta, cif, msa, template, tc, bc, bseq, tsel -> tuple(meta, cif, msa, template, tc, bc, bseq, tsel) }

    ch_rf3_batched = ch_rf3_input.collate(rf3_batch_int, true).map { batch ->
        def metas = batch.collect { it[0] }
        def cifs = batch.collect { it[1] }
        def row0 = batch[0]
        tuple(metas, cifs, row0[2], row0[3], row0[4], row0[5], row0[6], row0[7])
    }

    GENERATE_RF3_INPUT_JSON(ch_rf3_batched)
    ROSETTAFOLD3(
        GENERATE_RF3_INPUT_JSON.out.with_json,
        Channel.value(rfd3TargetChain),
        Channel.value(rfd3BinderChain),
    )

    def ch_rf3_per_design = ROSETTAFOLD3.out.per_design_bundle.flatMap { metas, pddir ->
        def ml = metas instanceof List ? metas : [metas]
        ml.collect { m -> tuple(m, file("${pddir}/${m.id}/model.cif"), file("${pddir}/${m.id}/scores.tsv")) }
    }
    def ch_rf3_refolded_cif = ch_rf3_per_design.map { m, mod, sc -> tuple(m, mod) }
    def ch_rf3_scores_with_meta = ch_rf3_per_design.map { m, mod, sc -> tuple(m, sc) }

    ch_rmsd_input = ch_mpnn_with_meta
        .join(ch_rf3_refolded_cif)
        .combine(Channel.value(rfd3TargetChain))
        .combine(Channel.value(rfd3BinderChain))
    RFD3_RMSD(ch_rmsd_input)

    ch_refold_scores = Channel.empty()

    // RF3's per-design score row, which is what --refold_filter and --refold_sort_by
    // select on: pair_pae_min, plddt and the rest of the first cofold's verdict.
    def ch_rf3_scored = ch_rf3_refolded_cif
        .join(ch_rf3_scores_with_meta)
        .map { meta, cif, score_tsv ->
            def rows = score_tsv.splitCsv(header: true, sep: '\t')
            tuple(meta, cif, (rows ? rows[0] : [:]) as Map)
        }

    if (refold_engines) {
        REFOLD(
            ch_rf3_scored,
            rfd3BinderChain,
            rfd3TargetChain,
            refold_selection,
            refold_engines,
            params.msa_method,
        )
        ch_refold_scores = REFOLD.out.scores
    }

    // FoldSeek structural search (optional), on the same basis as the refold: the same
    // filter and ranking, with its own cap when --foldseek_search_max is set.
    if (params.do_foldseek_search) {
        SELECT_FOR_FOLDSEEK(
            ch_rf3_scored,
            refold_selection.filter,
            refold_selection.sort_by,
            params.foldseek_search_max ?: refold_selection.max,
            'rfd3/foldseek',
        )

        // RF3 outputs model.cif per design; rename in task work dir (not launchDir)
        RFD3_STAGE_FOLDSEEK_CIF(SELECT_FOR_FOLDSEEK.out.selected)

        // Extract design (shorter/binder) chain from complexes
        FOLDSEEK_PREPARE_QUERIES(RFD3_STAGE_FOLDSEEK_CIF.out, 0)

        def ch_foldseek_pdbs = FOLDSEEK_PREPARE_QUERIES.out.design_chains
        ch_foldseek_meta = Channel.of([id: 'rfd3_foldseek'])

        FOLDSEEK_SEARCH(ch_foldseek_pdbs, ch_foldseek_meta)
    }

    def ch_rmsd_target_aligned_binder = RFD3_RMSD.out.rmsd_target_aligned_binder
        .map { meta, tsv -> tsv }
        .collectFile(name: 'rmsd_target_aligned_binder.tsv', storeDir: "${params.outdir}/rfd3/rosettafold3/rmsd", keepHeader: true, skip: 1)
    def ch_rmsd_complex = RFD3_RMSD.out.rmsd_complex
        .map { meta, tsv -> tsv }
        .collectFile(name: 'rmsd_complex.tsv', storeDir: "${params.outdir}/rfd3/rosettafold3/rmsd", keepHeader: true, skip: 1)
    def ch_rmsd_binder_aligned_binder = RFD3_RMSD.out.rmsd_binder_aligned_binder
        .map { meta, tsv -> tsv }
        .collectFile(name: 'rmsd_binder_aligned_binder.tsv', storeDir: "${params.outdir}/rfd3/rosettafold3/rmsd", keepHeader: true, skip: 1)
    def ch_rmsd_target_aligned_target = RFD3_RMSD.out.rmsd_target_aligned_target
        .map { meta, tsv -> tsv }
        .collectFile(name: 'rmsd_target_aligned_target.tsv', storeDir: "${params.outdir}/rfd3/rosettafold3/rmsd", keepHeader: true, skip: 1)

    ch_rmsd_tuple = ch_rmsd_target_aligned_binder
        .combine(ch_rmsd_complex)
        .combine(ch_rmsd_binder_aligned_binder)
        .combine(ch_rmsd_target_aligned_target)
        .ifEmpty(Channel.of(tuple(
            file("${projectDir}/assets/dummy_files/combine_placeholder_rmsd_target_aligned_binder"),
            file("${projectDir}/assets/dummy_files/combine_placeholder_rmsd_complex"),
            file("${projectDir}/assets/dummy_files/combine_placeholder_rmsd_binder_aligned_binder"),
            file("${projectDir}/assets/dummy_files/combine_placeholder_rmsd_target_aligned_target"),
        )))

    def ch_rf3_scores_merged = ROSETTAFOLD3.out.scores
        .collectFile(name: 'rf3_scores.tsv', keepHeader: true, skip: 1, sort: false)
    def ch_rfd3_scores_merged = RFDIFFUSION3.out.scores
        .collectFile(name: 'rfd3_scores.tsv', keepHeader: true, skip: 1, sort: false)

    // collect() yields a List; combine() flattens Lists into the parent tuple. Wrap as [list]
    // so one combine element stays a List for path(mpnn_cifs, stageAs: 'cifs/*').
    ch_mpnn_cifs = MPNN.out.cifs.flatten()
        .ifEmpty(Channel.of(file("${projectDir}/assets/dummy_files/empty")))
        .collect()
        .map { cifs -> [cifs as List] }

    ch_combine_core = ch_rf3_scores_merged
        .combine(ch_rfd3_scores_merged)
        .combine(ch_rmsd_tuple)
        .combine(ch_refold_scores.ifEmpty(file("${projectDir}/assets/dummy_files/combine_placeholder_refold")))
        .map { rf3, rfd3, rmsd_tab, rmsd_cx, rmsd_bab, rmsd_tat, refold ->
            [
                rf3_scores: rf3,
                rfd3_scores: rfd3,
                rmsd_target_aligned_binder: rmsd_tab,
                rmsd_complex: rmsd_cx,
                rmsd_binder_aligned_binder: rmsd_bab,
                rmsd_target_aligned_target: rmsd_tat,
                refold_scores: refold,
            ]
        }

    // Binder chain before MPNN list: matches COMBINE_RFD3_SCORES input order (val then path list).
    ch_combine_input = ch_combine_core
        .combine(Channel.value(rfd3BinderChain))
        .combine(ch_mpnn_cifs)
        .map { scores, binder_chain, mpnn_cifs_wrapped ->
            // collect()+wrap may arrive as [[cif, ...]] or flat [cif, ...]; [0] alone kept only one CIF.
            def mpnn_cifs = (mpnn_cifs_wrapped instanceof List && mpnn_cifs_wrapped.size() == 1 && mpnn_cifs_wrapped[0] instanceof List)
                ? mpnn_cifs_wrapped[0]
                : mpnn_cifs_wrapped
            tuple(
                scores.rf3_scores,
                scores.rfd3_scores,
                scores.rmsd_target_aligned_binder,
                scores.rmsd_complex,
                scores.rmsd_binder_aligned_binder,
                scores.rmsd_target_aligned_target,
                scores.refold_scores,
                binder_chain,
                mpnn_cifs,
            )
        }
    COMBINE_RFD3_SCORES(ch_combine_input)

    emit:
    rfd3_cifs = RFDIFFUSION3.out.cifs
    mpnn_cifs = MPNN.out.cifs
    mpnn_fastas = MPNN.out.fastas
    rf3_results = ROSETTAFOLD3.out.results
    combined_scores = COMBINE_RFD3_SCORES.out.combined_scores
    binders_fasta = COMBINE_RFD3_SCORES.out.binders_fasta
    rfd3_rmsd_target_aligned_binder = ch_rmsd_target_aligned_binder
    rfd3_rmsd_complex = ch_rmsd_complex
    rfd3_rmsd_binder_aligned_binder = ch_rmsd_binder_aligned_binder
    rfd3_rmsd_target_aligned_target = ch_rmsd_target_aligned_target
    foldseek_tsv = params.do_foldseek_search ? FOLDSEEK_SEARCH.out.tsv : Channel.empty()
    foldseek_tsv_annotated = params.do_foldseek_search ? FOLDSEEK_SEARCH.out.tsv_annotated : Channel.empty()
    foldseek_html = params.do_foldseek_search ? FOLDSEEK_SEARCH.out.html : Channel.empty()
}
