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
params.rfd4_recycle_strategy = 'medium/front_loaded'
params.rfd4_inference_experiment = false
// Which metrics run inline, as a comma-separated tag list. The defaults run RFD4-Proteina's own
// refolding of the complex (ppi_refolding) plus the hotspot conditioning checks, using the
// OpenFold3 weights baked into the container (AlphaFold 3 weights, if wanted, must be mounted at
// runtime). Set false to drop the metrics group entirely so no refolding oracle is built. Ignored
// by --rfd4_inference_experiment quickstart, which does not compose the metrics group.
//
// apo_refolding (binder alone) is available but off by default: it is a second refold per design,
// and it only yields refolded_apo_holo_rmsd when paired with a tag refolding the untransformed
// complex, as ppi_refolding does.
params.rfd4_metrics_tags = 'ppi_refolding,hotspot_conditioning'
// false runs the trunk eagerly (compile_lengths=null), skipping the torch.compile warm-up
params.rfd4_compile_lengths = true
params.rfd4_partial_trajectory_fraction = false
params.rfd4_inverse_folding = false
params.rfd4_extra_args = ''
params.rfd4_filters = false
// A Hydra config directory prepended to rfd4-generate's config search path. The default
// repairs the target-aligned binder RMSD metrics for C_CRDGRP-conditioned targets; see
// assets/rfd4_proteina/hydra_overrides. Point it at your own directory to add or override
// other configs (pass_criteria/, refolding_metrics/, ...), or set false for none.
params.rfd4_metrics_config_dir = "${projectDir}/assets/rfd4_proteina/hydra_overrides"

// Chain ids in the generated complex. Derived from --contigs in params mode; set these when
// using a condition spec whose target chain is not A.
params.rfd4_target_chain = false
params.rfd4_binder_chain = false

// Refold scoring. Its knobs - --refold_with, --refold_filter, --refold_sort_by,
// --refold_max and what it reports - are shared by every design workflow and live in
// conf/refold.config; nextflow.config derives --refold_publish_dir and --fold_publish_dir
// from --method.

// The fold stack REFOLD drives - MSA options, --n_predictions, the per-engine refold_*
// copies, AF2, Protenix, RF3, EnGens and FoldSeek - is shared by every design workflow
// and lives in conf/refold_fold_stack.config.

params.require_gpu = true
params.gpu_devices = ''
params.gpu_slots_per_device = 1
params.gpu_lock_timeout = 14400

include { UNIQUE_ID } from '../modules/local/common/unique_id'
include { PARSE_RFD4_CONFIG } from '../modules/local/rfd4_proteina/parse_rfd4_config'
include { GENERATE_RFD4_CONFIG } from '../modules/local/rfd4_proteina/generate_rfd4_config'
include { PREPARE_RFD4_INPUT } from '../modules/local/rfd4_proteina/prepare_rfd4_input'
include { RFD4_GENERATE; rfd4IsFalse } from '../modules/local/rfd4_proteina/rfd4_generate'
include { FILTER_DESIGNS as RFD4_FILTER_DESIGNS } from '../modules/local/rfd4_proteina/filter_designs'
include { COMBINE_RFD4_SCORES } from '../modules/local/rfd4_proteina/combine_rfd4_scores'
include { resolveRfd4Chains } from '../modules/local/rfd4_proteina/rfd4_utils'
include {
    checkRefoldParamRenames
    refoldEngines
    refoldHelpText
    refoldSelection
} from '../modules/local/common/refold_select_utils'
// Needs the refold / fold-stack params declared above.
include { REFOLD } from '../subworkflows/local/refold'
include { FOLDSEEK_SEARCH } from '../subworkflows/local/foldseek_search'
include { FOLDSEEK_PREPARE_QUERIES } from '../modules/local/foldseek/foldseek_prepare_queries'
include { foldValidationValidate } from '../modules/local/fold/common/fold_validation'

workflow RFD4_PROTEINA {

    main:

    if (params.rfd4_config && params.input_pdb) {
        throw new Exception('--rfd4_config and --input_pdb are mutually exclusive. Use --rfd4_config (input paths are read from the spec) or --input_pdb with --contigs.')
    }

    checkRefoldParamRenames(params, 'rfd4_proteina')

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
                                  with residue numbers is preserved from the input. The designed
                                  segment is named BINDER in the generated spec, which excludes
                                  cysteines from it; name it yourself with "BINDER: 50-120, ..."
                                  to pick a different name. [default: ${params.contigs}]
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
            --rfd4_recycle_strategy  Trunk recycle schedule, eg thinking_medium, high/midpeak [default: ${params.rfd4_recycle_strategy}]
            --rfd4_inference_experiment  Experiment preset, eg quickstart [default: inference_on_user_inputs]
            --rfd4_metrics_tags   Comma-separated inline metric tags; false drops the metrics group.
                                  Refolding uses the OpenFold3 weights baked into the container.
                                  Add apo_refolding to also refold the binder alone; it needs
                                  ppi_refolding alongside it for the apo/holo RMSD.
                                  Ignored with the quickstart preset. [default: ${params.rfd4_metrics_tags}]
            --rfd4_compile_lengths  Set false to run the trunk eagerly, skipping the torch.compile
                                  warm-up. Worth it for a smoke test, not for a real run. [default: ${params.rfd4_compile_lengths}]
            --rfd4_partial_trajectory_fraction  Partial diffusion: run only this final fraction of
                                  nsteps, refining the input rather than building from noise [default: disabled]
            --rfd4_inverse_folding  Hold the input backbone fixed and resample sequence only [default: ${params.rfd4_inverse_folding}]
            --rfd4_extra_args     Additional Hydra key=value overrides for rfd4-generate [default: ${params.rfd4_extra_args}]
            --rfd4_metrics_config_dir  Hydra config directory prepended to rfd4-generate's search
                                  path, for dropping in custom refolding_metrics/, pass_criteria/
                                  and so on. The default fixes the target-aligned binder RMSD
                                  metrics for C_CRDGRP targets; false uses the container's configs
                                  unchanged. [default: assets/rfd4_proteina/hydra_overrides]
            --rfd4_refold_oracle  Names the subdirectory the inline refolds are published to,
                                  refolded/<oracle>/. Set alphafold3 if you mounted AlphaFold 3
                                  weights; the container refolds with OpenFold3 otherwise.
                                  [default: ${params.rfd4_refold_oracle}]
            --rfd4_gzip_batch_cifs  Gzip the .cif files under batches/ [default: ${params.rfd4_gzip_batch_cifs}]
            --rfd4_filters        Semicolon-separated filters for generated designs, eg "rg<25" [default: disabled]
            --rfd4_target_chain   Target chain id in the generated complex [default: derived from --contigs or --rfd4_config, else A]
            --rfd4_binder_chain   Binder chain id in the generated complex [default: derived from --contigs or --rfd4_config, else B]

        ${refoldHelpText(params)}

        FoldSeek structural search (optional, enabled with --do_foldseek):
            --do_foldseek                 Enable FoldSeek search on the same designs selected for full refolding
            --foldseek_search_max         Maximum designs to search with FoldSeek [default: same as --refold_max]
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

    refold_engines = refoldEngines(params, 'rfd4_proteina', log)
    refold_selection = refoldSelection(params, 'rfd4_proteina', log)

    if (refold_engines) {
        // Every design pair here is a 2-chain complex and the binder is folded query-only,
        // exactly as in --method fold_pulldown, so pulldown: true applies for the same reasons.
        refold_validation = foldValidationValidate(params, refold_engines, [
            hasMultimer: true,
            af2DbPath: params.af2_db_path,
            pulldown: true,
        ])
        refold_errors = refold_validation[0]
        refold_warnings = refold_validation[1]
        refold_warnings.each { log.warn("rfd4_proteina: --refold_with: ${it}") }
        if (refold_errors) {
            error("rfd4_proteina: --refold_with: ${refold_errors.join('\nrfd4_proteina: --refold_with: ')}")
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

    chain_pair = resolveRfd4Chains(params)
    rfd4TargetChain = chain_pair[0]
    rfd4BinderChain = chain_pair[1]

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

    def ch_hydra_overrides = Channel.value(
        file(
            rfd4IsFalse(params.rfd4_metrics_config_dir)
                ? "${projectDir}/assets/rfd4_proteina/hydra_overrides_none"
                : params.rfd4_metrics_config_dir,
            checkIfExists: true,
        )
    )

    RFD4_GENERATE(
        ch_config_spec,
        ch_input_structures,
        ch_hydra_overrides,
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

    // What --refold_filter and --refold_sort_by get to select on: the inline RFD4-Proteina
    // metrics every design carries (rfd4_ppi_pass, rfd4_oracle_ppi_overall_iptm and the
    // rest, from --rfd4_metrics_tags), plus the backbone metrics --rfd4_filters computed
    // (rg, compactness) when that flag is set. RFD4-Proteina has no cofold of its own
    // before refolding, so those inline metrics are the closest thing to one.
    //
    // The design table is read here rather than taken from the collected rfd4_designs.tsv
    // because that one is assembled after refolding. It is emitted by the same task as the
    // structures, so selecting on it costs no extra wait.
    def ch_rfd4_metrics = RFD4_GENERATE.out.design_table
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.id.toString(), row as Map) }

    // The structural scores are keyed on the file name rather than a column in the table:
    // the scores TSV is named after the design it describes, which holds whatever
    // filter_designs.py chose to call its identifier column.
    def ch_design_scores = ch_rfd4_metrics
    if (params.rfd4_filters) {
        def ch_structural_scores = RFD4_FILTER_DESIGNS.out.scores.map { tsv ->
            def id = tsv.name
                .replaceFirst(/\.scores\.tsv$/, '')
                .replaceFirst(/\.(cif|cif\.gz|cif\.zst)$/, '')
            def rows = tsv.splitCsv(header: true, sep: '\t')
            tuple(id, (rows ? rows[0] : [:]) as Map)
        }
        ch_design_scores = ch_rfd4_metrics
            .join(ch_structural_scores)
            .map { id, metrics, structural -> tuple(id, metrics + structural) }
    }

    def ch_scored_designs = ch_designs_with_meta
        .map { meta, cif -> tuple(meta.id.toString(), meta, cif) }
        .join(ch_design_scores)
        .map { _id, meta, cif, scores -> tuple(meta, cif, scores) }

    // Refold scoring: fold each design's binder against its target with every engine in
    // --refold_with, compare the result with the design, and score the interface.
    ch_refold_scores = Channel.empty()
    if (refold_engines) {
        REFOLD(
            ch_scored_designs,
            rfd4BinderChain,
            rfd4TargetChain,
            refold_selection,
            refold_engines,
            params.msa_method,
        )
        ch_refold_scores = REFOLD.out.scores
    }

    if (params.do_foldseek) {
        def foldseek_max = (params.foldseek_search_max ?: refold_selection.max)
        def ch_foldseek_input = foldseek_max
            ? ch_designs.take(foldseek_max as int)
            : ch_designs

        FOLDSEEK_PREPARE_QUERIES(ch_foldseek_input, 0)
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

    ch_combine_input = ch_design_table
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
