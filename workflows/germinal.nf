#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Germinal binder design workflow

Usage via main.nf:
  nextflow run main.nf --method germinal --germinal_config configs/pdl1_vhh_protenix.yaml

*/

params.outdir = 'results'
params.germinal_config = false
params.germinal_pdb_dir = false
params.germinal_experiment_name = 'germinal_design'
params.germinal_n_traj = 10
params.germinal_batch_size = 1
params.germinal_max_passing_designs = 10000
params.germinal_max_hallucinated_trajectories = 10000
// Refold scoring. Its knobs - --refold_with, --refold_filter, --refold_sort_by,
// --refold_max and what it reports - are shared by every design workflow and live in
// conf/refold.config; nextflow.config derives --refold_publish_dir and --fold_publish_dir
// from --method.

params.require_gpu = true
params.gpu_devices = ''
params.gpu_slots_per_device = 1
params.gpu_lock_timeout = 14400

include { GERMINAL as GERMINAL_PROCESS } from '../modules/local/germinal/germinal'
include { GERMINAL_MERGE } from '../modules/local/germinal/germinal_merge'
include { REFOLD } from '../subworkflows/local/refold'
include {
    checkRefoldParamRenames
    refoldEngines
    refoldHelpText
    refoldIdColumn
    refoldSelection
} from '../modules/local/common/refold_select_utils'
include { foldValidationValidate } from '../modules/local/fold/common/fold_validation'

def resolveGerminalPdbDir(config_path) {
    if (params.germinal_pdb_dir) {
        return file(params.germinal_pdb_dir, checkIfExists: true)
    }

    def inferred = file("${config_path.parent}/../pdbs").normalize()
    if (!inferred.exists()) {
        error "Could not infer --germinal_pdb_dir from config location. Expected ${inferred} or pass --germinal_pdb_dir explicitly."
    }
    return inferred
}

workflow GERMINAL {

    main:

    if (!params.germinal_config) {
        log.info"""
        ==================================================================
        GERMINAL NEXTFLOW WRAPPER
        ==================================================================

        Required arguments:
            --germinal_config              Germinal Hydra config YAML (combined or partial)

        Optional arguments:
            --outdir                       Output directory [default: ${params.outdir}]
            --germinal_pdb_dir             Directory with target PDB and optional nb.pdb scaffold
                                           [default: ../pdbs relative to config directory]
            --germinal_experiment_name     Experiment name for output subdirectory [default: ${params.germinal_experiment_name}]
            --germinal_n_traj              Total number of trajectory attempts [default: ${params.germinal_n_traj}]
            --germinal_batch_size          Trajectories per batch (parallel task) [default: ${params.germinal_batch_size}]
            --germinal_max_passing_designs Max accepted designs per batch [default: ${params.germinal_max_passing_designs}]
            --germinal_max_hallucinated_trajectories
                                           Max hallucinated trajectories per batch [default: ${params.germinal_max_hallucinated_trajectories}]

            --require_gpu                  Fail if no GPU detected [default: ${params.require_gpu}]
            --gpu_devices                  GPU devices for allocation [default: ${params.gpu_devices}]
            --gpu_slots_per_device         Concurrent tasks allowed per GPU [default: ${params.gpu_slots_per_device}]
            --gpu_lock_timeout             Seconds a task waits for a free GPU [default: ${params.gpu_lock_timeout}]

        ${refoldHelpText(params)}

        """.stripIndent()
        exit 1
    }

    // Refolding is the held-out second opinion after Germinal's own AF2 pass has
    // already cofolded every accepted design.
    checkRefoldParamRenames(params, 'germinal')
    refold_engines = refoldEngines(params, 'germinal', log)
    refold_selection = refoldSelection(params, 'germinal', log)

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
        refold_warnings.each { log.warn("germinal: --refold_with: ${it}") }
        if (refold_errors) {
            error("germinal: --refold_with: ${refold_errors.join('\ngerminal: --refold_with: ')}")
        }
    }

    if (params.germinal_n_traj < 1) {
        error "germinal_n_traj must be >= 1"
    }
    if (params.germinal_batch_size < 1) {
        error "germinal_batch_size must be >= 1"
    }

    ch_config = Channel.fromPath(params.germinal_config, checkIfExists: true).first()
    ch_pdb_dir = ch_config.map { config ->
        resolveGerminalPdbDir(config)
    }

    def design_indices = 0..(params.germinal_n_traj - 1)
    def batches = design_indices.collate(params.germinal_batch_size)

    ch_batch_info = Channel.from(batches.withIndex())
        .map { batch, index -> [index, batch.size()] }

    GERMINAL_PROCESS(
        ch_batch_info,
        ch_config,
        ch_pdb_dir,
        params.germinal_experiment_name,
        params.germinal_max_passing_designs,
        params.germinal_max_hallucinated_trajectories
    )

    GERMINAL_MERGE(
        GERMINAL_PROCESS.out.batch_dir.collect(),
    )

    // Refold scoring (optional). Germinal has already cofolded every accepted design
    // with its own AF2 pass, so this is the held-out second opinion: fold the binder
    // against the target again with the engines --refold_with names, and score the
    // result the way every other workflow does.
    //
    // accepted_designs.csv is Germinal's own verdict on each design, and what
    // --refold_filter and --refold_sort_by select on. GERMINAL_MERGE produces it only
    // once every batch has finished, so refolding cannot overlap designing here.
    ch_refold_scores = Channel.empty()
    if (refold_engines) {
        def ch_germinal_score_rows = GERMINAL_MERGE.out.accepted_designs
            // quote: '"' because these are real CSVs, not comma-joined columns:
            // BindCraft writes InterfaceResidues as "B16,B19,B20,..." and several
            // amino-acid-count fields as quoted Python dicts. Nextflow's splitCsv does
            // no quote handling unless asked, so without this every embedded comma
            // starts a new column and each row's later fields shift left by however
            // many commas that row happened to contain. Nothing errors: the metrics
            // just silently come from the wrong columns, and a sort on Average_i_pTM
            // ranks designs on Average_pTM instead.
            .splitCsv(header: true, quote: '"')
            .map { row ->
                def id_column = refoldIdColumn(
                    row as Map,
                    ['design_name', 'design', 'Design', 'description', 'id', 'name'],
                    'germinal',
                )
                tuple(row[id_column].toString(), row as Map)
            }

        // remainder: true so a design with no row in the merged table still reaches the
        // refold rather than vanishing from it; the filter drops the mirror case, a row
        // whose structure was not among the accepted PDBs.
        def ch_germinal_scored = GERMINAL_PROCESS.out.accepted_pdbs
            .flatten()
            .map { pdb -> tuple(pdb.baseName, pdb) }
            .join(ch_germinal_score_rows, remainder: true)
            .filter { _id, pdb, _scores -> pdb != null }
            .map { id, pdb, scores -> tuple([id: id], pdb, (scores ?: [:]) as Map) }

        REFOLD(
            ch_germinal_scored,
            params.refold_binder_chain ?: 'B',
            params.refold_target_chain ?: 'A',
            refold_selection,
            refold_engines,
            params.msa_method,
        )
        ch_refold_scores = REFOLD.out.scores
    }

    emit:
    all_trajectories = GERMINAL_MERGE.out.all_trajectories
    accepted_designs = GERMINAL_MERGE.out.accepted_designs
    failure_counts = GERMINAL_MERGE.out.failure_counts
    accepted_pdbs = GERMINAL_PROCESS.out.accepted_pdbs
    refold_scores = ch_refold_scores
}
