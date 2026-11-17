/*
REFOLD: score binder designs by folding them again, with any engine.

Refolding is not a Boltz feature. A design is scored by folding its binder sequence
against its target with a chosen predictor and asking the same questions whichever
predictor that is: how confident is the interface, did the binder come back where it
was designed to sit, does it hold its fold on its own, and what does the interface
look like once rebuilt. This subworkflow is the single path for all of that, and
`--refold_with` is the single parameter that chooses the engines.

It replaces the pair of paths that came before: BOLTZ_REFOLD_CORE, which folded and
compared inside one Boltz-specific process, and REFOLD_ORACLES, which folded with any
engine but only reported confidence. Nothing here is specific to one design method,
so every design workflow can use it: it takes designs as tuple(meta, structure, scores)
and gives back one row per design.

Refolding is the held-out second opinion, not a workflow's only cofold. Each design
has already been folded once by the pipeline's own predictor - AF2 initial guess in
rfd and rfd_partial, RF3 in rfd3, BindCraft's and Germinal's own AF2 passes - and the
scores from that pass ride along with each design so that SELECT_FOR_REFOLD can choose
which of them are worth the second, more expensive look.

The pieces:

  SELECT_FOR_REFOLD       which designs are worth refolding, from the first cofold
  REFOLD_PREPARE_INPUTS   target and binder sequences, read out of the designs
  FOLD_PULLDOWN_CORE      fold every (target, binder) pair with every engine
  FOLD (monomer)          fold each binder on its own with the same engines
  REFOLD_COMPARE_COMPLEX  RMSD of the refolded complex against the design
  REFOLD_COMPARE_MONOMER  RMSD of the lone binder against design and complex
  BINDCRAFT_SCORING       interface metrics on the refolded complex
  REFOLD_SCORES_WIDE      all of it pivoted to one row per design
  REFOLD_SAMPLE_SCORES    the same, unreduced: one row per sample

Parameters are NOT declared here: Nextflow does not evaluate top-level `params.x = ...`
statements in an included script. They live in conf/refold.config (the --refold_* knobs,
for every entrypoint) and conf/refold_fold_stack.config (the fold stack this drives, for
the design workflows only). A workflow adopting this subworkflow declares just two of
its own, --refold_publish_dir and --fold_publish_dir, which name a different directory
in each, and calls checkRefoldParamRenames / refoldEngines / refoldSelection from
modules/local/common/refold_select_utils.nf to read the rest.
*/

include { FOLD_PULLDOWN_CORE } from './fold_pulldown_core'
include { FOLD as FOLD_BINDER_MONOMER } from './fold'
include { SELECT_FOR_REFOLD } from './refold_select'
include { REFOLD_PREPARE_INPUTS } from '../../modules/local/refold/refold_prepare_inputs'
include { REFOLD_COMPARE_COMPLEX } from '../../modules/local/refold/refold_compare_complex'
include { REFOLD_COMPARE_MONOMER } from '../../modules/local/refold/refold_compare_monomer'
include { REFOLD_STRUCTURE_TO_PDB } from '../../modules/local/refold/refold_structure_to_pdb'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_REFOLD } from '../../modules/local/rfd/bindcraft_scoring'
include { REFOLD_SCORES_WIDE } from '../../modules/local/refold/refold_scores_wide'
include { REFOLD_SAMPLE_SCORES } from '../../modules/local/refold/refold_sample_scores'

workflow REFOLD {

    take:
    ch_designs      // tuple(meta, structure, Map scores): designs and their first-cofold scores
    binder_chain    // String: binder chain id in the design
    target_chain    // String: target chain id in the design
    selection       // Map: [filter: String, sort_by: String, max: Integer] - see SELECT_FOR_REFOLD
    engines         // List<String>: FOLD_PREDICT method names
    msa_method      // 'jackhmmer_af2' | 'mmseqs2_colabfold'

    main:

    SELECT_FOR_REFOLD(
        ch_designs,
        selection.filter,
        selection.sort_by,
        selection.max,
        params.refold_publish_dir ?: 'refold',
    )
    ch_selected = SELECT_FOR_REFOLD.out.selected

    // One task for the whole set: the target is read once and every binder sequence
    // comes out of the same pass, which also makes the "do the designs agree on the
    // target" check possible at all.
    // Ids and structures as two lists built from one collect, so they stay in step:
    // the process links each structure under its design id, because a workflow's own
    // filenames are not always distinct (RF3 writes every design to model.cif).
    ch_refold_inputs = ch_selected
        .collect(flat: false)
        .map { entries ->
            tuple(
                entries.collect { entry -> entry[1] },
                entries.collect { entry -> entry[0].id.toString() },
            )
        }

    // The target to fold against. Unset means the one the designs were built against,
    // read back out of their target chain. --refold_target_fasta replaces it, which is
    // how a design made against a truncated construct is refolded against the
    // full-length protein; see REFOLD_PREPARE_INPUTS for the subsequence check.
    def ch_target_fasta = params.refold_target_fasta
        ? Channel.value(file(params.refold_target_fasta, checkIfExists: true))
        : Channel.value(file("${projectDir}/assets/dummy_files/no_refold_target.fasta"))

    REFOLD_PREPARE_INPUTS(
        ch_refold_inputs.map { structures, _ids -> structures },
        ch_refold_inputs.map { _structures, ids -> ids },
        binder_chain,
        target_chain,
        params.refold_target_id,
        ch_target_fasta,
    )

    // Re-read the FASTAs as (meta, single-record file) pairs, the shape
    // FOLD_PULLDOWN_CORE takes. This mirrors how --method fold_pulldown reads
    // --targets and --binders.
    ch_targets = REFOLD_PREPARE_INPUTS.out.targets
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }
        .merge(REFOLD_PREPARE_INPUTS.out.targets.splitFasta(file: true))

    ch_binders = REFOLD_PREPARE_INPUTS.out.binders
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }
        .merge(REFOLD_PREPARE_INPUTS.out.binders.splitFasta(file: true))

    FOLD_PULLDOWN_CORE(
        ch_targets,
        ch_binders,
        engines,
        msa_method,
        // Query FASTAs for --templates matching; binders only with --binder_templates,
        // as in --method fold_pulldown.
        params.binder_templates
            ? REFOLD_PREPARE_INPUTS.out.targets.mix(REFOLD_PREPARE_INPUTS.out.binders)
            : REFOLD_PREPARE_INPUTS.out.targets,
    )

    // The predicted model structures, out of each engine's own output layout. AF2 is
    // the one engine whose choice of structure depends on run options, so they come
    // along (see FoldNaming.isModelStructure).
    def af2_opts = [
        no_relax: params.af2_no_relax instanceof Boolean
            ? params.af2_no_relax
            : params.af2_no_relax.toString().toBoolean(),
    ]

    // One entry per predicted sample, not per fold job. A job may hold several samples
    // - --n_predictions > 1 with a larger batch size, or an engine handed the whole
    // seed list at once - and each has to be compared and interface-scored on its own.
    // (fold_batch, fold_sample) is the key throughout: the batch separates jobs, and
    // the engine's own sample label separates samples within one.
    ch_complex_models = FOLD_PULLDOWN_CORE.out.predictions
        .flatMap { meta, files ->
            def opts = af2_opts + [keep_best: meta.af2_keep_models == 'best']
            FoldNaming.modelStructures(meta.fold_tool, files, opts).collect { model ->
                tuple(
                    [
                        id: meta.binder,
                        fold_tool: meta.fold_tool,
                        pair_id: meta.id,
                        fold_batch: meta.fold_batch ?: 1,
                        fold_sample: FoldNaming.sampleLabel(meta.fold_tool, model),
                    ],
                    [model],
                )
            }
        }

    // Pair each engine's models back up with the design they came from. combine(by: 0)
    // rather than join(): one design feeds as many rows as there are engines, and join
    // would consume the design after the first.
    ch_design_by_id = ch_selected.map { meta, structure -> tuple(meta.id.toString(), structure) }

    ch_compare_in = ch_complex_models
        .map { meta, models -> tuple(meta.id.toString(), meta, models) }
        .combine(ch_design_by_id, by: 0)
        .map { _id, meta, models, design -> tuple(meta, design, models) }

    REFOLD_COMPARE_COMPLEX(
        ch_compare_in,
        binder_chain,
        target_chain,
        // FOLD_PULLDOWN_MSA always builds its pairs target first, binder second, so a
        // refolded complex is chain A = target, chain B = binder whatever the design's
        // own chain ids are.
        'B',
        'A',
    )

    ch_rmsd_tables = REFOLD_COMPARE_COMPLEX.out.rmsd
        .map { _meta, tsv -> tsv }
        .collectFile(
            name: 'refold_rmsd.tsv',
            storeDir: "${params.outdir}/${params.refold_publish_dir ?: 'refold'}",
            keepHeader: true,
            skip: 1,
        )

    // --- interface scores on the refolded complex ---
    ch_interface_scores = Channel.empty()
    if (params.refold_interface_scoring) {
        // BindCraft reads PDB only, and most engines write mmCIF. Every sample is
        // converted and scored, so the interface metrics exist for each of them in
        // refold_sample_scores.tsv; the wide table then reports the selected one.
        REFOLD_STRUCTURE_TO_PDB(ch_complex_models.map { meta, models -> tuple(meta, models.first()) })

        BINDCRAFT_SCORING_REFOLD(
            REFOLD_STRUCTURE_TO_PDB.out.pdb.map { _meta, pdb -> pdb },
            // The binder in the refolded complex, not in the design: passing the design's
            // chain id here is what made the old Boltz path score the target as if it
            // were the binder whenever the two differed.
            'B',
            'default_4stage_multimer',
            "${params.refold_publish_dir ?: 'refold'}/interface_scores",
        )

        ch_interface_scores = BINDCRAFT_SCORING_REFOLD.out.scores.collectFile(
            name: 'refold_interface_scores.tsv',
            storeDir: "${params.outdir}/${params.refold_publish_dir ?: 'refold'}",
            keepHeader: true,
            skip: 1,
        )
    }

    // --- the binder on its own ---
    ch_monomer_rmsd = Channel.empty()
    if (params.refold_binder_monomer) {
        FOLD_BINDER_MONOMER(
            ch_binders.map { meta, fasta -> tuple(meta + [n_chains: 1], fasta) },
            engines,
            msa_method,
            'monomer',
            // No MSA for the lone binder, matching how the complex refold treats the
            // binder chain and how the Boltz monomer refold did it before: a de novo
            // binder has no homologs to find.
            params.refold_monomer_single_sequence,
        )

        ch_monomer_models = FOLD_BINDER_MONOMER.out.predictions
            .flatMap { meta, files ->
                def opts = af2_opts + [keep_best: meta.af2_keep_models == 'best']
                FoldNaming.modelStructures(meta.fold_tool, files, opts).collect { model ->
                    tuple(
                        [
                            id: meta.id,
                            fold_tool: meta.fold_tool,
                            fold_batch: meta.fold_batch ?: 1,
                            fold_sample: FoldNaming.sampleLabel(meta.fold_tool, model),
                        ],
                        [model],
                    )
                }
            }

        // Key on (design, engine, batch, sample): a monomer is only comparable with the
        // complex the same engine built, and the batch and sample keep the comparison
        // one-to-one. The monomer and the complex folds run the same batching and the
        // same seeds, so a sample label matches its counterpart. Without this every
        // monomer sample would pair with every complex sample, giving a grid of N x N
        // comparisons rather than N.
        def sample_key = { meta ->
            [
                meta.id.toString(),
                meta.fold_tool.toString(),
                meta.fold_batch.toString(),
                meta.fold_sample.toString(),
            ]
        }
        ch_monomer_in = ch_monomer_models
            .map { meta, models -> tuple(sample_key(meta), meta, models) }
            .combine(
                ch_complex_models.map { meta, models -> tuple(sample_key(meta), models.first()) },
                by: 0,
            )
            .map { key, meta, monomers, complex_model -> tuple(key[0], meta, monomers, complex_model) }
            .combine(ch_design_by_id, by: 0)
            .map { _id, meta, monomers, complex_model, design ->
                tuple(meta, design, complex_model, monomers)
            }

        REFOLD_COMPARE_MONOMER(
            ch_monomer_in,
            binder_chain,
            'B',
            params.refold_monomer_chain,
        )

        ch_monomer_rmsd = REFOLD_COMPARE_MONOMER.out.rmsd
            .map { _meta, tsv -> tsv }
            .collectFile(
                name: 'refold_monomer_rmsd.tsv',
                storeDir: "${params.outdir}/${params.refold_publish_dir ?: 'refold'}",
                keepHeader: true,
                skip: 1,
            )
    }

    def empty_table = file("${projectDir}/assets/dummy_files/empty")

    // Collected once and fed to both tables: the wide one reduces each design's
    // samples to a single row, the sample one keeps them all.
    ch_long_tables = ch_rmsd_tables.mix(ch_monomer_rmsd).collect().ifEmpty([empty_table])
    ch_interface_collected = ch_interface_scores.collect().ifEmpty([empty_table])

    REFOLD_SAMPLE_SCORES(
        FOLD_PULLDOWN_CORE.out.scores,
        ch_long_tables,
        ch_interface_collected,
    )

    // The sample table comes first: the wide table reads it to decide which sample
    // each design's structural columns describe.
    REFOLD_SCORES_WIDE(
        FOLD_PULLDOWN_CORE.out.summary,
        ch_long_tables,
        ch_interface_collected,
        REFOLD_SAMPLE_SCORES.out.scores,
        params.refold_metrics,
    )

    emit:
    scores = REFOLD_SCORES_WIDE.out.scores
    sample_scores = REFOLD_SAMPLE_SCORES.out.scores
    fold_scores = FOLD_PULLDOWN_CORE.out.scores
    summary = FOLD_PULLDOWN_CORE.out.summary
    predictions = FOLD_PULLDOWN_CORE.out.predictions
    rmsd = ch_rmsd_tables
    selection = SELECT_FOR_REFOLD.out.report
}
