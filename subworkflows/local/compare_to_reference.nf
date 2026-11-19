/*
COMPARE_TO_REFERENCE: compare predicted complexes against the reference they came from.

Shared by REFOLD (the reference is the design the pipeline just built) and
fold_pulldown's --complex_structures path (the reference is a complex supplied on the
command line). Both ask the same question of the same kind of input, so the model
expansion and the prediction-to-reference join live here rather than in each caller.

Two things make one join serve both:

  - Predictions are keyed by binder id, and `combine` rather than `join` means one
    reference feeds as many rows as there are predictions for its binder. In REFOLD
    that is one row per engine; with several --targets it is one per (engine, target),
    with no change here.
  - FOLD_PULLDOWN_MSA always builds its pairs target first, binder second, so a
    predicted complex is chain A = target, chain B = binder for every engine. Only the
    reference-side chain ids vary, and those ride in each reference's own entry.
*/

include { REFOLD_COMPARE_COMPLEX } from '../../modules/local/refold/refold_compare_complex'
include { foldNamingModelStructures; foldNamingSampleLabel } from '../../modules/local/fold/common/fold_naming'

workflow COMPARE_TO_REFERENCE {

    take:
    ch_predictions  // FOLD_PULLDOWN_CORE.out.predictions: tuple(meta, files)
    ch_references   // tuple(String binder_id, path structure, String binder_chain, String target_chain)
    publish_dir     // String: directory under params.outdir for the collected table
    table_name      // String: filename for the collected table
    label_target    // Boolean: add a 'target' column. Only meaningful when one binder
                    // can be folded against several targets, as in fold_pulldown with
                    // --complex_structures; a refold has one target, where the column
                    // would hold the same value on every row.

    main:

    // The predicted model structures, out of each engine's own output layout. AF2 is
    // the one engine whose choice of structure depends on run options, so they come
    // along (see foldNamingIsModelStructure).
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
    ch_complex_models = ch_predictions
        .flatMap { meta, files ->
            def opts = af2_opts + [keep_best: meta.af2_keep_models == 'best']
            foldNamingModelStructures(meta.fold_tool, files, opts).collect { model ->
                tuple(
                    [
                        id: meta.binder,
                        fold_tool: meta.fold_tool,
                        pair_id: meta.id,
                        target_id: label_target ? meta.target : null,
                        fold_batch: meta.fold_batch ?: 1,
                        fold_sample: foldNamingSampleLabel(meta.fold_tool, model),
                    ],
                    [model],
                )
            }
        }

    // Pair each engine's models back up with the reference they came from. combine(by: 0)
    // rather than join(): one reference feeds as many rows as there are engines (and
    // targets), and join would consume it after the first.
    ch_reference_by_id = ch_references
        .map { ref_id, structure, binder_chain, target_chain ->
            tuple(ref_id.toString(), structure, binder_chain, target_chain)
        }

    ch_compare_in = ch_complex_models
        .map { meta, models -> tuple(meta.id.toString(), meta, models) }
        .combine(ch_reference_by_id, by: 0)
        .map { _id, meta, models, structure, binder_chain, target_chain ->
            tuple(meta + [binder_chain: binder_chain, target_chain: target_chain], structure, models)
        }

    REFOLD_COMPARE_COMPLEX(ch_compare_in, 'B', 'A')

    ch_table = REFOLD_COMPARE_COMPLEX.out.rmsd
        .map { _meta, tsv -> tsv }
        .collectFile(
            name: table_name,
            storeDir: "${params.outdir}/${publish_dir}",
            keepHeader: true,
            skip: 1,
        )

    emit:
    rmsd = REFOLD_COMPARE_COMPLEX.out.rmsd
    models = ch_complex_models
    table = ch_table
}
