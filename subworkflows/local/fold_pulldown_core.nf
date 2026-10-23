/*
FOLD_PULLDOWN_CORE: co-fold every binder against every target with any set of
folding oracles, then summarise and rank interface scores.

This is the engine behind `--method fold_pulldown`, lifted out of that workflow so
other workflows can reuse it. It takes sequences as channels rather than reading
`params.targets` / `params.binders`, so a caller that has structures rather than
FASTA files (a binder design workflow, say) can feed it sequences it extracted
itself. See REFOLD_ORACLES for that case.

Oracle choice is entirely the caller's: whatever FOLD_PREDICT supports
(af2, af2_mono, boltz, rf3, protenix, af3, openfold3, chai1, esmfold2,
esmfold2_fast) can be passed in `methods`.

Publishing still goes to `${params.outdir}/${params.fold_publish_dir}`, because the
processes underneath read that param directly. A caller that wants its own
subdirectory sets `params.fold_publish_dir` in its own param block.
*/

include { FOLD_PULLDOWN_MSA } from './fold_pulldown_msa'
include { FOLD_PREDICT } from './fold_predict'
include { FOLD_TEMPLATES } from './fold_templates'
include { FOLD_PULLDOWN_MERGE_SCORES } from '../../modules/local/fold/common/fold_pulldown_merge_scores'

workflow FOLD_PULLDOWN_CORE {

    take:
    ch_targets      // tuple(meta, fasta)  meta: [id, seq]
    ch_binders      // tuple(meta, fasta)
    methods         // List<String>
    msa_method      // 'jackhmmer_af2' | 'mmseqs2_colabfold'
    ch_template_src // channel of paths to search for templates (may be empty)

    main:

    FOLD_PULLDOWN_MSA(ch_targets, ch_binders, methods, msa_method)

    FOLD_TEMPLATES(ch_template_src)

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
        FOLD_TEMPLATES.out.templates,
        methods,
    )

    FOLD_PULLDOWN_MERGE_SCORES(FOLD_PREDICT.out.scores, ch_pairs_tsv)

    emit:
    scores = FOLD_PULLDOWN_MERGE_SCORES.out.scores
    summary = FOLD_PULLDOWN_MERGE_SCORES.out.summary
    pairs = ch_pairs_tsv
    predictions = FOLD_PREDICT.out.predictions
}
