#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
AF2IG_SCORING: the primary cofold's score table for rfd and rfd_partial.

AF2 initial guess is these workflows' own cofold - every design goes through it - and
this adds BindCraft's interface metrics on top, then merges the lot into
combined_scores.tsv.

It used to be BOLTZ_REFOLD_SCORING, which did this and refolding both, choosing between
them on whether --refold_af2ig_filters was set. Refolding is REFOLD's job now, in every
workflow and with any engine, so what is left here is the part that always runs. The
refold's own table is merged in when there is one, which is what ch_refold_scores
carries: an empty placeholder file when --refold_with names no engines.
*/

include { COMBINE_SCORES } from '../../modules/local/common/combine_scores'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_AF2IG } from '../../modules/local/rfd/bindcraft_scoring'

workflow AF2IG_SCORING {

    take:
    ch_af2ig_pdbs      // channel of PDB files
    ch_af2ig_scores    // channel of AF2 initial guess score files
    binder_chain       // String: binder chain id in the AF2IG complex
    ch_refold_scores   // REFOLD's wide table, or an empty file when nothing was refolded
    outdir             // String: params.outdir

    main:

    BINDCRAFT_SCORING_AF2IG(
        ch_af2ig_pdbs,
        binder_chain,
        'default_4stage_multimer',
        'rfd/af2_initial_guess/extra_scores/',
    )

    def ch_extra_scores = BINDCRAFT_SCORING_AF2IG.out.scores.collectFile(
        name: "${outdir}/rfd/af2_initial_guess/af2ig_extra_scores.tsv",
        keepHeader: true,
        skip: 1,
    )

    COMBINE_SCORES(
        ch_af2ig_scores.collect(),
        ch_extra_scores,
        ch_refold_scores,
        ch_af2ig_pdbs.collect(),
    )

    emit:
    scores = COMBINE_SCORES.out.combined_scores
    af2_initial_guess_scores = COMBINE_SCORES.out.af2_initial_guess_scores
    extra_scores = COMBINE_SCORES.out.extra_scores
    shape_scores = COMBINE_SCORES.out.shape_scores
    binders_fasta = COMBINE_SCORES.out.binders_fasta
}
