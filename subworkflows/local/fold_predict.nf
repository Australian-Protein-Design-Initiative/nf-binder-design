/*
FOLD_PREDICT: run selected structure predictors and merge scores.

Takes per-engine MSA-ready channels from FOLD_MSA (or FOLD_PULLDOWN_MSA) and
dispatches ALPHAFOLD2 / BOLTZ_FOLD / ROSETTAFOLD3_FOLD / PROTENIX_FOLD, then
merges per-tool score TSVs into fold_scores.tsv. Optional EnGens clustering.
*/

include { ALPHAFOLD2 } from './alphafold2'
include { ALPHAFOLD2 as ALPHAFOLD2_MONO } from './alphafold2'
include { BOLTZ_FOLD } from './boltz_fold'
include { ROSETTAFOLD3_FOLD } from './rosettafold3_fold'
include { PROTENIX_FOLD } from './protenix_fold'
include { ENGENS_CLUSTER } from './engens'
include { FOLD_MERGE_SCORES } from '../../modules/fold/common/fold_merge_scores'

workflow FOLD_PREDICT {
    take:
    ch_af2_in      // tuple(meta, fasta, msas_dir, a3m) - may be empty if af2 not selected
    ch_for_boltz   // from FOLD_MSA / FOLD_PULLDOWN_MSA
    ch_for_rf3
    ch_for_protenix
    methods        // List<String>

    main:
    ch_af2_pred = Channel.empty()
    ch_af2_mono_pred = Channel.empty()
    ch_boltz_pred = Channel.empty()
    ch_rf3_pred = Channel.empty()
    ch_protenix_pred = Channel.empty()
    ch_scores = Channel.empty()

    if ('af2' in methods) {
        ALPHAFOLD2(ch_af2_in, 'af2')
        ch_af2_pred = ALPHAFOLD2.out.predictions
        ch_scores = ch_scores.mix(ALPHAFOLD2.out.tsv)
    }
    // af2_mono reuses the same per-chain MSA directories as af2 - it only assembles
    // them into features.pkl differently (block diagonal, one chain-break jump).
    if ('af2_mono' in methods) {
        ALPHAFOLD2_MONO(ch_af2_in, 'af2_mono')
        ch_af2_mono_pred = ALPHAFOLD2_MONO.out.predictions
        ch_scores = ch_scores.mix(ALPHAFOLD2_MONO.out.tsv)
    }
    if ('boltz' in methods) {
        BOLTZ_FOLD(ch_for_boltz)
        ch_boltz_pred = BOLTZ_FOLD.out.predictions
        ch_scores = ch_scores.mix(BOLTZ_FOLD.out.tsv)
    }
    if ('rf3' in methods) {
        ROSETTAFOLD3_FOLD(ch_for_rf3)
        ch_rf3_pred = ROSETTAFOLD3_FOLD.out.predictions
        ch_scores = ch_scores.mix(ROSETTAFOLD3_FOLD.out.tsv)
    }
    if ('protenix' in methods) {
        PROTENIX_FOLD(ch_for_protenix)
        ch_protenix_pred = PROTENIX_FOLD.out.predictions
        ch_scores = ch_scores.mix(PROTENIX_FOLD.out.tsv)
    }

    FOLD_MERGE_SCORES(ch_scores.collect())

    ch_predictions = ch_af2_pred.mix(ch_af2_mono_pred, ch_boltz_pred, ch_rf3_pred, ch_protenix_pred)

    if (!params.skip_engens) {
        ENGENS_CLUSTER(ch_predictions)
    }

    emit:
    scores = FOLD_MERGE_SCORES.out.tsv
    predictions = ch_predictions
}
