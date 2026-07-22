/*
ROSETTAFOLD3_FOLD: generic RF3 folding for fold.nf, decoupled from the
binder-design-coupled ROSETTAFOLD3 (modules/local/rfd3/rosettafold3.nf).

Phase 1 scope: monomer only (GENERATE_RF3_FOLD_INPUT / RF3_FOLD are written
generically for N chains, but fold.nf itself rejects multi-record FASTA
before this subworkflow ever runs - see fold.nf's monomer-only guard).
*/

include { GENERATE_RF3_FOLD_INPUT } from '../../modules/fold/rf3/generate_rf3_fold_input'
include { RF3_FOLD } from '../../modules/fold/rf3/rf3_fold'

// See boltz_fold.nf - same --n_predictions / --*_batch_size split semantics.
def foldPredictionBatches(batch_size_param, int default_batch, n_predictions) {
    def bs = (batch_size_param != null && !(batch_size_param instanceof Boolean)) \
        ? (batch_size_param as int) : null
    if (n_predictions) {
        def n = n_predictions as int
        def chunk = bs != null ? bs : n
        def n_batches = ((n + chunk - 1).intdiv(chunk)) as int
        return (0..<n_batches).collect { i ->
            def remaining = n - (i * chunk)
            remaining < chunk ? remaining : chunk
        }
    }
    return [bs != null ? bs : default_batch]
}

workflow ROSETTAFOLD3_FOLD {
    take:
    ch_for_rf3 // tuple(meta, fasta, a3m)

    main:
    GENERATE_RF3_FOLD_INPUT(ch_for_rf3)

    def batches = foldPredictionBatches(params.rf3_batch_size, 5, params.n_predictions)
    def namespaced = batches.size() > 1
    def base_seed = params.rf3_seed ? (params.rf3_seed as int) : null

    ch_batched = GENERATE_RF3_FOLD_INPUT.out.with_json.flatMap { meta, fasta, a3m, json ->
        batches.withIndex().collect { n_samples, i ->
            def m = meta + [
                fold_batch: i + 1,
                fold_batch_size: n_samples,
                fold_namespaced: namespaced,
            ]
            if (base_seed != null) {
                m = m + [rf3_seed: base_seed + i]
            }
            [m, fasta, a3m, json]
        }
    }

    RF3_FOLD(ch_batched)

    emit:
    predictions = RF3_FOLD.out.predictions
    confidence_json = RF3_FOLD.out.confidence_json
}
