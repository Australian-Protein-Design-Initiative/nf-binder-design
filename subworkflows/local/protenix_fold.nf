/*
PROTENIX_FOLD: generic Protenix (AF3-style) folding for fold.nf, the 4th
--methods engine (see plans/fold-nf-multi-method-folding.md Phase 3).

Phase 3 scope (still monomer only, consistent with Phase 1): GENERATE_PROTENIX_INPUT /
PROTENIX_FOLD are written N-chain-generic, but fold.nf itself rejects
multi-record FASTA before this subworkflow ever runs - see fold.nf's
monomer-only guard.
*/

include { GENERATE_PROTENIX_INPUT } from '../../modules/fold/protenix/generate_protenix_input'
include { PROTENIX_FOLD as PROTENIX_FOLD_PROCESS } from '../../modules/fold/protenix/protenix_fold'

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

workflow PROTENIX_FOLD {
    take:
    ch_for_protenix // tuple(meta, fasta, a3m)

    main:
    GENERATE_PROTENIX_INPUT(ch_for_protenix)

    def batches = foldPredictionBatches(params.protenix_batch_size, 5, params.n_predictions)
    def namespaced = batches.size() > 1

    if (namespaced && params.protenix_seeds && params.protenix_seeds.toString().contains(',')) {
        error(
            "fold.nf: --protenix_seeds with multiple comma-separated values cannot be " +
            "auto-offset across --protenix_batch_size jobs; pin a single seed (batches " +
            "use seed, seed+1, ...) or leave --protenix_seeds unset."
        )
    }

    def base_seed = params.protenix_seeds ? (params.protenix_seeds.toString().split(',')[0].trim() as int) : null

    ch_batched = GENERATE_PROTENIX_INPUT.out.with_json.flatMap { meta, fasta, a3m, json ->
        batches.withIndex().collect { n_samples, i ->
            def m = meta + [
                fold_batch: i + 1,
                fold_batch_size: n_samples,
                fold_namespaced: namespaced,
            ]
            if (base_seed != null) {
                m = m + [protenix_seeds: "${base_seed + i}"]
            }
            [m, fasta, a3m, json]
        }
    }

    PROTENIX_FOLD_PROCESS(ch_batched)

    emit:
    predictions = PROTENIX_FOLD_PROCESS.out.predictions
    confidence_json = PROTENIX_FOLD_PROCESS.out.confidence_json
}
