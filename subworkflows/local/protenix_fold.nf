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

workflow PROTENIX_FOLD {
    take:
    ch_for_protenix // tuple(meta, fasta, a3m)

    main:
    GENERATE_PROTENIX_INPUT(ch_for_protenix)
    PROTENIX_FOLD_PROCESS(GENERATE_PROTENIX_INPUT.out.with_json)

    emit:
    predictions = PROTENIX_FOLD_PROCESS.out.predictions
    confidence_json = PROTENIX_FOLD_PROCESS.out.confidence_json
}
