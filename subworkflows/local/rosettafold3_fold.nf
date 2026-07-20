/*
ROSETTAFOLD3_FOLD: generic RF3 folding for fold.nf, decoupled from the
binder-design-coupled ROSETTAFOLD3 (modules/local/rfd3/rosettafold3.nf).

Phase 1 scope: monomer only (GENERATE_RF3_FOLD_INPUT / RF3_FOLD are written
generically for N chains, but fold.nf itself rejects multi-record FASTA
before this subworkflow ever runs - see fold.nf's monomer-only guard).
*/

include { GENERATE_RF3_FOLD_INPUT } from '../../modules/fold/rf3/generate_rf3_fold_input'
include { RF3_FOLD } from '../../modules/fold/rf3/rf3_fold'

workflow ROSETTAFOLD3_FOLD {
    take:
    ch_for_rf3 // tuple(meta, fasta, a3m)

    main:
    GENERATE_RF3_FOLD_INPUT(ch_for_rf3)
    RF3_FOLD(GENERATE_RF3_FOLD_INPUT.out.with_json)

    emit:
    predictions = RF3_FOLD.out.predictions
    confidence_json = RF3_FOLD.out.confidence_json
}
