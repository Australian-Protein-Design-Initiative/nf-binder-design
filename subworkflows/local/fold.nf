/*
FOLD: shared multi-method structure folding for --method fold.

Runs FOLD_MSA then FOLD_PREDICT. Used by workflows/fold.nf; fold_pulldown
builds its own MSA channel and calls FOLD_PREDICT directly.
*/

include { FOLD_MSA } from './fold_msa'
include { FOLD_PREDICT } from './fold_predict'

workflow FOLD {
    take:
    ch_input   // tuple(meta, fasta)
    methods    // List<String>
    msa_method // 'jackhmmer_af2' | 'mmseqs2_colabfold'

    main:
    FOLD_MSA(ch_input, methods, msa_method)

    def a3m_stub = file("${projectDir}/assets/dummy_files/empty")
    if (('af2' in methods) || ('af2_mono' in methods)) {
        if (MsaSubsample.isEnabled(params.msa_subsample)) {
            ch_af2_in = FOLD_MSA.out.af2_msas
                .join(FOLD_MSA.out.a3m.map { meta, fasta, a3m -> [meta, a3m] })
                .map { meta, fasta, msas, a3m -> [meta, fasta, msas, a3m] }
        }
        else {
            ch_af2_in = FOLD_MSA.out.af2_msas.map { meta, fasta, msas -> [meta, fasta, msas, a3m_stub] }
        }
    }
    else {
        ch_af2_in = Channel.empty()
    }

    FOLD_PREDICT(
        ch_af2_in,
        FOLD_MSA.out.for_boltz,
        FOLD_MSA.out.for_rf3,
        FOLD_MSA.out.for_protenix,
        methods,
    )

    emit:
    scores = FOLD_PREDICT.out.scores
    predictions = FOLD_PREDICT.out.predictions
}
