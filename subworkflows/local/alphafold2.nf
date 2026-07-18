/*
Thin subworkflow wrapper around the existing GPU ALPHAFOLD2 predict module, so
af2.nf and fold.nf share identical predict-stage wiring. MSA generation is
deliberately NOT part of this subworkflow: af2.nf always uses the native
jackhmmer/hhblits stage (ALPHAFOLD2_JACKHMMER_MSA), while fold.nf's FOLD_MSA
(subworkflows/local/fold_msa.nf) can also bridge a ColabFold a3m into the same
precomputed-MSA directory contract this subworkflow's `take` expects.
*/

include { ALPHAFOLD2 as ALPHAFOLD2_PREDICT } from '../../modules/local/af2/alphafold2'

workflow ALPHAFOLD2 {
    take:
    // tuple(meta, fasta, msas_dir) - msas_dir MUST be a directory named
    // exactly meta.id containing features.pkl (the only file the predict
    // module actually reads under --use_precomputed_msas=true; see
    // modules/local/af2/colabfold_a3m_to_af2_msas.nf for why).
    ch_af2_msas

    main:
    ALPHAFOLD2_PREDICT(ch_af2_msas)

    emit:
    predictions = ALPHAFOLD2_PREDICT.out.predictions
}
