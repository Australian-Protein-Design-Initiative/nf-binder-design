/*
Thin subworkflow wrapper around the GPU ALPHAFOLD2 predict module, keeping the
predict-stage wiring separate from MSA generation. MSA generation is
deliberately NOT part of this subworkflow: fold.nf's FOLD_MSA
(subworkflows/local/fold_msa.nf) feeds this stage either the native
jackhmmer/hhblits MSAs (--msa_method jackhmmer_af2) or a ColabFold a3m bridged
into the same precomputed-MSA directory contract this subworkflow's `take`
expects (--msa_method mmseqs2_colabfold).
*/

include { ALPHAFOLD2 as ALPHAFOLD2_PREDICT } from '../../modules/fold/af2/alphafold2'

workflow ALPHAFOLD2 {
    take:
    // tuple(meta, fasta, msas_dir) - msas_dir MUST be a directory named
    // exactly meta.id containing features.pkl (the only file the predict
    // module actually reads under --use_precomputed_msas=true; see
    // modules/fold/af2/colabfold_a3m_to_af2_msas.nf for why).
    ch_af2_msas

    main:
    ALPHAFOLD2_PREDICT(ch_af2_msas)

    emit:
    predictions = ALPHAFOLD2_PREDICT.out.predictions
}
