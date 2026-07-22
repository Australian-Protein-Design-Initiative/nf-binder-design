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
    // AF2 always generates its 5 trained models per predict. --af2_keep_models
    // chooses which we retain toward --n_predictions (no free --af2_batch_size):
    //   all : keep all 5/run -> ceil(n_predictions / 5) runs
    //   best: keep only the top-ranked/run -> n_predictions runs
    def keep_models = params.af2_keep_models
    def kept_per_run = (keep_models == 'best') ? 1 : 5
    def n_pred = params.n_predictions ? (params.n_predictions as int) : 0
    def n_runs = n_pred ? ((n_pred + kept_per_run - 1).intdiv(kept_per_run)) : 1

    // Seeds. Pinned --af2_random_seed fans out deterministically
    // (base, base+1, ...) for reproducibility. Unset -> leave each run's seed
    // out entirely so AF2 draws its own random seed AND the task command stays
    // seed-free (drawing a random int here would bust -resume on every re-parse).
    // Runs stay distinct via meta.af2_run, which is part of the task input hash,
    // so omitting the seed does not collapse the N runs into one.
    def base = params.af2_random_seed ? (params.af2_random_seed as int) : null
    def seeds = (0..<n_runs).collect { i -> base != null ? base + i : null }
    def namespaced = n_runs > 1

    ch_runs = ch_af2_msas.flatMap { meta, fasta, msas ->
        seeds.withIndex().collect { seed, i ->
            def m = meta + [af2_run: i + 1, af2_keep_models: keep_models, af2_namespaced: namespaced]
            if (seed != null) { m = m + [af2_seed: seed] }
            [m, fasta, msas]
        }
    }

    ALPHAFOLD2_PREDICT(ch_runs)

    emit:
    predictions = ALPHAFOLD2_PREDICT.out.predictions
}
