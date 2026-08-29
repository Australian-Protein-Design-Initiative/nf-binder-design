/*
Thin subworkflow wrapper around the GPU ALPHAFOLD2 predict module, keeping the
predict-stage wiring separate from MSA generation. MSA generation is
deliberately NOT part of this subworkflow: fold.nf's FOLD_MSA
(subworkflows/local/fold_msa.nf) feeds this stage either the native
jackhmmer/hhblits MSAs (--msa_method jackhmmer_af2) or a ColabFold a3m bridged
into the same precomputed-MSA directory contract this subworkflow's `take`
expects (--msa_method mmseqs2_colabfold).

When --msa_subsample is on, shallow depth jobs rebuild features.pkl from a
subsampled a3m (empty templates); full-depth jobs reuse the precomputed
features.pkl (templates kept for jackhmmer).
*/

include { ALPHAFOLD2 as ALPHAFOLD2_PREDICT } from '../../modules/fold/af2/alphafold2'
include { FOLD_SCORE_AF2 } from '../../modules/fold/af2/fold_score_af2'

workflow ALPHAFOLD2 {
    take:
    // tuple(meta, fasta, msas_dir, a3m)
    // a3m may be a dummy stub when --msa_subsample is off.
    ch_af2_input
    // 'af2' (multimer) or 'af2_mono' (monomer weights + residue_index chain break).
    // Both engines can run in the same pipeline; the tag keeps their published
    // predictions, score TSVs and FoldNaming prefixes apart.
    tool

    main:
    def complex_mode = (tool == 'af2_mono') ? 'chainbreak' : 'multimer'
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

    ch_runs = ch_af2_input.flatMap { meta, fasta, msas, a3m ->
        def n_seq = MsaSubsample.isEnabled(params.msa_subsample) \
            ? MsaSubsample.countA3mSequences(a3m) : null
        def depth_jobs = MsaSubsample.depthJobs(
            params.msa_subsample, params.msa_subsample_include_full, n_seq
        )
        def namespaced = n_runs > 1 || depth_jobs.size() > 1
        def jobs = []
        seeds.withIndex().each { seed, i ->
            depth_jobs.each { depth ->
                def m = meta + [af2_run: i + 1, af2_keep_models: keep_models,
                                af2_namespaced: namespaced,
                                af2_tool: tool, af2_complex_mode: complex_mode]
                if (seed != null) { m = m + [af2_seed: seed] }
                if (depth != null) {
                    def s = MsaSubsample.stableSeed(meta.id.toString(), i + 1, depth[0], depth[1])
                    m = m + [
                        msa_max_seq: depth[0],
                        msa_max_extra_seq: depth[1],
                        msa_subsample_seed: s,
                        msa_depth_tag: MsaSubsample.depthTag(depth[0], depth[1]),
                    ]
                }
                else if (MsaSubsample.isEnabled(params.msa_subsample)) {
                    m = m + [msa_depth_tag: 'full']
                }
                jobs << [m, fasta, msas, a3m]
            }
        }
        jobs
    }

    ALPHAFOLD2_PREDICT(ch_runs)

    // Score each run: FoldNaming.af2Prefix(meta) is the exact fold/predictions/
    // filename prefix the module's saveAs uses, so predictions_file lines up.
    // The predictions glob stages flat into the scoring task's work dir, so
    // FOLD_SCORE_AF2 reads it with --run-dir . (see fold_score_af2.nf).
    // Drop the MSA intermediates first: AF2 multimer writes per-chain
    // msas/<chain>/*.sto (+ .a3m) with identical basenames across chains, which
    // collide when the predictions glob stages flat. score_af2_run.py only needs
    // the top-level structures + ranking_debug.json / pae_model_*.json /
    // result_model_*.pkl, so excluding msas/ is safe (and tidier for monomer too).
    ch_score_in = ALPHAFOLD2_PREDICT.out.predictions
        .map { meta, files ->
            def scoring = (files instanceof List ? files : [files]).findAll {
                !it.toString().contains('/msas/')
            }
            [meta, scoring, FoldNaming.af2Prefix(meta)]
        }
    FOLD_SCORE_AF2(ch_score_in)

    ch_tsv = FOLD_SCORE_AF2.out.collectFile(
        name: "${tool}_fold_scores.tsv",
        storeDir: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/${tool}",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = ALPHAFOLD2_PREDICT.out.predictions
    tsv = ch_tsv
}
