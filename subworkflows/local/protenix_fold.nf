/*
PROTENIX_FOLD: generic Protenix (AF3-style) folding for fold.nf, the 4th
--methods engine (see plans/fold-nf-multi-method-folding.md Phase 3).

Monomer inputs (meta.n_chains == 1) use GENERATE_PROTENIX_INPUT with a single
unpairedMsaPath a3m; multimer inputs (n_chains > 1) use
GENERATE_PROTENIX_INPUT_COMPLEX with per-chain paired + unpaired a3m bundles
from FOLD_MSA (Protenix pairs by species mnemonic; see
plans/fold-nf-multimer-paired-msa.md). Both feed the same PROTENIX_FOLD predict.
*/

include { GENERATE_PROTENIX_INPUT } from '../../modules/fold/protenix/generate_protenix_input'
include { GENERATE_PROTENIX_INPUT_COMPLEX } from '../../modules/fold/protenix/generate_protenix_input_complex'
include { PROTENIX_FOLD as PROTENIX_FOLD_PROCESS } from '../../modules/fold/protenix/protenix_fold'
include { FOLD_PARSE_CONFIDENCE } from '../../modules/fold/common/fold_parse_confidence'

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
    ch_for_protenix // monomer: tuple(meta, fasta, a3m); multimer: tuple(meta, fasta, [paired...+unpaired...])

    main:
    ch_mono = ch_for_protenix.filter { meta, fasta, msa -> (meta.n_chains ?: 1) == 1 }
    ch_multi = ch_for_protenix.filter { meta, fasta, msa -> (meta.n_chains ?: 1) > 1 }
    GENERATE_PROTENIX_INPUT(ch_mono)
    GENERATE_PROTENIX_INPUT_COMPLEX(ch_multi)
    ch_with_json = GENERATE_PROTENIX_INPUT.out.with_json.mix(GENERATE_PROTENIX_INPUT_COMPLEX.out.with_json)

    def batches = foldPredictionBatches(params.protenix_batch_size, 5, params.n_predictions)

    if (params.protenix_seeds && params.protenix_seeds.toString().contains(',')) {
        // namespaced is only known after per-a3m depth filtering; reject multi-seed
        // whenever batching or subsample could fan out (same rule as before).
        def maybe_depths = MsaSubsample.depthJobs(params.msa_subsample, params.msa_subsample_include_full)
        if (batches.size() > 1 || maybe_depths.size() > 1) {
            error(
                "fold.nf: --protenix_seeds with multiple comma-separated values cannot be " +
                "auto-offset across --protenix_batch_size jobs; pin a single seed (batches " +
                "use seed, seed+1, ...) or leave --protenix_seeds unset."
            )
        }
    }

    def base_seed = params.protenix_seeds ? (params.protenix_seeds.toString().split(',')[0].trim() as int) : null

    ch_batched = ch_with_json.flatMap { meta, fasta, a3m, json ->
        def n_seq = MsaSubsample.isEnabled(params.msa_subsample) \
            ? MsaSubsample.countA3mSequences(a3m) : null
        def depth_jobs = MsaSubsample.depthJobs(
            params.msa_subsample, params.msa_subsample_include_full, n_seq
        )
        def namespaced = batches.size() > 1 || depth_jobs.size() > 1
        def jobs = []
        batches.withIndex().each { n_samples, i ->
            depth_jobs.each { depth ->
                def m = meta + [
                    fold_batch: i + 1,
                    fold_batch_size: n_samples,
                    fold_namespaced: namespaced,
                ]
                if (base_seed != null) {
                    m = m + [protenix_seeds: "${base_seed + i}"]
                }
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
                jobs << [m, fasta, a3m, json]
            }
        }
        jobs
    }

    PROTENIX_FOLD_PROCESS(ch_batched)

    // One score row per Protenix sample. model/struct/predictions_file are
    // derived from the summary filename (complex_summary_confidence_sample_N.json
    // -> complex_sample_N.cif); predictions_file uses the same FoldNaming prefix
    // as the module's flat-gather saveAs.
    ch_conf = PROTENIX_FOLD_PROCESS.out.confidence_json.flatMap { meta, jsons ->
        def files = (jsons instanceof List) ? jsons : [jsons]
        files.collect { j ->
            def mm = (j.name =~ /_summary_confidence_sample_(\d+)\.json$/)
            def idx = mm ? mm[0][1] : '0'
            def struct = j.name.replaceFirst(/_summary_confidence_sample_(\d+)\.json$/, '_sample_$1.cif')
            def pred = "${FoldNaming.flatPrefix('protenix', meta)}${struct}"
            [meta, 'protenix', "sample_${idx}", struct, pred, j]
        }
    }
    FOLD_PARSE_CONFIDENCE(ch_conf)
    ch_tsv = FOLD_PARSE_CONFIDENCE.out.collectFile(
        name: 'protenix_fold_scores.tsv',
        storeDir: "${params.outdir}/fold/protenix",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = PROTENIX_FOLD_PROCESS.out.predictions
    confidence_json = PROTENIX_FOLD_PROCESS.out.confidence_json
    tsv = ch_tsv
}
