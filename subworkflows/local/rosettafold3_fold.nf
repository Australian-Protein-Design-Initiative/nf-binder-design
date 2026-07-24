/*
ROSETTAFOLD3_FOLD: generic RF3 folding for fold.nf, decoupled from the
binder-design-coupled ROSETTAFOLD3 (modules/local/rfd3/rosettafold3.nf).

Monomer inputs (meta.n_chains == 1) use GENERATE_RF3_FOLD_INPUT with the shared
single a3m; multimer inputs (n_chains > 1) use GENERATE_RF3_FOLD_INPUT_COMPLEX
with the per-chain TaxID=-annotated a3m bundle from FOLD_MSA (see
plans/fold-nf-multimer-paired-msa.md). Both feed the same RF3_FOLD predict.
*/

include { GENERATE_RF3_FOLD_INPUT } from '../../modules/fold/rf3/generate_rf3_fold_input'
include { GENERATE_RF3_FOLD_INPUT_COMPLEX } from '../../modules/fold/rf3/generate_rf3_fold_input_complex'
include { RF3_FOLD } from '../../modules/fold/rf3/rf3_fold'
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

workflow ROSETTAFOLD3_FOLD {
    take:
    ch_for_rf3 // monomer: tuple(meta, fasta, a3m); multimer: tuple(meta, fasta, [a3m...])

    main:
    ch_mono = ch_for_rf3.filter { meta, fasta, msa -> (meta.n_chains ?: 1) == 1 }
    ch_multi = ch_for_rf3.filter { meta, fasta, msa -> (meta.n_chains ?: 1) > 1 }
    GENERATE_RF3_FOLD_INPUT(ch_mono)
    GENERATE_RF3_FOLD_INPUT_COMPLEX(ch_multi)
    ch_with_json = GENERATE_RF3_FOLD_INPUT.out.with_json.mix(GENERATE_RF3_FOLD_INPUT_COMPLEX.out.with_json)

    def batches = foldPredictionBatches(params.rf3_batch_size, 5, params.n_predictions)
    def base_seed = params.rf3_seed ? (params.rf3_seed as int) : null

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
                    m = m + [rf3_seed: base_seed + i]
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

    RF3_FOLD(ch_batched)

    // One score row per diffusion sample. RF3 also writes a top-level
    // best-model *_summary_confidences.json with no sample index - skip it (it
    // isn't gathered into fold/predictions/). model/struct/predictions_file are
    // derived from the summary filename; predictions_file uses the same
    // FoldNaming prefix as the module's flat-gather saveAs.
    ch_conf = RF3_FOLD.out.confidence_json.flatMap { meta, jsons ->
        def files = (jsons instanceof List) ? jsons : [jsons]
        files.findAll { it.name ==~ /.*_seed-\d+_sample-\d+_summary_confidences\.json/ }
            .collect { j ->
                def mm = (j.name =~ /_(seed-\d+_sample-\d+)_summary_confidences\.json$/)
                def model = mm ? mm[0][1] : j.baseName
                def struct = j.name.replaceFirst(/_summary_confidences\.json$/, '_model.cif')
                def pred = "${FoldNaming.flatPrefix('rf3', meta)}${struct}"
                [meta, 'rf3', model, struct, pred, j]
            }
    }
    FOLD_PARSE_CONFIDENCE(ch_conf)
    ch_tsv = FOLD_PARSE_CONFIDENCE.out.collectFile(
        name: 'rf3_fold_scores.tsv',
        storeDir: "${params.outdir}/fold/rf3",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = RF3_FOLD.out.predictions
    confidence_json = RF3_FOLD.out.confidence_json
    tsv = ch_tsv
}
