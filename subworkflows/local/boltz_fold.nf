/*
BOLTZ_FOLD: generic single-target Boltz-2 folding for fold.nf.

Phase 1 scope: monomer only. Multimer (one `protein:` entry per FASTA record,
each chain's own `msa:`) is Phase 2 work (see plans/
fold-nf-multi-method-folding.md §5).
*/

include { FOLD_CREATE_BOLTZ_YAML } from '../../modules/fold/boltz/fold_create_boltz_yaml'
include { BOLTZ } from '../../modules/local/common/boltz'
include { FOLD_PARSE_BOLTZ_CONFIDENCE } from '../../modules/fold/boltz/fold_parse_boltz_confidence'

// Split --n_predictions across jobs of at most batch_size samples each.
// batch_size unset (Boolean false / null) + n_predictions set => one job of N.
// Neither set => one job with default_batch samples.
// Note: do not use Groovy truthiness on batch_size (0 == false).
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

workflow BOLTZ_FOLD {
    take:
    ch_for_boltz // tuple(meta, fasta, a3m)

    main:
    FOLD_CREATE_BOLTZ_YAML(ch_for_boltz)

    ch_templates = params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates")

    // BOLTZ's process signature is shared with boltz_pulldown.nf's
    // target+binder complex mode, so it always expects two MSA paths to
    // stage (their names never appear in the script body - they only need to
    // be physically present alongside the YAML, since the YAML's `msa:`
    // field is what boltz predict actually reads). fold.nf has one MSA, so
    // the unused "target" slot gets the same empty-placeholder file
    // boltz_pulldown.nf uses for its own target-msa-less branches.
    def batches = foldPredictionBatches(params.boltz_batch_size, 1, params.n_predictions)
    def namespaced = batches.size() > 1
    def base_seed = params.boltz_seed ? (params.boltz_seed as int) : null

    ch_boltz_input = FOLD_CREATE_BOLTZ_YAML.out.yaml.flatMap { meta, yaml, msa ->
        batches.withIndex().collect { n_samples, i ->
            def m = meta + [
                fold_batch: i + 1,
                fold_batch_size: n_samples,
                fold_namespaced: namespaced,
            ]
            if (base_seed != null) {
                m = m + [boltz_seed: base_seed + i]
            }
            [m, yaml, file("${projectDir}/assets/dummy_files/empty_target_msa"), msa]
        }
    }

    // step_name is BOLTZ's publishDir subdir under params.outdir; use 'fold/boltz'
    // so predictions land at <outdir>/fold/boltz/ alongside boltz_fold_scores.tsv
    // (matching af2's <outdir>/fold/af2/ and rf3's <outdir>/fold/rf3/ layout).
    // fold.nf always requests mmCIF so fold/predictions/ stays format-uniform
    // with af2/rf3/protenix (boltz_pulldown hardcodes pdb).
    BOLTZ(ch_boltz_input, ch_templates, 'fold/boltz', 'cif')

    // Fan out to one row per diffusion sample so all models are scored (not
    // just the top-ranked model_0). The all-samples emit is a list when >1
    // sample, a single path when 1; normalise, then tag each with its model
    // index parsed from the confidence_<id>_model_<N>.json filename.
    ch_conf_per_model = BOLTZ.out.confidence_json_all
        .flatMap { meta, jsons ->
            def files = (jsons instanceof List) ? jsons : [jsons]
            files.collect { j ->
                def m = (j.name =~ /_model_(\d+)\.json$/)
                def idx = m ? m[0][1] : '0'
                [meta, idx, j]
            }
        }

    FOLD_PARSE_BOLTZ_CONFIDENCE(ch_conf_per_model)

    ch_tsv = FOLD_PARSE_BOLTZ_CONFIDENCE.out.collectFile(
        name: 'boltz_fold_scores.tsv',
        storeDir: "${params.outdir}/fold/boltz",
        keepHeader: true,
        skip: 1,
    )

    emit:
    // All per-sample structures (model_0..N-1) for Engens / fold/predictions
    // gather; not just model_0.
    predictions = BOLTZ.out.structure_all
    confidence_json = BOLTZ.out.confidence_json
    tsv = ch_tsv
}
