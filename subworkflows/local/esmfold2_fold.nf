/*
ESMFOLD2_FOLD: ESMFold2 folding for fold.nf / fold_pulldown.nf (--methods esmfold2,
and --methods esmfold2_fast via the ESMFOLD2_FAST_FOLD alias with tool 'esmfold2_fast').

Consumes the same MSA bundle as PROTENIX_FOLD / ALPHAFOLD3_FOLD / OPENFOLD3_FOLD -
monomer (meta, fasta, a3m); multimer (meta, fasta, [paired...+unpaired...]).
Multimers go through GENERATE_ESMFOLD2_INPUT_COMPLEX to get key=<taxid> headers;
monomers need no pairing, so their a3m is passed straight through.
esmfold2_fast (biohub/ESMFold2-Fast) has no MSA encoder: FOLD_MSA hands it a
placeholder instead of an MSA and it always folds from sequence alone. Each tool
gets its own process alias so the configs can give them different containers.

ESMFold2's seed is a plain fold() argument, so batch i gets seed base+i
(base = --esmfold2_seeds, else 42) to keep batches distinct while -resume hashes
stay stable.
*/

include { GENERATE_ESMFOLD2_INPUT_COMPLEX } from '../../modules/local/fold/esmfold2/generate_esmfold2_input_complex'
include { ESMFOLD2 as ESMFOLD2_PROCESS } from '../../modules/local/fold/esmfold2/esmfold2'
include { ESMFOLD2 as ESMFOLD2_FAST_PROCESS } from '../../modules/local/fold/esmfold2/esmfold2'
include { FOLD_PARSE_CONFIDENCE } from '../../modules/local/fold/common/fold_parse_confidence'

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

workflow ESMFOLD2_FOLD {
    take:
    ch_for_esmfold2 // monomer: tuple(meta, fasta, a3m); multimer: tuple(meta, fasta, [paired...+unpaired...])
    tool            // 'esmfold2' | 'esmfold2_fast'

    main:
    def batches = foldPredictionBatches(params.esmfold2_batch_size, 5, params.n_predictions)
    def base_seed = params.esmfold2_seeds ? (params.esmfold2_seeds.toString().split(',')[0].trim() as int) : 42

    def is_fast = tool == 'esmfold2_fast'
    def single_sequence = is_fast || params.esmfold2_single_sequence
    def hf_model = is_fast ? 'biohub/ESMFold2-Fast' : 'biohub/ESMFold2'

    // Only multimers with real MSAs need the key= pairing render; in
    // single-sequence mode FOLD_MSA hands over a placeholder instead of per-chain
    // a3ms, so the generator would find nothing to render and fail on its output.
    ch_needs_pairing = ch_for_esmfold2.branch { meta, _fasta, _msa ->
        pair: (meta.n_chains ?: 1) > 1 && !single_sequence
        passthrough: true
    }
    ch_passthrough = ch_needs_pairing.passthrough
    GENERATE_ESMFOLD2_INPUT_COMPLEX(ch_needs_pairing.pair)
    ch_with_msa = ch_passthrough.mix(GENERATE_ESMFOLD2_INPUT_COMPLEX.out.with_msa)

    ch_batched = ch_with_msa.flatMap { meta, fasta, a3ms ->
        def is_mono = (meta.n_chains ?: 1) == 1
        def main_a3m = a3ms instanceof List ? a3ms[0] : a3ms
        def subsample = !single_sequence && MsaSubsample.isEnabled(params.msa_subsample)
        def n_seq = (is_mono && subsample) ? MsaSubsample.countA3mSequences(main_a3m) : null
        // Single-sequence runs get a placeholder, not an MSA, so there is nothing to
        // subsample (and its zero row count would reject every depth).
        def depth_jobs = subsample \
            ? MsaSubsample.depthJobs(params.msa_subsample, params.msa_subsample_include_full, n_seq) \
            : [null]
        def namespaced = batches.size() > 1 || depth_jobs.size() > 1
        def jobs = []
        batches.withIndex().each { n_samples, i ->
            depth_jobs.each { depth ->
                def m = meta + [
                    fold_batch: i + 1,
                    fold_batch_size: n_samples,
                    fold_namespaced: namespaced,
                    esmfold2_seed: base_seed + i,
                    esmfold2_tool: tool,
                    esmfold2_model: hf_model,
                    esmfold2_single_sequence: single_sequence,
                ]
                if (depth != null) {
                    def s = MsaSubsample.stableSeed(meta.id.toString(), i + 1, depth[0], depth[1])
                    m = m + [
                        msa_max_seq: depth[0],
                        msa_max_extra_seq: depth[1],
                        msa_subsample_seed: s,
                        msa_depth_tag: MsaSubsample.depthTag(depth[0], depth[1]),
                    ]
                }
                else if (subsample) {
                    m = m + [msa_depth_tag: 'full']
                }
                jobs << [m, fasta, a3ms]
            }
        }
        jobs
    }

    if (is_fast) {
        ESMFOLD2_FAST_PROCESS(ch_batched)
        ch_predictions = ESMFOLD2_FAST_PROCESS.out.predictions
        ch_confidence_json = ESMFOLD2_FAST_PROCESS.out.confidence_json
    }
    else {
        ESMFOLD2_PROCESS(ch_batched)
        ch_predictions = ESMFOLD2_PROCESS.out.predictions
        ch_confidence_json = ESMFOLD2_PROCESS.out.confidence_json
    }

    // One score row per ESMFold2 sample (seed_S_sample_N). run_esmfold2.py writes
    // AF3-shaped confidence files, so ipSAE and PAE come from the sibling
    // *_confidences.json; see protenix_fold.nf for why missing files fall back to
    // distinct dummy names.
    def no_ipsae_pae = file("${projectDir}/assets/dummy_files/no_ipsae_pae.json")
    def no_ipsae_struct = file("${projectDir}/assets/dummy_files/no_ipsae_structure.cif")

    ch_conf = ch_predictions.flatMap { meta, files ->
        def all = files instanceof List ? files : [files]
        def byName = [:]
        all.each { f -> byName[f.name] = f }
        all.findAll { it.name ==~ /.*_seed_\d+_sample_\d+_summary_confidences\.json$/ }
            .collect { j ->
                def mm = (j.name =~ /_(seed_\d+_sample_\d+)_summary_confidences\.json$/)
                def model = mm ? mm[0][1] : 'sample'
                def struct = j.name.replaceFirst(/_summary_confidences\.json$/, '_model.cif')
                def paeName = j.name.replaceFirst(/_summary_confidences\.json$/, '_confidences.json')
                def pae = byName[paeName]
                def cif = byName[struct]
                def doIpsae = (pae != null && cif != null)
                def pred = "${FoldNaming.flatPrefix(tool, meta)}${struct}"
                [meta, tool, model, struct, pred, j, doIpsae ? pae : no_ipsae_pae, doIpsae ? cif : no_ipsae_struct, doIpsae]
            }
    }
    FOLD_PARSE_CONFIDENCE(ch_conf)
    ch_tsv = FOLD_PARSE_CONFIDENCE.out.collectFile(
        name: "${tool}_fold_scores.tsv",
        storeDir: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/${tool}",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = ch_predictions
    confidence_json = ch_confidence_json
    tsv = ch_tsv
}
