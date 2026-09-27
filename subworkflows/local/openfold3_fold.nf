/*
OPENFOLD3_FOLD: OpenFold3 folding for fold.nf / fold_pulldown.nf (--methods openfold3).

Consumes the same MSA bundle as PROTENIX_FOLD / ALPHAFOLD3_FOLD - monomer
(meta, fasta, a3m); multimer (meta, fasta, [paired...+unpaired...]). Seeds are
not part of OpenFold3's query JSON, so one JSON is generated per input and the
batch / MSA-depth fan-out happens afterwards. OpenFold3's default seed is a
fixed 42, so batch i gets seed base+i (base = --openfold3_seeds, else 42) to keep
batches distinct while -resume hashes stay stable.
*/

include { GENERATE_OPENFOLD3_INPUT } from '../../modules/fold/openfold3/generate_openfold3_input'
include { GENERATE_OPENFOLD3_INPUT_COMPLEX } from '../../modules/fold/openfold3/generate_openfold3_input_complex'
include { OPENFOLD3 as OPENFOLD3_PROCESS } from '../../modules/fold/openfold3/openfold3'
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

workflow OPENFOLD3_FOLD {
    take:
    ch_for_openfold3 // monomer: tuple(meta, fasta, a3m); multimer: tuple(meta, fasta, [paired...+unpaired...])

    main:
    def batches = foldPredictionBatches(params.openfold3_batch_size, 5, params.n_predictions)
    def base_seed = params.openfold3_seeds ? (params.openfold3_seeds.toString().split(',')[0].trim() as int) : 42

    GENERATE_OPENFOLD3_INPUT(ch_for_openfold3.filter { meta, _fasta, _msa -> (meta.n_chains ?: 1) == 1 })
    GENERATE_OPENFOLD3_INPUT_COMPLEX(ch_for_openfold3.filter { meta, _fasta, _msa -> (meta.n_chains ?: 1) > 1 })
    ch_with_json = GENERATE_OPENFOLD3_INPUT.out.with_json.mix(GENERATE_OPENFOLD3_INPUT_COMPLEX.out.with_json)

    ch_batched = ch_with_json.flatMap { meta, fasta, msa_dirs, query_json ->
        def is_mono = (meta.n_chains ?: 1) == 1
        def main_a3m = (msa_dirs instanceof List ? msa_dirs[0] : msa_dirs).resolve('colabfold_main.a3m')
        def n_seq = (is_mono && MsaSubsample.isEnabled(params.msa_subsample)) \
            ? MsaSubsample.countA3mSequences(main_a3m) : null
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
                    openfold3_seed: base_seed + i,
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
                else if (MsaSubsample.isEnabled(params.msa_subsample)) {
                    m = m + [msa_depth_tag: 'full']
                }
                jobs << [m, fasta, msa_dirs, query_json]
            }
        }
        jobs
    }

    OPENFOLD3_PROCESS(ch_batched)

    // One score row per OpenFold3 sample (seed_S_sample_N). ipSAE and PAE need
    // the sibling *_confidences.json; see protenix_fold.nf for why missing files
    // fall back to distinct dummy names.
    def no_ipsae_pae = file("${projectDir}/assets/dummy_files/no_ipsae_pae.json")
    def no_ipsae_struct = file("${projectDir}/assets/dummy_files/no_ipsae_structure.cif")

    ch_conf = OPENFOLD3_PROCESS.out.predictions.flatMap { meta, files ->
        def all = files instanceof List ? files : [files]
        def byName = [:]
        all.each { f -> byName[f.name] = f }
        all.findAll { it.name ==~ /.*_seed_\d+_sample_\d+_confidences_aggregated\.json$/ }
            .collect { j ->
                def mm = (j.name =~ /_(seed_\d+_sample_\d+)_confidences_aggregated\.json$/)
                def model = mm ? mm[0][1] : 'sample'
                def struct = j.name.replaceFirst(/_confidences_aggregated\.json$/, '_model.cif')
                def paeName = j.name.replaceFirst(/_confidences_aggregated\.json$/, '_confidences.json')
                def pae = byName[paeName]
                def cif = byName[struct]
                def doIpsae = (pae != null && cif != null)
                def pred = "${FoldNaming.flatPrefix('openfold3', meta)}${struct}"
                [meta, 'openfold3', model, struct, pred, j, doIpsae ? pae : no_ipsae_pae, doIpsae ? cif : no_ipsae_struct, doIpsae]
            }
    }
    FOLD_PARSE_CONFIDENCE(ch_conf)
    ch_tsv = FOLD_PARSE_CONFIDENCE.out.collectFile(
        name: 'openfold3_fold_scores.tsv',
        storeDir: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/openfold3",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = OPENFOLD3_PROCESS.out.predictions
    confidence_json = OPENFOLD3_PROCESS.out.confidence_json
    tsv = ch_tsv
}
