/*
ALPHAFOLD3_FOLD: AlphaFold3 folding for fold.nf / fold_pulldown.nf (--methods af3).

Consumes the same MSA bundle as PROTENIX_FOLD - monomer (meta, fasta, a3m);
multimer (meta, fasta, [paired...+unpaired...]) - and by default runs AF3
inference only (--run_data_pipeline=false). AF3 requires modelSeeds in the input
JSON and has no CLI seed flag, so batches are fanned out BEFORE JSON generation
and each gets seed base+i (base = --af3_seeds, else a fixed 1 so -resume hashes
stay stable).

With a comma list in --af3_seeds every job instead carries that whole seed list,
which is how Germinal runs AF3; see lib/AF3Input.groovy. --af3_run_data_pipeline
turns AF3's data pipeline on so it searches for templates (the Germinal
behaviour), which additionally requires --af3_db_dir.
*/

include { GENERATE_AF3_INPUT } from '../../modules/local/fold/af3/generate_af3_input'
include { GENERATE_AF3_INPUT_COMPLEX } from '../../modules/local/fold/af3/generate_af3_input_complex'
include { ALPHAFOLD3 as ALPHAFOLD3_PROCESS } from '../../modules/local/fold/af3/alphafold3'
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

workflow ALPHAFOLD3_FOLD {
    take:
    ch_for_af3 // monomer: tuple(meta, fasta, a3m); multimer: tuple(meta, fasta, [paired...+unpaired...])
    ch_templates // value: FOLD_TEMPLATES directory (or placeholder)

    main:
    // A comma list in --af3_seeds means "all of these seeds, in one job" (Germinal
    // style); a single value keeps the one-seed-per-batch fan-out. The seed list
    // IS the replication in that mode, so the --n_predictions fan-out is
    // suppressed - otherwise every batch would repeat the same seeds and the run
    // would produce len(seeds) x n_predictions structures instead of
    // len(seeds) x --af3_batch_size.
    def literal_seeds = AF3Input.literalSeeds(params)
    def batches = literal_seeds \
        ? [(params.af3_batch_size ?: 5) as int] \
        : foldPredictionBatches(params.af3_batch_size, 5, params.n_predictions)

    ch_batched = ch_for_af3.flatMap { meta, fasta, msa ->
        def is_mono = (meta.n_chains ?: 1) == 1
        def n_seq = (is_mono && MsaSubsample.isEnabled(params.msa_subsample)) \
            ? MsaSubsample.countA3mSequences(msa) : null
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
                    af3_seeds: AF3Input.seedsForBatch(params, i),
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
                jobs << [m, fasta, msa]
            }
        }
        jobs
    }

    GENERATE_AF3_INPUT(ch_batched.filter { meta, _fasta, _msa -> (meta.n_chains ?: 1) == 1 }, ch_templates)
    GENERATE_AF3_INPUT_COMPLEX(ch_batched.filter { meta, _fasta, _msa -> (meta.n_chains ?: 1) > 1 }, ch_templates)
    ch_with_json = GENERATE_AF3_INPUT.out.with_json.mix(GENERATE_AF3_INPUT_COMPLEX.out.with_json)

    // The databases are only read when the data pipeline runs; otherwise stage a
    // placeholder so the process signature stays the same.
    def af3_dbs = AF3Input.runDataPipeline(params) \
        ? file(params.af3_db_dir, checkIfExists: true) \
        : file("${projectDir}/assets/dummy_files", checkIfExists: true)

    ALPHAFOLD3_PROCESS(ch_with_json, file(params.af3_model_dir, checkIfExists: true), af3_dbs)

    // One score row per AF3 sample (seed-S_sample-N). ipSAE and pLDDT both need
    // the sibling *_confidences.json; see protenix_fold.nf for why missing files
    // fall back to distinct dummy names.
    def no_ipsae_pae = file("${projectDir}/assets/dummy_files/no_ipsae_pae.json")
    def no_ipsae_struct = file("${projectDir}/assets/dummy_files/no_ipsae_structure.cif")

    ch_conf = ALPHAFOLD3_PROCESS.out.predictions.flatMap { meta, files ->
        def all = files instanceof List ? files : [files]
        def byName = [:]
        all.each { f -> byName[f.name] = f }
        all.findAll { it.name ==~ /.*_seed-\d+_sample-\d+_summary_confidences\.json$/ }
            .collect { j ->
                def mm = (j.name =~ /_(seed-\d+_sample-\d+)_summary_confidences\.json$/)
                def model = mm ? mm[0][1] : 'sample'
                def struct = j.name.replaceFirst(/_summary_confidences\.json$/, '_model.cif')
                def paeName = j.name.replaceFirst(/_summary_confidences\.json$/, '_confidences.json')
                def pae = byName[paeName]
                def cif = byName[struct]
                def doIpsae = (pae != null && cif != null)
                def pred = "${FoldNaming.flatPrefix('af3', meta)}${struct}"
                [meta, 'af3', model, struct, pred, j, doIpsae ? pae : no_ipsae_pae, doIpsae ? cif : no_ipsae_struct, doIpsae]
            }
    }
    FOLD_PARSE_CONFIDENCE(ch_conf)
    ch_tsv = FOLD_PARSE_CONFIDENCE.out.collectFile(
        name: 'af3_fold_scores.tsv',
        storeDir: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/af3",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = ALPHAFOLD3_PROCESS.out.predictions
    confidence_json = ALPHAFOLD3_PROCESS.out.confidence_json
    tsv = ch_tsv
}
