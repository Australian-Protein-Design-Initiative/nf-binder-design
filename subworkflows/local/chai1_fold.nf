/*
CHAI1_FOLD: Chai-1 folding for fold.nf / fold_pulldown.nf (--methods chai1).

Consumes the same MSA bundle as PROTENIX_FOLD / OPENFOLD3_FOLD - monomer
(meta, fasta, a3m); multimer (meta, fasta, [paired...+unpaired...]). Chai-1 pairs
chains itself from per-row taxids, so the bundle needs no per-engine render and
goes straight to the CHAI1 process. With --chai1_seeds unset no seed is passed:
run_chai1.py draws one per task, so batches differ while -resume hashes stay
stable; with a base seed, batch i gets base+i.
*/

include { CHAI1 as CHAI1_PROCESS } from '../../modules/local/fold/chai1/chai1'
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

workflow CHAI1_FOLD {
    take:
    ch_for_chai1 // monomer: tuple(meta, fasta, a3m); multimer: tuple(meta, fasta, [paired...+unpaired...])
    ch_templates // value: FOLD_TEMPLATES directory (or placeholder)

    main:
    def batches = foldPredictionBatches(params.chai1_batch_size, 5, params.n_predictions)
    def base_seed = params.chai1_seeds ? (params.chai1_seeds.toString().split(',')[0].trim() as int) : null

    ch_batched = ch_for_chai1.flatMap { meta, fasta, a3ms ->
        def is_mono = (meta.n_chains ?: 1) == 1
        def main_a3m = a3ms instanceof List ? a3ms[0] : a3ms
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
                    chai1_seed: base_seed != null ? base_seed + i : null,
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
                jobs << [m, fasta, a3ms]
            }
        }
        jobs
    }

    CHAI1_PROCESS(ch_batched, ch_templates)

    // One score row per Chai-1 sample (seed_S_sample_N). run_chai1.py writes
    // AF3-shaped confidence files, as run_esmfold2.py does; see protenix_fold.nf
    // for why missing files fall back to distinct dummy names.
    def no_ipsae_pae = file("${projectDir}/assets/dummy_files/no_ipsae_pae.json")
    def no_ipsae_struct = file("${projectDir}/assets/dummy_files/no_ipsae_structure.cif")

    ch_conf = CHAI1_PROCESS.out.predictions.flatMap { meta, files ->
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
                def pred = "${FoldNaming.flatPrefix('chai1', meta)}${struct}"
                [meta, 'chai1', model, struct, pred, j, doIpsae ? pae : no_ipsae_pae, doIpsae ? cif : no_ipsae_struct, doIpsae]
            }
    }
    FOLD_PARSE_CONFIDENCE(ch_conf)
    ch_tsv = FOLD_PARSE_CONFIDENCE.out.collectFile(
        name: 'chai1_fold_scores.tsv',
        storeDir: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/chai1",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = CHAI1_PROCESS.out.predictions
    tsv = ch_tsv
}
