/*
BOLTZ_FOLD: generic single-target Boltz-2 folding for fold.nf.

Monomer inputs (meta.n_chains == 1) use FOLD_CREATE_BOLTZ_YAML (one protein
entry, shared a3m); multimer inputs (n_chains > 1) use
FOLD_CREATE_BOLTZ_YAML_COMPLEX (one protein entry per chain, each with its own
key,sequence CSV from FOLD_MSA; see plans/fold-nf-multimer-paired-msa.md). Both
feed the same BOLTZ predict.
*/

include { FOLD_CREATE_BOLTZ_YAML } from '../../modules/local/fold/boltz/fold_create_boltz_yaml'
include { FOLD_CREATE_BOLTZ_YAML_COMPLEX } from '../../modules/local/fold/boltz/fold_create_boltz_yaml_complex'
include { BOLTZ } from '../../modules/local/common/boltz'
include { FOLD_PARSE_BOLTZ_CONFIDENCE } from '../../modules/local/fold/boltz/fold_parse_boltz_confidence'
include { foldNamingFlatPrefix } from '../../modules/local/fold/common/fold_naming'
include { foldParamsEngine } from '../../modules/local/fold/common/fold_params'
include { msaSubsampleCountA3mSequences; msaSubsampleDepthJobs; msaSubsampleDepthTag; msaSubsampleIsEnabled; msaSubsampleStableSeed } from '../../modules/local/fold/common/msa_subsample'

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
    ch_for_boltz // monomer: tuple(meta, fasta, a3m); multimer: tuple(meta, fasta, [csv...])
    ch_templates // value: FOLD_TEMPLATES directory (or placeholder)
    // publish_label keeps two fold runs that share an output directory apart: the refold
    // folds the complexes and then the binder monomers, and both write their per-engine
    // score table to the same storeDir. Without a label the second overwrites the first,
    // and the merge downstream then sees rows for complexes it never folded.
    publish_label

    main:
    ch_mono = ch_for_boltz.filter { meta, fasta, msa -> (meta.n_chains ?: 1) == 1 }
    ch_multi = ch_for_boltz.filter { meta, fasta, msa -> (meta.n_chains ?: 1) > 1 }

    FOLD_CREATE_BOLTZ_YAML(ch_mono, ch_templates)
    FOLD_CREATE_BOLTZ_YAML_COMPLEX(ch_multi, ch_templates)
    ch_yaml = FOLD_CREATE_BOLTZ_YAML.out.yaml.mix(FOLD_CREATE_BOLTZ_YAML_COMPLEX.out.yaml)

    // BOLTZ's process signature dates from a target+binder complex mode, so
    // it always expects two MSA paths to stage (their names never appear in
    // the script body - they only need to be physically present alongside the
    // YAML, since the YAML's `msa:` field is what boltz predict actually
    // reads). fold.nf has one MSA, so the unused "target" slot gets an
    // empty-placeholder file.
    // default_batch = 5: with --n_predictions unset, emit 5 samples like RF3 /
    // Protenix (Boltz's own native default is 1; we lift it for cross-engine parity).
    // Via foldParamsEngine so a refolding workflow can namespace these away from its own
    // same-named design-time parameters (see modules/local/fold/common/fold_params.nf).
    def batches = foldPredictionBatches(foldParamsEngine(params, 'boltz_batch_size'), 5, params.n_predictions)
    def boltz_seed = foldParamsEngine(params, 'boltz_seed')
    def base_seed = boltz_seed ? (boltz_seed as int) : null

    ch_boltz_input = ch_yaml.flatMap { meta, yaml, msa ->
        def n_seq = msaSubsampleIsEnabled(params.msa_subsample) \
            ? msaSubsampleCountA3mSequences(msa) : null
        def depth_jobs = msaSubsampleDepthJobs(
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
                    m = m + [boltz_seed: base_seed + i]
                }
                if (depth != null) {
                    def s = msaSubsampleStableSeed(meta.id.toString(), i + 1, depth[0], depth[1])
                    m = m + [
                        msa_max_seq: depth[0],
                        msa_max_extra_seq: depth[1],
                        msa_subsample_seed: s,
                        msa_depth_tag: msaSubsampleDepthTag(depth[0], depth[1]),
                    ]
                }
                else if (msaSubsampleIsEnabled(params.msa_subsample)) {
                    m = m + [msa_depth_tag: 'full']
                }
                jobs << [m, yaml, file("${projectDir}/assets/dummy_files/empty_target_msa"), msa]
            }
        }
        jobs
    }

    // step_name is BOLTZ's publishDir subdir under params.outdir. Honour
    // fold_publish_dir so fold_pulldown lands at <outdir>/fold_pulldown/boltz/
    // (fold.nf keeps <outdir>/fold/boltz/). mmCIF so predictions stay
    // format-uniform with af2/rf3/protenix.
    BOLTZ(ch_boltz_input, ch_templates, "${params.fold_publish_dir ?: 'fold'}/boltz", 'cif')

    // Fan out to one row per diffusion sample so all models are scored (not
    // just the top-ranked model_0). The all-samples emit is a list when >1
    // sample, a single path when 1; normalise, then tag each with its model
    // index parsed from the confidence_<id>_model_<N>.json filename.
    empty_ipsae_tsv = file("${projectDir}/assets/dummy_files/empty_ipsae.tsv")
    // Zero bytes, and parse_boltz_confidence.py reads that as "no PAE matrix" rather
    // than failing. Only reachable if Boltz stopped writing pae_*.npz.
    empty_pae_npz = file("${projectDir}/assets/dummy_files/empty_pae.npz")

    // Pair each model's confidence JSON with its own ipSAE TSV (ipsae.py writes
    // one <model>_10_10_ipsae.tsv per diffusion sample - see boltz.nf's per-pae
    // loop) so B2's ipsae_d0chn/ipsae_d0dom/pdockq/pdockq2/lis columns get
    // filled for Boltz too, not just RF3/Protenix/AF3/OpenFold3.
    ch_conf_per_model = BOLTZ.out.confidence_json_all
        .join(BOLTZ.out.ipsae_tsv)
        .join(BOLTZ.out.pae_npz)
        .flatMap { meta, jsons, ipsae_tsvs, pae_npzs ->
            def files = (jsons instanceof List) ? jsons : [jsons]
            def ipsae_files = (ipsae_tsvs instanceof List) ? ipsae_tsvs : [ipsae_tsvs]
            def pae_files = (pae_npzs instanceof List) ? pae_npzs : [pae_npzs]
            files.collect { j ->
                def m = (j.name =~ /_model_(\d+)\.json$/)
                def idx = m ? m[0][1] : '0'
                // confidence_<id>_model_N.json -> native <id>_model_N.cif; the
                // flat-gather name uses the same foldNamingFlatPrefix as the module.
                def struct = j.name.replaceFirst(/^confidence_/, '').replaceFirst(/\.json$/, '.cif')
                def pred = "${foldNamingFlatPrefix('boltz', meta)}${struct}"
                def stem = struct.replaceFirst(/\.(cif|pdb)$/, '')
                def ipsae = ipsae_files.find { it.name.startsWith("${stem}_") && it.name.endsWith('_ipsae.tsv') }
                // Boltz names them pae_<stem>.npz, one per diffusion sample.
                def pae_npz = pae_files.find { it.name == "pae_${stem}.npz" }
                [meta, idx, struct, pred, j, ipsae ?: empty_ipsae_tsv, pae_npz ?: empty_pae_npz]
            }
        }

    FOLD_PARSE_BOLTZ_CONFIDENCE(ch_conf_per_model)

    ch_tsv = FOLD_PARSE_BOLTZ_CONFIDENCE.out.collectFile(
        name: "${publish_label ? publish_label + '_' : ''}boltz_fold_scores.tsv",
        storeDir: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/boltz",
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
