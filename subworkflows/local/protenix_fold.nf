/*
PROTENIX_FOLD: generic Protenix (AF3-style) folding for fold.nf, the 4th
--methods engine (see plans/fold-nf-multi-method-folding.md Phase 3).

Monomer inputs (meta.n_chains == 1) use GENERATE_PROTENIX_INPUT with a single
unpairedMsaPath a3m; multimer inputs (n_chains > 1) use
GENERATE_PROTENIX_INPUT_COMPLEX with per-chain paired + unpaired a3m bundles
from FOLD_MSA (Protenix pairs by species mnemonic; see
plans/fold-nf-multimer-paired-msa.md). Both feed the same PROTENIX_FOLD predict.
*/

include { GENERATE_PROTENIX_INPUT } from '../../modules/local/fold/protenix/generate_protenix_input'
include { GENERATE_PROTENIX_INPUT_COMPLEX } from '../../modules/local/fold/protenix/generate_protenix_input_complex'
include { PROTENIX_FOLD as PROTENIX_FOLD_PROCESS } from '../../modules/local/fold/protenix/protenix_fold'
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

workflow PROTENIX_FOLD {
    take:
    ch_for_protenix // monomer: tuple(meta, fasta, a3m); multimer: tuple(meta, fasta, [paired...+unpaired...])
    ch_templates // value: FOLD_TEMPLATES directory (or placeholder)
    // publish_label keeps two fold runs that share an output directory apart: the refold
    // folds the complexes and then the binder monomers, and both write their per-engine
    // score table to the same storeDir. Without a label the second overwrites the first,
    // and the merge downstream then sees rows for complexes it never folded.
    publish_label

    main:
    ch_mono = ch_for_protenix.filter { meta, fasta, msa -> (meta.n_chains ?: 1) == 1 }
    ch_multi = ch_for_protenix.filter { meta, fasta, msa -> (meta.n_chains ?: 1) > 1 }
    GENERATE_PROTENIX_INPUT(ch_mono, ch_templates)
    GENERATE_PROTENIX_INPUT_COMPLEX(ch_multi, ch_templates)
    ch_with_json = GENERATE_PROTENIX_INPUT.out.with_json.mix(GENERATE_PROTENIX_INPUT_COMPLEX.out.with_json)

    // Via FoldParams so a refolding workflow can namespace these away from its own
    // same-named design-time parameters (see lib/FoldParams.groovy).
    def batches = foldPredictionBatches(FoldParams.engine(params, 'protenix_batch_size'), 5, params.n_predictions)

    // --protenix_seeds may be a comma list; the first value is the base the per-batch
    // fan-out counts from, as for the other engines.
    def protenix_seeds = FoldParams.engine(params, 'protenix_seeds')
    def base_seed = protenix_seeds \
        ? (protenix_seeds.toString().split(',')[0].trim() as int) : null

    // runner/inference.py loops `for seed in seeds:` with the model resident, so the
    // one-job-per-seed fan-out collapses into a single job (lib/FoldSeeds.groovy).
    // This is only safe because the seed is put back into the published and scored
    // names above - Protenix itself separates seeds by directory alone.
    def seed_list = FoldSeeds.collapse('protenix', batches, base_seed, params.fold_seed_list)
    if (seed_list) {
        batches = [batches[0]]
    }

    ch_batched = ch_with_json.flatMap { meta, fasta, a3m, json, protenix_templates ->
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
                    m = m + [protenix_seeds: seed_list ? seed_list.join(',') : "${base_seed + i}"]
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
                jobs << [m, fasta, a3m, json, protenix_templates]
            }
        }
        jobs
    }

    PROTENIX_FOLD_PROCESS(ch_batched)

    // One score row per Protenix sample. model/struct/predictions_file are
    // derived from the summary filename (complex_summary_confidence_sample_N.json
    // -> complex_sample_N.cif); predictions_file uses the same FoldNaming prefix
    // as the module's flat-gather saveAs. ipSAE needs the sibling
    // *_full_data_sample_N.json (token-pair PAE) and the mmCIF.
    // When a sample's PAE JSON or mmCIF is missing (e.g. a prediction job was
    // interrupted mid-write) ipSAE is skipped -- but FOLD_PARSE_CONFIDENCE still
    // declares three separate path inputs, so they must be three DISTINCT file
    // names or Nextflow aborts the whole run with "input file name collision".
    // Passing the summary JSON in all three slots (the old behaviour) crashed the
    // pipeline on the first incomplete sample.
    def no_ipsae_pae = file("${projectDir}/assets/dummy_files/no_ipsae_pae.json")
    def no_ipsae_struct = file("${projectDir}/assets/dummy_files/no_ipsae_structure.cif")

    // Keyed by seed as well as basename: a multi-seed job writes one
    // <name>_sample_<r>.cif per seed, all with the same basename, under different
    // seed_<N>/ directories. Keying on the basename alone would make every seed's
    // row point at whichever file was indexed last.
    ch_conf = PROTENIX_FOLD_PROCESS.out.predictions.flatMap { meta, files ->
        def all = files instanceof List ? files : [files]
        def byName = [:]
        all.each { f -> byName["${FoldNaming.protenixSeedBit(f)}${f.name}".toString()] = f }
        all.findAll { it.name ==~ /.*_summary_confidence_sample_\d+\.json$/ }
            .collect { j ->
                def seed_bit = FoldNaming.protenixSeedBit(j)
                def mm = (j.name =~ /_summary_confidence_sample_(\d+)\.json$/)
                def idx = mm ? mm[0][1] : '0'
                def struct = j.name.replaceFirst(/_summary_confidence_sample_(\d+)\.json$/, '_sample_$1.cif')
                def paeName = j.name.replaceFirst(/_summary_confidence_sample_(\d+)\.json$/, '_full_data_sample_$1.json')
                def pae = byName["${seed_bit}${paeName}".toString()]
                def cif = byName["${seed_bit}${struct}".toString()]
                def doIpsae = (pae != null && cif != null)
                def pred = "${FoldNaming.flatPrefix('protenix', meta)}${seed_bit}${struct}"
                // The model label carries the seed so one design's samples stay distinct
                // in the score table, as they do for every other engine.
                def model = "${seed_bit}sample_${idx}".toString()
                [meta, 'protenix', model, struct, pred, j, doIpsae ? pae : no_ipsae_pae, doIpsae ? cif : no_ipsae_struct, doIpsae]
            }
    }
    FOLD_PARSE_CONFIDENCE(ch_conf)
    ch_tsv = FOLD_PARSE_CONFIDENCE.out.collectFile(
        name: "${publish_label ? publish_label + '_' : ''}protenix_fold_scores.tsv",
        storeDir: "${params.outdir}/${params.fold_publish_dir ?: 'fold'}/protenix",
        keepHeader: true,
        skip: 1,
    )

    emit:
    predictions = PROTENIX_FOLD_PROCESS.out.predictions
    confidence_json = PROTENIX_FOLD_PROCESS.out.confidence_json
    tsv = ch_tsv
}
