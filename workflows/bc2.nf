#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
BindCraft2 (BC2) binder design workflow

Usage via main.nf:
  nextflow run main.nf --method bc2 --bc2_config campaign.json
  nextflow run main.nf --method bc2 --input_pdb target.pdb --target_chains A

*/

params.outdir = 'results'

// Config mode: a BindCraft2 campaign settings JSON.
params.bc2_config = false

// Params mode: generate a campaign JSON from flat params (names shared with
// --method bindcraft so the two are interchangeable at the call site).
params.input_pdb = false
params.target_chains = false
params.hotspot_res = false
params.binder_length_range = '60-150'
params.design_name = 'bc2_design'

// Fraction of each target's hotspots to keep per batch. Applies in config mode
// as well, since it is done when the campaign JSON is staged into the task.
params.hotspot_subsample = 1.0

// Batching. BC2 stops at whichever comes first of number_of_final_designs
// accepted designs or max_trajectories attempts; this workflow always drives
// batching by trajectory budget, for predictable runtimes. See
// plans/i-d-like-to-add-sparkling-whale.md for why there is no --bc2_n_designs.
params.bc2_n_traj = 10
params.bc2_batch_size = 1
params.bc2_seed = 0

// Passthrough to `bindcraft design`.
params.bc2_core = false
params.bc2_modality = false
params.bc2_set = false
params.bc2_metadata = false
params.bc2_extra_args = ''

// Optional: only needed for a container that does not bake the AlphaFold
// parameters. The v1.0.1 image does.
params.bc2_af2_params = false

// Leave BC2's own multi-GPU fan-out on. Off by default so Nextflow controls
// parallelism and the gpu_lock.sh allocator is not second-guessed.
params.bc2_auto_multi_gpu = false

// Design property flags, one per settings/property/*.json in BindCraft2.
params.bc2_bigbang = false
params.bc2_disulfide_staple = false
params.bc2_forced_targeting = false
params.bc2_humanize = false
params.bc2_initial_guess = false
params.bc2_mixed_topology = false
params.bc2_protease_stable = false
params.bc2_termini_accessible = false
params.bc2_termini_together = false

// FoldSeek structural novelty search over the ranked designs (optional).
params.do_foldseek_search = false
params.foldseek_database = 'CATH50'
params.foldseek_databases_path = false
params.foldseek_use_webserver = false
params.foldseek_mode = '3diaa'
params.foldseek_maxaccept = 1
params.foldseek_gzip_output = false
params.foldseek_include_html_output = false
params.foldseek_cath_names_path = false

params.require_gpu = true
params.gpu_devices = ''
params.gpu_slots_per_device = 1
params.gpu_lock_timeout = 14400

include { BC2_GENERATE_CONFIG } from '../modules/local/bindcraft2/bc2_generate_config'
include { BINDCRAFT2 } from '../modules/local/bindcraft2/bindcraft2'
include { BC2_CAMPAIGN_OUTPUT } from '../modules/local/bindcraft2/bc2_campaign_output'
include { BC2_REPORTING } from '../modules/local/bindcraft2/bc2_reporting'
include {
    bc2AdaptiveNote;
    bc2BatchOverrides;
    bc2BinderChainCount;
    bc2BatchSizes;
    bc2CampaignSeed;
    bc2SoftDefaults;
    bc2StageArgs;
    buildBc2Args;
    extractBc2InputPaths;
    makeBc2BatchId;
    resolveBc2SetList;
    validateBc2SetList;
} from '../modules/local/bindcraft2/bindcraft2_utils'
include { resolveInputPdbs } from '../modules/local/bindcraft/bindcraft_utils'
include { FOLDSEEK_SEARCH } from '../subworkflows/local/foldseek_search'
include { FOLDSEEK_PREPARE_QUERIES } from '../modules/local/foldseek/foldseek_prepare_queries'

/*
 * BC2 hotspots are laxer than BindCraft v1's: chain prefixes are optional (an
 * unprefixed residue means the first selected chain) and ranges are allowed, eg
 * "A54,A56,B12-16". So this validates the shape without requiring a chain ID on
 * every residue, unlike parseHotspotResidues() in workflows/bindcraft.nf.
 */
def parseBc2HotspotResidues(hotspot_res) {
    // Omitting --hotspot_res leaves the default false; --hotspot_res "" parses as true.
    if (hotspot_res == null || hotspot_res == false || hotspot_res == true) {
        return false
    }
    if (!(hotspot_res instanceof CharSequence)) {
        error 'Invalid hotspot_res type: expected a string, or omitted/empty for no hotspots.'
    }

    def residues = hotspot_res.split(',').collect { residue -> residue.trim() }.findAll { residue -> !residue.isEmpty() }
    if (residues.isEmpty()) {
        return false
    }

    // {chain}{resnum}, {resnum}, {chain}{start}-{end} or {start}-{end}
    def pattern = ~/^[A-Za-z]?\d+(-\d+)?$/
    residues.each { residue ->
        if (!(residue =~ pattern)) {
            error "Invalid hotspot residue '${residue}'. BindCraft2 hotspots look like 'A54,A56' or 'B12-16' (chain prefix optional, ranges allowed)."
        }
    }

    def normalised = residues.join(',')
    if (normalised != hotspot_res) {
        log.warn("Normalised --hotspot_res from '${hotspot_res}' to '${normalised}'")
    }
    return normalised
}

def validateBc2BinderLengthRange(binder_length_range) {
    if (binder_length_range == null || binder_length_range == false) {
        return true
    }
    def text = binder_length_range.toString().trim()
    if (text.isEmpty()) {
        error 'Invalid binder_length_range: cannot be empty.'
    }
    // '80' exact, '60-150' inclusive range, '60,80,100' discrete choices
    if (!(text ==~ /^\d+(-\d+)?$/) && !(text ==~ /^\d+(,\d+)+$/)) {
        error "Invalid binder_length_range '${binder_length_range}'. Use '60-150' for a range, '80' for an exact length, or '60,80,100' for discrete choices."
    }
    if (text.contains('-')) {
        def parts = text.split('-')
        if (Integer.parseInt(parts[0]) > Integer.parseInt(parts[1])) {
            error "Invalid binder_length_range '${binder_length_range}': the first value must not exceed the second."
        }
    }
    return true
}

workflow BC2 {

    main:

    if (!params.bc2_config && !params.input_pdb) {
        log.info(
            """
        ==================================================================
        BINDCRAFT2 NEXTFLOW WRAPPER
        ==================================================================

        Runs BindCraft2 (https://github.com/PacesaLab/BindCraft2) campaigns in
        parallel batches, then gathers them with `bindcraft campaign_output`.

        One of these is required:
            --bc2_config              BindCraft2 campaign settings JSON.
                                      Files it references (targets[].target_path,
                                      binder_scaffold) are staged automatically.
            --input_pdb               Target PDB/mmCIF/FASTA. A campaign JSON is
                                      generated from the params below.
                                      A directory or a quoted glob runs each structure.

        Params mode options:
            --target_chains           Target chain(s), eg "A" or "A,B"
            --hotspot_res             Hotspots, eg "A54,A56" or "B12-16" [default: none]
            --binder_length_range     "60-150" range, "80" exact, or "60,80,100" choices [default: ${params.binder_length_range}]
            --design_name             campaign_name for the generated config [default: ${params.design_name}]

        Both modes:
            --hotspot_subsample       Fraction of each target's hotspots to keep per batch,
                                      ceil(N*fraction) with a minimum of 1. Each batch draws its
                                      own subset from its campaign_seed [default: ${params.hotspot_subsample}]

        Batching:
            --bc2_n_traj              Total trajectories to run [default: ${params.bc2_n_traj}]
            --bc2_batch_size          Trajectories per parallel task [default: ${params.bc2_batch_size}]
            --bc2_seed                Base campaign_seed; batch i uses bc2_seed + i [default: ${params.bc2_seed}]

        BindCraft2 settings:
            --bc2_set                 Any BindCraft2 setting as KEY=VALUE, ';'-separated,
                                      eg --bc2_set 'binder_lengths=[70,90];filters.i_pTM.threshold=0.8'
                                      Dotted keys nest. Values are JSON, falling back to a string.
                                      See settings/core/reference.json upstream for all 235 settings.
            --bc2_modality            Binder format, eg "VHH" or "VHH,induced_fit"
                                      (binder, large_binder, peptide, cyclic_peptide, homo_oligomer,
                                       multidomain, VHH, scFv, Fab, ARP, induced_fit, fold_switch)
            --bc2_core                Core profile, eg "benchmark" (pins the seed, disables
                                      autotune and desperation for reproducible runs)
            --bc2_metadata            JSON file of fields to record with the campaign
            --bc2_extra_args          Extra arguments appended to `bindcraft design`

        Design properties (boolean flags):
            --bc2_humanize            Human-like sequence preferences, MHC-II anchor objective
            --bc2_bigbang             Start gradient stages from coordinates on hand
            --bc2_disulfide_staple    Allow cysteine, require a detected disulfide
            --bc2_forced_targeting    Concentrate the paratope on declared hotspots
            --bc2_initial_guess       Re-predict the candidate from the trajectory pose
            --bc2_mixed_topology      At most 50% helix, at least 20% beta sheet
            --bc2_protease_stable     Penalise cleavage motifs and exposed loops/termini
            --bc2_termini_accessible  Point both chain ends away from the target
            --bc2_termini_together    Keep the N and C termini within 10 A

        GPU:
            --require_gpu           Fail fast if no GPU is detected [default: ${params.require_gpu}]
            --gpu_devices           GPU devices to use (comma-separated or 'all') [default: ${params.gpu_devices}]
            --gpu_slots_per_device  Concurrent tasks per GPU [default: ${params.gpu_slots_per_device}]
            --gpu_lock_timeout      Seconds a task waits for a free GPU [default: ${params.gpu_lock_timeout}]
            --bc2_auto_multi_gpu    Let BindCraft2 fan out across GPUs itself [default: false]
            --bc2_af2_params        AlphaFold parameter directory (only for a no-weights image)

        FoldSeek (optional, enabled with --do_foldseek_search):
            --do_foldseek_search          Enable FoldSeek structural similarity search on ranked designs
            --foldseek_database           Database to search [default: ${params.foldseek_database}]
                                           Local search (default):
                                             - CATH50 (default) - CATH domain DB, combined AF2+PDB at 50% seq.id.
                                             - PDB              - Protein Data Bank
                                             - Alphafold/UniProt50 - AF2 clustered at 50% seq.id.
                                             - ESMAtlas30       - ESM metagenomic atlas at 30% seq.id.
                                             - BFMD             - Big Fantastic Multimer Database
                                             - BFVD             - Big Fantastic Virus Database
                                           Remote search (--foldseek_use_webserver true):
                                             - pdb100, afdb50, afdb-swissprot, afdb-proteome,
                                               cath50, bfmd, BFVD, mgnify_esm30
            --foldseek_databases_path     Path to local databases directory (auto-downloaded if unset)
            --foldseek_use_webserver      Use FoldSeek web API instead of local search [default: ${params.foldseek_use_webserver}]
            --foldseek_mode               Search mode [default: ${params.foldseek_mode}]
            --foldseek_maxaccept          Max accepted alignments per query [default: ${params.foldseek_maxaccept}]
            --foldseek_gzip_output        Gzip TSV output [default: ${params.foldseek_gzip_output}]
            --foldseek_include_html_output  Include HTML report [default: ${params.foldseek_include_html_output}]
            --foldseek_cath_names_path    Path to local cath-names.txt (auto-downloaded if unset)
                                          CATH annotation is automatic when the database name starts with "CATH"

        Example:
            nextflow run main.nf --method bc2 \\
                --bc2_config campaign.json \\
                --bc2_n_traj 100 --bc2_batch_size 5 \\
                --bc2_modality VHH --bc2_humanize \\
                -profile local

        """.stripIndent()
        )
        exit(1)
    }

    if (params.bc2_config && params.input_pdb) {
        error '--bc2_config and --input_pdb are mutually exclusive. Use --bc2_config (the target is named in the config) or --input_pdb with --target_chains.'
    }

    def bc2_set_list = resolveBc2SetList(params.bc2_set)
    validateBc2SetList(bc2_set_list)

    // Not `?:` -- Groovy treats 0 as falsy, so `--bc2_n_traj 0` would silently
    // become the fallback instead of failing validation below.
    def n_traj = (params.bc2_n_traj == null || params.bc2_n_traj == false) ? 0 : (params.bc2_n_traj as int)
    def batch_size = (params.bc2_batch_size == null || params.bc2_batch_size == false) ? 0 : (params.bc2_batch_size as int)
    if (n_traj < 1) {
        error "--bc2_n_traj must be at least 1 (got ${params.bc2_n_traj})"
    }
    if (batch_size < 1) {
        error "--bc2_batch_size must be at least 1 (got ${params.bc2_batch_size})"
    }

    def hotspot_subsample = (params.hotspot_subsample == null || params.hotspot_subsample == false)
        ? 1.0
        : (params.hotspot_subsample as Double)
    if (hotspot_subsample <= 0.0 || hotspot_subsample > 1.0) {
        error "--hotspot_subsample must be greater than 0.0 and at most 1.0 (got ${params.hotspot_subsample})"
    }

    def batch_sizes = bc2BatchSizes(n_traj, batch_size)
    def soft_defaults = bc2SoftDefaults(params)

    // Metadata is optional, but a Nextflow path input cannot be conditionally
    // absent; a sentinel file keeps the process signature fixed and is detected
    // by name in the module.
    def ch_metadata = params.bc2_metadata
        ? channel.fromPath(params.bc2_metadata, checkIfExists: true).first()
        : channel.fromPath("${projectDir}/assets/dummy_files/NO_BC2_METADATA").first()

    def ch_jobs
    def ch_config

    if (params.bc2_config) {
        def config_file = file(params.bc2_config, checkIfExists: true)
        // targets[].target_path and binder_scaffold are the only file-path
        // settings, and resolve against the config's own directory. An empty
        // list means the campaign names only shipped presets, whose paths
        // resolve inside the container.
        def input_paths = extractBc2InputPaths(config_file.toString(), config_file.parent.toString())
        input_paths.each { p ->
            if (!file(p).exists()) {
                error "File referenced by --bc2_config not found: ${p}"
            }
        }
        log.info("BindCraft2: campaign ${config_file.name}" +
                 (input_paths ? ", staging ${input_paths.size()} referenced file(s)" : ', no files to stage (shipped presets only)'))

        // Subsampling rewrites targets[].hotspots in the staged copy of the
        // campaign JSON, so hotspots that come from a shipped target preset
        // (eg {"target": "hPDL1"}) are out of reach -- BindCraft2 merges those
        // in itself, after staging.
        if (hotspot_subsample < 1.0) {
            log.warn('--hotspot_subsample only thins hotspots written in the campaign JSON\'s own targets[], not those inherited from a shipped "target" preset.')
        }

        ch_config = channel.fromPath(config_file).first()
        def ch_input_files = channel.of(input_paths.collect { path -> file(path) }).first()

        ch_jobs = channel
            .fromList(batch_sizes.withIndex().collect { size, index -> [index, size] })
            .map { index, size ->
                def batch_id = "batch_${index}"
                def args = buildBc2Args(params, soft_defaults, bc2BatchOverrides(params, batch_id, index, size))
                tuple(batch_id, args, bc2StageArgs(hotspot_subsample, bc2CampaignSeed(params, index)))
            }
            .combine(ch_config)
            .combine(ch_input_files)
    }
    else {
        def hotspot_res = parseBc2HotspotResidues(params.hotspot_res)
        validateBc2BinderLengthRange(params.binder_length_range)

        if (hotspot_subsample < 1.0 && !hotspot_res) {
            log.warn('--hotspot_subsample has no effect without --hotspot_res')
        }

        def input_pdbs = resolveInputPdbs(params.input_pdb)
        log.info("BindCraft2: ${input_pdbs.size()} input structure(s) x ${n_traj} trajectories " +
                 "= ${input_pdbs.size() * n_traj} total, in ${input_pdbs.size() * batch_sizes.size()} batch task(s)")

        BC2_GENERATE_CONFIG(
            channel.fromList(input_pdbs),
            params.design_name,
            params.target_chains ?: '',
            hotspot_res ?: '',
            params.binder_length_range,
        )

        // Pair each generated config back with its structure so the batch id
        // carries the structure name, and stage the structure alongside.
        def ch_generated = channel.fromList(input_pdbs).merge(BC2_GENERATE_CONFIG.out.config_json)

        ch_jobs = ch_generated
            .combine(channel.fromList(batch_sizes.withIndex().collect { size, index -> [index, size] }))
            .map { pdb, config, index, size ->
                def batch_id = makeBc2BatchId(pdb, index)
                // seedIndex must be unique across the whole run, not per
                // structure, so batches of different targets do not share an
                // RNG stream.
                def seed_index = (input_pdbs.indexOf(pdb) * batch_sizes.size()) + index
                def args = buildBc2Args(params, soft_defaults, bc2BatchOverrides(params, batch_id, seed_index, size))
                def stage_args = bc2StageArgs(hotspot_subsample, bc2CampaignSeed(params, seed_index))
                tuple(batch_id, args, stage_args, config, [pdb])
            }
    }

    def adaptive_note = bc2AdaptiveNote(params, batch_size)
    if (adaptive_note) {
        log.info(adaptive_note)
    }

    BINDCRAFT2(ch_jobs, ch_metadata)

    // Stage every campaign under one parent so `bindcraft campaign_output`
    // discovers them all and writes a combined summary.
    BC2_CAMPAIGN_OUTPUT(BINDCRAFT2.out.batch_dir.map { _id, dir -> dir }.collect())

    BC2_REPORTING(
        BC2_CAMPAIGN_OUTPUT.out.summary,
        BC2_CAMPAIGN_OUTPUT.out.ranked.ifEmpty(file("${projectDir}/assets/dummy_files/NO_BC2_RANKED.csv")),
        BC2_CAMPAIGN_OUTPUT.out.refolded.ifEmpty(file("${projectDir}/assets/dummy_files/NO_BC2_REFOLDED.csv")),
        BC2_CAMPAIGN_OUTPUT.out.trajectories.ifEmpty(file("${projectDir}/assets/dummy_files/NO_BC2_TRAJECTORIES.csv")),
    )

    if (params.do_foldseek_search) {
        // ranked_cifs rather than the batch directories: it is populated only
        // after `bindcraft campaign_output` has re-ranked every campaign, which
        // is the point at which 3_Ranked/ is authoritative. It is also
        // `optional: true`, so a run that accepts nothing skips the search
        // instead of failing.
        //
        // BindCraft2 writes every target chain before the binder
        // (written_chains() in bindcraft/protein.py sorts on is_binder_chain),
        // so the binder is the last N chains. N comes from the campaign
        // metadata rather than the usual smallest-chain guess, which picks the
        // wrong chain for a multi-chain binder (Fab, scFv, multidomain, an
        // oligomer) or a target smaller than its binder.
        ch_n_binder_chains = BINDCRAFT2.out.campaign_metadata
            .first()
            .map { metadata -> bc2BinderChainCount(metadata) }
            .ifEmpty(0)

        FOLDSEEK_PREPARE_QUERIES(
            BC2_CAMPAIGN_OUTPUT.out.ranked_cifs.flatten().collect(),
            ch_n_binder_chains,
        )

        FOLDSEEK_SEARCH(
            FOLDSEEK_PREPARE_QUERIES.out.design_chains,
            channel.of([id: 'bc2_foldseek']),
        )
    }

    emit:
    summary = BC2_CAMPAIGN_OUTPUT.out.summary
    ranked = BC2_CAMPAIGN_OUTPUT.out.ranked
    report = BC2_REPORTING.out
    foldseek_tsv = params.do_foldseek_search ? FOLDSEEK_SEARCH.out.tsv : channel.empty()
    foldseek_tsv_annotated = params.do_foldseek_search ? FOLDSEEK_SEARCH.out.tsv_annotated : channel.empty()
}
