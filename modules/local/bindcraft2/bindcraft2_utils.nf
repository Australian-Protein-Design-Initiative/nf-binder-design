// Utility functions for the BindCraft2 (bc2) workflow.
//
// Kept free of `import` and `@Field` so `nextflow lint` stays quiet on this
// file; groovy.json.JsonSlurper is referenced by its full name instead.

// Names of settings/property/*.json in BindCraft2. Each has a CLI flag whose
// name is the setting with underscores as dashes (bindcraft/cli.py:60-61), and
// each is also a plain boolean setting, so `--humanize` == `--set humanize=true`.
def bc2PropertyNames() {
    return [
        'bigbang',
        'disulfide_staple',
        'forced_targeting',
        'humanize',
        'initial_guess',
        'mixed_topology',
        'protease_stable',
        'termini_accessible',
        'termini_together',
    ]
}

/**
 * True when a param was actually set by the user (not null, not Groovy false,
 * not the string "false"). Mirrors mpnnWorkflowParamSet() in rfd3_utils.nf.
 */
boolean bc2ParamSet(def v) {
    if (v == null || v == false) {
        return false
    }
    def s = v.toString().trim()
    return s && !s.equalsIgnoreCase('false')
}

/**
 * Normalise --bc2_set into a list of KEY=VALUE strings.
 *
 * Nextflow does not accumulate a repeated `--bc2_set` on the command line (the
 * last one wins), so a single semicolon-separated string is the CLI form, while
 * a real list works from -params-file or a config file. Semicolon rather than
 * comma because BC2 values are JSON and routinely contain commas, eg
 * `binder_lengths=[70,90]`.
 */
List<String> resolveBc2SetList(def v) {
    if (!bc2ParamSet(v)) {
        return []
    }
    def items = (v instanceof List || v.getClass().isArray())
        ? (v as List).collect { item -> item?.toString() }
        : v.toString().split(';') as List
    return items.collect { item -> item?.trim() }.findAll { item -> item }
}

// Settings the workflow computes per batch. A user --bc2_set for one of these
// would either be silently clobbered or silently break batching, so they are
// rejected with a pointer to the param that does control them.
def bc2ReservedSettings() {
    return [
        'project_folder'  : 'the batch id is used, so batches cannot collide in the gather step',
        'campaign_seed'   : 'set the base with --bc2_seed; each batch needs a distinct seed',
        'max_trajectories': 'set --bc2_n_traj and --bc2_batch_size instead',
    ]
}

/**
 * Validate KEY=VALUE assignments in the driver, so a typo or a reserved setting
 * fails before a GPU task starts rather than after. BC2 raises the malformed
 * case itself, but only once the campaign is launching.
 */
def validateBc2SetList(List<String> assignments) {
    def reserved = bc2ReservedSettings()
    assignments.each { a ->
        def idx = a.indexOf('=')
        if (idx < 1) {
            throw new Exception(
                "--bc2_set entries are written KEY=VALUE, as in 'binder_lengths=[70,90]' or " +
                "'filters.i_pTM.threshold=0.8'; got '${a}'. Separate multiple entries with ';'."
            )
        }
        def key = a.substring(0, idx).trim()
        if (reserved.containsKey(key)) {
            throw new Exception(
                "--bc2_set ${key}=... is not allowed: the workflow sets it per batch (${reserved[key]})."
            )
        }
    }
    return true
}

/**
 * Single-quote a value for the shell, escaping any embedded single quotes.
 * BC2 values such as [70,90] and {"a":1} would otherwise be mangled.
 */
String bc2ShellQuote(String s) {
    return "'" + s.replace("'", "'\\''") + "'"
}

/**
 * Build the full `bindcraft design` argument string (everything after the
 * campaign JSON).
 *
 * Order matters. BC2 concatenates preset assignments before --set assignments
 * and applies them in order, so later ones win (bindcraft/cli.py:251-256).
 * Three tiers are emitted:
 *
 *   1. soft defaults  - values the pipeline prefers but the user may override
 *   2. user --bc2_set - overrides tier 1
 *   3. batch overrides - pipeline-managed settings that must always win
 *
 * Tier 3 keys are rejected in --bc2_set by validateBc2SetList(), so tier 2 and
 * tier 3 never actually contend; the ordering just makes that explicit.
 */
String buildBc2Args(params, List<String> softDefaults, List<String> batchOverrides) {
    def args = []

    if (bc2ParamSet(params.bc2_core)) {
        args << "--core ${params.bc2_core}"
    }
    if (bc2ParamSet(params.bc2_modality)) {
        args << "--modality ${params.bc2_modality}"
    }

    bc2PropertyNames().each { name ->
        if (bc2ParamSet(params["bc2_${name}"])) {
            args << "--${name.replace('_', '-')}"
        }
    }

    (softDefaults ?: []).each { assignment ->
        args << "--set ${bc2ShellQuote(assignment)}"
    }
    resolveBc2SetList(params.bc2_set).each { assignment ->
        args << "--set ${bc2ShellQuote(assignment)}"
    }
    (batchOverrides ?: []).each { assignment ->
        args << "--set ${bc2ShellQuote(assignment)}"
    }

    return args.join(' ')
}

// A campaign stops at whichever comes first of number_of_final_designs accepted
// designs or max_trajectories attempts. This workflow always drives batching by
// trajectory budget, so the accepted-design target has to be taken out of play.
//
// It cannot simply be left unset: number_of_final_designs is absent from
// settings/core/default.json, so BC2 falls back to
// settings.get('number_of_final_designs', 1) == 1 (campaign.py:135) and a batch
// would stop at its first acceptance, throwing away the rest of its budget.
def bc2UnlimitedDesigns() { return 1000000 }

/**
 * Settings the pipeline prefers but a user --bc2_set may override.
 */
List<String> bc2SoftDefaults(params) {
    def defaults = []

    defaults << "number_of_final_designs=${bc2UnlimitedDesigns()}".toString()

    // BC2 fans workers out across every visible card and packs several per card
    // by default (design_workers.py:239), which would fight the pipeline's own
    // GPU allocator. One process per task unless explicitly asked otherwise.
    if (!bc2ParamSet(params.bc2_auto_multi_gpu)) {
        defaults << 'auto_multi_gpu=false'
    }

    return defaults
}

/**
 * Per-batch --set overrides. These always win.
 *
 * project_folder is the batch id rather than a bare index, so a run over
 * several input structures cannot produce two campaigns of the same name; the
 * gather step stages every batch into one parent directory and would otherwise
 * collide.
 *
 * campaign_seed must differ per batch. BC2 derives each trajectory's RNG as
 * fold_in(PRNGKey(campaign_seed), trajectory_number) and every batch claims
 * trajectory numbers from its own .campaign_state.json starting at 1, so
 * batches sharing a seed would produce byte-identical designs. campaign_seed is
 * also part of the design-recipe hash that names output files, so distinct
 * seeds keep filenames distinct and make pooling the batches safe. seedIndex is
 * unique across the whole run, not per input structure.
 */
List<String> bc2BatchOverrides(params, String batchId, int seedIndex, int nTraj) {
    return [
        "project_folder=${batchId}".toString(),
        "campaign_seed=${bc2CampaignSeed(params, seedIndex)}".toString(),
        "max_trajectories=${nTraj}".toString(),
    ]
}

/**
 * The campaign_seed a batch runs with: the base plus its run-wide index.
 */
int bc2CampaignSeed(params, int seedIndex) {
    return ((params.bc2_seed ?: 0) as int) + seedIndex
}

/**
 * Arguments for `stage_bc2_config.py stage`, beyond the input and output paths.
 *
 * Hotspot subsampling happens at staging time rather than when the campaign
 * JSON is generated, so it applies in config mode too, where the pipeline never
 * writes the JSON. Each batch draws with its own campaign_seed, so the subsets
 * differ between batches while the run as a whole still reproduces from
 * --bc2_seed.
 *
 * Returns an empty string at the default fraction of 1.0, keeping the staging
 * command identical to what it was before subsampling existed.
 */
String bc2StageArgs(def hotspotSubsample, def seed) {
    def fraction = (hotspotSubsample as Double)
    if (fraction >= 1.0) {
        return ''
    }
    return "--hotspot-subsample ${fraction} --subsample-seed ${seed as Integer}".toString()
}

// BC2's two adaptive mechanisms both read the campaign's own trajectory table,
// so they only engage once a single batch has run enough trajectories:
//   autotune    reviews every AUTOTUNE_REVIEW_TRAJECTORIES = 10 trajectories
//               (parameter_sweep.py:23, a module constant, not a setting)
//   desperation engages after desperation_trajectories = 750 fruitless ones,
//               then one ladder rung per 50 after that (desperation.py:5-6)
// Splitting a run into small batches therefore turns both off. That is usually
// what you want for a fixed-budget comparison, but it is a real difference from
// a single long vanilla BC2 run, so it is reported rather than left implicit.
def bc2AutotuneReviewTrajectories() { return 10 }
def bc2DesperationDefault() { return 750 }

/**
 * Human-readable note on whether autotune and desperation can engage.
 */
String bc2AdaptiveNote(params, int batchSize) {
    def inert = []
    def benchmarkCore = bc2ParamSet(params.bc2_core) && params.bc2_core.toString().trim().equalsIgnoreCase('benchmark')
    if (!benchmarkCore && batchSize < bc2AutotuneReviewTrajectories()) {
        inert << "autotune (needs ${bc2AutotuneReviewTrajectories()} trajectories in one batch)"
    }
    def desperationAfter = bc2DesperationThreshold(params)
    if (batchSize < desperationAfter) {
        inert << "desperation (needs ${desperationAfter} fruitless trajectories in one batch)"
    }
    if (!inert) {
        return ''
    }
    return 'Note - in this mode adaptive tuning and desperation mode will not engage. ' +
           "BC2 adaptive tuning is inert at --bc2_batch_size ${batchSize}: ${inert.join(' and ')}. " +
           'Both read only the batch\'s own trajectory table. Raise --bc2_batch_size to let them engage, ' +
           'or lower the threshold with --bc2_set \'desperation_trajectories=N\'.'
}

/**
 * The effective desperation threshold, honouring a user --bc2_set.
 */
int bc2DesperationThreshold(params) {
    def override = resolveBc2SetList(params.bc2_set).find { entry -> entry.trim().startsWith('desperation_trajectories=') }
    if (override) {
        def raw = override.substring(override.indexOf('=') + 1).trim()
        if (raw ==~ /^\d+$/) {
            return raw.toInteger()
        }
    }
    return bc2DesperationDefault()
}

/**
 * Split a total across n batches, giving the remainder to the earliest batches.
 */
List<Integer> bc2BatchSizes(int total, int batchSize) {
    if (total < 1) {
        total = 1
    }
    if (batchSize < 1) {
        batchSize = 1
    }
    def nBatches = ((total + batchSize - 1).intdiv(batchSize))
    return (0..<nBatches).collect { i -> Math.min(batchSize, total - (i * batchSize)) as Integer }
}

/**
 * Resolved file paths referenced by a BC2 campaign JSON.
 *
 * BC2 resolves `targets[].target_path` and `binder_scaffold` against the
 * campaign JSON's own directory (bindcraft/settings.py read_settings()); those
 * are the only two file-path settings. Returns an empty list for a campaign
 * that only names shipped presets (eg {"target": "hPDL1"}), whose paths resolve
 * inside the container and need no staging.
 *
 * Driver-side JSON parsing is fine here, unlike BoltzGen's YAML configs which
 * need PyYAML and therefore a process.
 */
List<String> extractBc2InputPaths(String configPath, String configDir = null) {
    def pathObj = file(configPath)
    if (!pathObj.exists()) {
        throw new Exception("BC2 campaign config not found: ${configPath}")
    }
    def baseDir = configDir ? file(configDir) : pathObj.parent

    def config = new groovy.json.JsonSlurper().parse(pathObj)
    if (!(config instanceof Map)) {
        throw new Exception("BC2 campaign config must be a JSON object: ${configPath}")
    }

    def raw = []
    (config.targets ?: []).each { target ->
        if (target instanceof Map && target.target_path) {
            raw << target.target_path.toString()
        }
    }
    if (config.binder_scaffold) {
        raw << config.binder_scaffold.toString()
    }

    def seen = [] as LinkedHashSet
    def result = []
    raw.each { p ->
        def resolved = (new File(p).isAbsolute() ? file(p) : file("${baseDir}/${p}")).toString()
        if (!seen.contains(resolved)) {
            seen << resolved
            result << resolved
        }
    }
    return result
}

/**
 * Batch directory / tag name. Mirrors makeBindcraftBatchId() for the v1
 * workflow, but is keyed on the campaign rather than always on an input PDB,
 * since a BC2 campaign may name a shipped target and have no input file.
 */
String makeBc2BatchId(def structure, int index) {
    if (!structure) {
        return "batch_${index}"
    }
    def key = file(structure.toString()).simpleName.replaceAll(/[^A-Za-z0-9]+/, '_')
    return "${key}_${index}"
}

/**
 * How many trailing chains of a BC2 output structure make up the binder.
 *
 * BindCraft2 writes every target chain before the binder: written_chains() in
 * bindcraft/protein.py sorts the complex on `(is_binder_chain(name), name)`,
 * and is_binder_chain is false for targets, so the binder always occupies the
 * last chain letters. The count is therefore enough to identify it, and is
 * more reliable than picking the chain with fewest atoms, which breaks for a
 * multi-chain binder (Fab, scFv, multidomain, a symmetric oligomer) or a
 * target smaller than the binder.
 *
 * `resolved.binder_chains` in campaign_metadata.json holds one entry per binder
 * chain. Returns 0 when it is missing, which the extractor reads as "fall back
 * to the smallest-chain heuristic".
 */
int bc2BinderChainCount(def metadataFile) {
    try {
        def metadata = new groovy.json.JsonSlurper().parse(file(metadataFile.toString()))
        def chains = (metadata instanceof Map) ? metadata.resolved?.binder_chains : null
        return (chains instanceof List && chains.size() > 0) ? chains.size() : 0
    }
    catch (Exception e) {
        log.warn("Could not read binder_chains from ${metadataFile}: ${e.message}. Falling back to the smallest-chain heuristic for FoldSeek queries.")
        return 0
    }
}
