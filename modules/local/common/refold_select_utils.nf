/*
Choosing which designs go through to refold scoring.

Refolding is the expensive, held-out second opinion: every design has already been
cofolded once by the pipeline's own predictor (AF2 initial guess in rfd, RF3 in rfd3,
BindCraft's and Germinal's own AF2 passes), and refolding asks a different set of
predictors the same question about a subset worth the GPU time. So what gets filtered
on is the primary cofold's scores, which is why selection takes a score map per design
and reads nothing out of the structure itself.

Three controls, and they compose:

  --refold_filter 'ipsae>=0.5;pae_interaction<=10'   a threshold gate: every design
                                                     that passes goes through
  --refold_sort_by '-ipsae'                          rank by one metric
  --refold_max 20                                    cap how many go through

The gate is per-design, so it costs nothing to wait for: a design is admitted as soon
as its own scores exist, and refolding starts while the rest are still being scored.
Ranking is not. "The best 20" is not knowable until the last design is in, so
--refold_sort_by makes the refold wait for the whole set. That trade-off is why they
are separate parameters rather than one.

--refold_max without --refold_sort_by keeps streaming and takes the first N to arrive.
That is arbitrary with respect to quality, but it is what a user who wants "some of
them, cheaply" is asking for, and it is what rfd4_proteina did before this existed.

These are plain functions rather than a lib/ class so that nf-test can drive them
directly (see tests/modules/local/common/refold_select_utils.nf.test); lib/ is not on
the classpath for a test script that lives outside the project root.
*/

// Operators, matching bin/filter_designs.py so one filter syntax covers both. Longest
// first: matched the other way round, '>' would claim the '>' of '>=' and leave '='.
def refoldFilterOperators() {
    return ['>=', '<=', '==', '!=', '>', '<', '=']
}

/*
Split a filter spec into [metric, op, value] maps.

Semicolon-separated, as everywhere else filters are written in this pipeline
(--rfd_filters, --rfd3_filters). An empty or falsy spec is no filters at all, not a
filter that rejects everything.
*/
def parseRefoldFilters(spec) {
    if (!spec || spec.toString().trim().isEmpty()) {
        return []
    }
    def expression = ~/^\s*([A-Za-z0-9_.]+)\s*(>=|<=|==|!=|>|<|=)\s*(.+?)\s*$/
    return spec.toString().split(';')
        .collect { it.trim() }
        .findAll { !it.isEmpty() }
        .collect { text ->
            def match = expression.matcher(text)
            if (!match.matches()) {
                throw new IllegalArgumentException(
                    "Cannot parse filter expression '${text}'. Expected <metric><operator><value>, " +
                    "with the operator one of ${refoldFilterOperators().join(', ')} - for example 'ipsae>=0.5'."
                )
            }
            return [metric: match.group(1), op: match.group(2), value: match.group(3)]
        }
}

/*
A score value as a Double, or null where it is not a number.

Score tables write missing values in several ways depending on which tool produced
them - an empty cell, the string 'None' from a pandas round trip, 'NA', 'nan' - and
all of them mean the same thing here.
*/
def refoldNumeric(value) {
    if (value == null) {
        return null
    }
    if (value instanceof Number) {
        def d = (value as Number).doubleValue()
        return Double.isNaN(d) ? null : d
    }
    def text = value.toString().trim()
    if (text.isEmpty() || text.toLowerCase() in ['none', 'na', 'nan', 'null']) {
        return null
    }
    try {
        def d = Double.parseDouble(text)
        return Double.isNaN(d) ? null : d
    }
    catch (NumberFormatException ignored) {
        return null
    }
}

// Compare one score against one parsed filter, coercing as bin/filter_designs.py does.
def refoldCompare(score, String op, String expected) {
    def wanted = refoldNumeric(expected)
    if (wanted != null) {
        def actual = refoldNumeric(score)
        // A design whose metric is missing fails every numeric comparison, '!=' included:
        // letting it pass != would admit exactly the designs nothing is known about.
        if (actual == null) {
            return false
        }
        switch (op) {
            case '>': return actual > wanted
            case '<': return actual < wanted
            case '>=': return actual >= wanted
            case '<=': return actual <= wanted
            case '==':
            case '=': return actual == wanted
            case '!=': return actual != wanted
        }
    }

    def lowered = expected.toLowerCase()
    if (lowered in ['true', 'yes', 'false', 'no']) {
        def wantedBool = lowered in ['true', 'yes']
        def actualBool = score == null
            ? false
            : (score instanceof Boolean ? score : score.toString().toLowerCase() in ['true', 'yes', '1'])
        switch (op) {
            case '==':
            case '=': return actualBool == wantedBool
            case '!=': return actualBool != wantedBool
            default:
                throw new IllegalArgumentException(
                    "Operator '${op}' does not apply to the boolean value '${expected}'."
                )
        }
    }

    def actual = score == null ? '' : score.toString()
    switch (op) {
        case '==':
        case '=': return actual == expected
        case '!=': return actual != expected
        case '>': return actual > expected
        case '<': return actual < expected
        case '>=': return actual >= expected
        case '<=': return actual <= expected
    }
    return false
}

/*
Whether a design's scores satisfy every filter.

A metric the score table does not have is a mistake in the filter rather than a design
that fails it, so it stops the run and names what is available. Silently rejecting
everything would look like a threshold that is merely too tight, and the run would
report "0 designs refolded" with nothing to say why.
*/
def refoldPasses(Map scores, List filters, String design_id = null) {
    return filters.every { f ->
        if (!scores.containsKey(f.metric)) {
            throw new IllegalArgumentException(
                "Filter metric '${f.metric}' is not a column in the scores" +
                (design_id ? " for design '${design_id}'" : '') + '. Available: ' +
                (scores.keySet().sort().join(', ') ?: '(no scores at all)')
            )
        }
        return refoldCompare(scores[f.metric], f.op, f.value)
    }
}

/*
Split a sort spec into its metric and direction.

A leading '-' means descending, as in --refold_sort_by '-ipsae' for "highest ipSAE
first". Metrics where low is good - pae_interaction, rmsd - need no prefix.
*/
def parseRefoldSortKey(spec) {
    if (!spec || spec.toString().trim().isEmpty()) {
        return null
    }
    def text = spec.toString().trim()
    def descending = text.startsWith('-')
    def key = descending ? text.substring(1).trim() : text
    if (key.isEmpty()) {
        throw new IllegalArgumentException("--refold_sort_by '${spec}' names no metric.")
    }
    return [key: key, descending: descending]
}

/*
Order two designs by a sort key.

Designs with no value for the metric go last whichever direction is asked for: an
unscored design is not a good one, and sorting it to the front would have --refold_max
select precisely the designs nothing is known about.
*/
def compareRefoldScores(Map a, Map b, Map sort) {
    def valA = refoldNumeric(a?.get(sort.key))
    def valB = refoldNumeric(b?.get(sort.key))
    if (valA == null && valB == null) {
        return 0
    }
    if (valA == null) {
        return 1
    }
    if (valB == null) {
        return -1
    }
    return sort.descending ? (valB <=> valA) : (valA <=> valB)
}

/*
Apply ranking and the cap to designs already through the gate.

Entries are [meta, structure, scores] triples; the return keeps that shape and that
order, so the caller can report on what came back.
*/
def rankRefoldDesigns(List entries, Map sort, Integer max) {
    def ordered = sort
        ? entries.toSorted { a, b -> compareRefoldScores((a[2] ?: [:]) as Map, (b[2] ?: [:]) as Map, sort) }
        : entries
    return (max != null && max > 0) ? ordered.take(max) : ordered
}

/*
Drop the trailing _model<N> that a predictor adds to a structure filename.

BindCraft writes an accepted design to <design>_model<N>.pdb, where N is which of its
five AF2 models produced the structure, while final_design_stats.csv names the same
design in its Design column without that suffix. Joining the structures to their
scores on the raw filename therefore matches nothing: every design arrives at the
refold with an empty score map, so --refold_filter admits all of them and
--refold_sort_by ranks none of them, with no error anywhere to say so.

The design id itself keeps the suffix - two models of one design are two structures and
have to stay distinguishable - so this is only for the join key.
*/
def refoldStripModelSuffix(name) {
    return name?.toString()?.replaceFirst(/_model\d+$/, '')
}

/*
Metric names a run could have sorted by, for the message when --refold_sort_by names
one that is not there.

Only columns with a numeric value somewhere in the set are offered: a design's scores
map also carries its sequence, its notes and its id, and suggesting those as sort keys
would send the reader down a second wrong turn. The list is capped because a refold
score table runs to sixty-odd columns and a warning nobody reads to the end helps
nobody.
*/
def refoldSortKeySuggestions(List entries, int limit = 12) {
    def numeric = [] as Set
    entries.each { entry ->
        ((entry[2] ?: [:]) as Map).each { key, value ->
            if (refoldNumeric(value) != null) {
                numeric << key.toString()
            }
        }
    }
    if (!numeric) {
        return 'No design carries any numeric score, so nothing here can be ranked.'
    }
    def names = numeric.toList().sort()
    def shown = names.take(limit)
    def tail = names.size() > limit ? ", and ${names.size() - limit} more" : ''
    return "Metrics available: ${shown.join(', ')}${tail}."
}

// --refold_max as a positive Integer, or null for "no cap".
def refoldMaxValue(max) {
    if (max == null || max == false) {
        return null
    }
    def value = max instanceof Number
        ? (max as Number).intValue()
        : Integer.parseInt(max.toString().trim())
    return value > 0 ? value : null
}

/*
Whether the whole design set has to be in hand before refolding can start.

The cap is what forces it, with or without a rank: "the best N" needs the last design
scored, and "the first N" needs a count shared across the channel. A threshold gate
admits designs one at a time, so refolding overlaps the scoring of the designs still
to come, and a rank with no cap changes only the order of a set that goes through
whole - which is why that case keeps streaming rather than collecting for nothing.
*/
def refoldNeedsCompleteSet(sort_by, max) {
    return refoldMaxValue(max) != null
}

// The selection report's header, matching refoldReportRow.
def refoldReportHeader() {
    return ['design_id', 'refold_selection', 'refold_rank', 'refold_sort_metric', 'refold_sort_value'].join('\t')
}

/*
The --refold_* section of a workflow's --help, so that six workflows describe the same
parameters the same way.

Indented to sit inside a .stripIndent() help string at the usual eight spaces; a
caller splices it in with ${refoldHelpText(params)}.
*/
def refoldHelpText(Map params) {
    return """Refold scoring (--refold_with):
            A held-out second opinion on the designs this run has already cofolded with
            its own predictor. Folds each design's binder sequence against its target
            with the engines you name, compares the result with the design, and scores
            the rebuilt interface. Every engine gives the same columns in
            combined_scores.tsv, prefixed with the engine name: confidence (iptm,
            ipSAE, plddt, ptm), RMSD against the design (target-aligned binder and
            whole complex), the binder folded on its own, and BindCraft interface
            metrics. A cross-engine consensus_z ranks designs when more than one engine
            is named.

            --refold_with          Comma-separated engines: af2, af2_mono, boltz, rf3,
                                   protenix, af3, openfold3, chai1, esmfold2,
                                   esmfold2_fast [default: ${params.refold_with}]

            Which designs are refolded, from the first cofold's scores. Both are
            optional, and refold_selection.tsv records what happened to every design:
            --refold_filter        Threshold gate, semicolon-separated, same syntax as
                                   --rfd_filters, eg 'ipsae>=0.5;pae_interaction<=10'.
                                   Per-design, so refolding starts while the rest of the
                                   set is still being scored [default: ${params.refold_filter}]
            --refold_sort_by       Rank by one metric before --refold_max takes the best
                                   of them; '-' prefix for descending. Write it with an
                                   '=', as --refold_sort_by='-ipsae', or Nextflow reads
                                   the leading '-' as one of its own options. Cannot
                                   stream: the refold waits for every design to be
                                   scored [default: ${params.refold_sort_by}]
            --refold_max           Cap on designs refolded. With --refold_sort_by, the
                                   best N; without it, the first N to be scored
                                   [default: ${params.refold_max}]

            --refold_metrics       Summary columns carried into combined_scores.tsv
                                   [default: ${params.refold_metrics}]
            --refold_structure_stat  How several samples per design reduce to one RMSD
                                   and one interface score: selected|min|first|mean|max
                                   [default: ${params.refold_structure_stat}, both read
                                   from the best-ranked sample]
            --refold_binder_monomer  Also fold each binder alone, for the monomer-vs-complex
                                   and monomer-vs-design comparisons [default: ${params.refold_binder_monomer}]
            --refold_monomer_single_sequence  Fold that monomer from sequence alone, with no
                                   MSA search [default: ${params.refold_monomer_single_sequence}]
            --refold_interface_scoring  BindCraft interface metrics on the refolded complex
                                   [default: ${params.refold_interface_scoring}]
            --output_rmsd_aligned  Keep the superimposed structures from the RMSD
                                   comparisons [default: ${params.output_rmsd_aligned}]
            --n_predictions        Structures per design per engine, run as that many
                                   seeds of one sample each [default: ${params.n_predictions}]
            --msa_method           jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]
            --create_target_msa    Build an MSA for the target [default: ${params.create_target_msa}]
            --use_msa_server       Let the engine fetch its own MSA [default: ${params.use_msa_server}]
            --consensus_metric / --z_stat / --z_scope / --min_pool
                                   Ranking, as in --method fold_pulldown
            --colabfold_envdb      ColabFold environment database (local MSA) [default: ${params.colabfold_envdb}]
            --uniref30             UniRef30 database (local MSA) [default: ${params.uniref30}]

            Every engine flag from --method fold --help (--af2_*, --boltz_*,
            --protenix_*, --af3_*, --openfold3_*, --chai1_*, --esmfold2_*) applies here
            too. RF3's are the exception: use --refold_rf3_* (ckpt_path, num_steps,
            n_recycles, batch_size, early_stopping_plddt_threshold, seed), since plain
            --rf3_* names are taken by the design-time RF3 step in other workflows.

            Renamed, still accepted with a warning: --refold_oracles and
            --refold_oracles_max are --refold_with and --refold_max;
            --refold_af2ig_filters is --refold_filter. Renamed and now an error:
            --refold_oracle_metrics is --refold_metrics, --refold_oracle_target_id is
            --refold_target_id, and --boltz_refold_batch_size is
            --refold_boltz_batch_size (samples per fold job, not designs per task).
            Refolding takes its MSA and template choices from the fold stack, so
            --refold_create_target_msa, --refold_use_msa_server, --refold_target_fasta,
            --refold_target_templates and --refold_alignment are gone; use the fold
            stack's own flags above."""
}

/*
Which column of a score table names the design.

Workflows that hand REFOLD a table written by the design tool rather than by this
pipeline - BindCraft's final_design_stats.csv, Germinal's accepted_designs.csv - do not
all spell the identifier the same way, and the spelling can change with the tool. The
candidates are tried in order and the first one the table has wins.

Failing here, naming the columns that are actually present, beats returning null: the
join that follows would drop every design, and the run would report nothing refolded
with no indication that the key was the problem.
*/
def refoldIdColumn(Map row, List candidates, String workflow_name) {
    def found = candidates.find { row.containsKey(it) }
    if (!found) {
        throw new IllegalArgumentException(
            "${workflow_name}: cannot tell which column of the score table names each design. " +
            "Looked for ${candidates.join(', ')}; the table has ${row.keySet().sort().join(', ')}. " +
            'Refolding needs it to match each score row to its structure.'
        )
    }
    return found
}

/*
Stop a run that uses a parameter name refolding no longer has.

Refolding used to be several paths with a parameter set each - --refold_oracles beside
--full_refold_with, --refold_af2ig_filters in rfd only, MSA and template options of its
own - and they are now one path with one set of names. Every old name is still declared
(conf/refold.config) so that setting one fails here with its replacement named, rather
than drawing Nextflow's generic unknown-parameter warning, which is easy to lose in a
long log while the run carries on with the default.
*/
def checkRefoldParamRenames(Map params, String workflow_name) {
    [
        refold_use_msa_server: '--use_msa_server',
        refold_create_target_msa: '--create_target_msa',
        refold_target_templates: '--templates',
        refold_alignment: '--create_target_msa with --msa_method mmseqs2_colabfold',
        refold_target_fasta: 'a target sequence in the design structures themselves',
    ].each { name, replacement ->
        if (params[name]) {
            throw new IllegalArgumentException(
                "${workflow_name}: --${name} has been replaced. Refolding takes its MSA and " +
                "template choices from the fold stack now; use ${replacement}."
            )
        }
    }

    [
        refold_oracle_metrics: '--refold_metrics',
        refold_oracle_target_id: '--refold_target_id',
    ].each { name, replacement ->
        if (params[name]) {
            throw new IllegalArgumentException(
                "${workflow_name}: --${name} has been renamed; use ${replacement}."
            )
        }
    }

    if (params.boltz_refold_batch_size) {
        throw new IllegalArgumentException(
            "${workflow_name}: --boltz_refold_batch_size is gone. It batched several designs " +
            'into one Boltz refold comparison task, and the refold is not Boltz-specific any ' +
            'more; use --refold_boltz_batch_size, which is samples per fold job.'
        )
    }
}

/*
The engines to refold with, from --refold_with.

--refold_oracles named the same thing and is still accepted, folded into the list
rather than running a second path beside it.
*/
def refoldEngines(Map params, String workflow_name, log_handle) {
    def engines = FoldValidation.parseMethods(params.refold_with ?: '')
    if (params.refold_oracles) {
        log_handle.warn(
            "${workflow_name}: --refold_oracles is deprecated and will be removed; it is the " +
            "same thing as --refold_with, which takes a list of engines. Adding " +
            "'${params.refold_oracles}' to it."
        )
        engines = (engines + FoldValidation.parseMethods(params.refold_oracles)).unique()
    }
    return engines
}

/*
The selection controls as one map, for SELECT_FOR_REFOLD.

Deprecated spellings fold in here so that each workflow's call site reads the same.
--refold_af2ig_filters was rfd's name for the threshold gate back when only rfd had
one; --refold_oracles_max was a second name for the cap.
*/
def refoldSelection(Map params, String workflow_name, log_handle) {
    def filter = params.refold_filter
    if (params.refold_af2ig_filters) {
        log_handle.warn(
            "${workflow_name}: --refold_af2ig_filters is deprecated; use --refold_filter, " +
            'which selects designs for refolding in every workflow.'
        )
        filter = filter ?: params.refold_af2ig_filters
    }

    def max = params.refold_max
    if (params.refold_oracles_max) {
        log_handle.warn("${workflow_name}: --refold_oracles_max is deprecated; use --refold_max.")
        max = max ?: params.refold_oracles_max
    }

    return [filter: filter, sort_by: params.refold_sort_by, max: max]
}

// One row of the selection report: what was known about a design and what happened to it.
def refoldReportRow(String design_id, Map scores, Map sort, boolean selected, Integer rank_position) {
    def sort_value = sort ? (refoldNumeric(scores?.get(sort.key))?.toString() ?: '') : ''
    return [
        design_id,
        selected ? 'selected' : 'rejected',
        rank_position == null ? '' : rank_position.toString(),
        sort ? sort.key : '',
        sort_value,
    ].join('\t')
}
