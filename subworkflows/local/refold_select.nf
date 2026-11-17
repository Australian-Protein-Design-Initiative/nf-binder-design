/*
SELECT_FOR_REFOLD: decide which designs are worth refolding.

Every design workflow cofolds its designs once with its own predictor before this
point - AF2 initial guess in rfd and rfd_partial, RF3 in rfd3, BindCraft's and
Germinal's own AF2 passes - and refolding is the held-out second opinion on a subset.
This is where that subset is chosen, from the primary cofold's scores.

    --refold_filter 'ipsae>=0.5;pae_interaction<=10'   threshold gate
    --refold_sort_by='-ipsae'                          rank by one metric
    --refold_max 20                                    cap how many go through

The sort key is written with an '=' because Nextflow reads a leading '-' in a separate
argument as one of its own options and stops with "Unknown option: -ipsae".

The gate is per-design and streams: a design is admitted the moment its own scores
exist, so refolding starts while the rest of the set is still being scored. The cap
does not, with or without a rank - "the first 20" needs a counter shared across the
channel and "the best 20" needs the last design in, and both make the refold wait for
the whole set. Only --refold_filter on its own overlaps the two stages.

Every design that reached here is accounted for in refold_selection.tsv, whether it
was dropped by the gate or fell below the cap, so a thin refold set can be traced to
the threshold or the ranking that caused it.
*/

include {
    parseRefoldFilters
    parseRefoldSortKey
    rankRefoldDesigns
    refoldMaxValue
    refoldNumeric
    refoldPasses
    refoldReportHeader
    refoldReportRow
    refoldSortKeySuggestions
} from '../../modules/local/common/refold_select_utils'

workflow SELECT_FOR_REFOLD {

    take:
    ch_scored       // tuple(meta, structure, Map scores): scores may be [:]
    filters         // String or false: 'ipsae>=0.5;pae_interaction<=10'
    sort_by         // String or false: metric, '-' prefix for descending
    max_designs     // Integer or false: cap
    publish_dir     // String: directory under params.outdir for refold_selection.tsv

    main:

    def parsed_filters = parseRefoldFilters(filters)
    def sort = parseRefoldSortKey(sort_by)
    def max = refoldMaxValue(max_designs)

    if (sort && !max) {
        // Not just cosmetic: with no cap every design goes through, so the ranking
        // decides nothing. Saying so beats silently ignoring a parameter the user
        // expected to narrow the set.
        log.warn(
            "--refold_sort_by '${sort_by}' ranks the designs but --refold_max is unset, " +
            'so every design that passes --refold_filter is refolded anyway and the ' +
            'ranking selects nothing. Set --refold_max to refold only the best of them.'
        )
    }
    if (max && !sort) {
        // Worth saying out loud: the cap looks like "the best N" and is not. rfd3 used
        // to rank by pair_pae_min unless told otherwise, so a run carried over from
        // before --refold_sort_by existed gets arbitrary designs rather than good ones
        // unless it names a metric.
        log.warn(
            "--refold_max ${max} caps the refold at the first ${max} designs to be scored, " +
            'in whatever order they finish. Add --refold_sort_by (for example ' +
            "'-ipsae', or 'pae_interaction' where low is good) to refold the best ${max} instead."
        )
    }

    // The gate. A design carries its own verdict, so this costs no wait.
    def ch_gated = ch_scored.branch { meta, _structure, scores ->
        passed: !parsed_filters || refoldPasses((scores ?: [:]) as Map, parsed_filters, meta.id?.toString())
        rejected: true
    }

    // The cap, with the rank when one is asked for.
    //
    // Only the cap collects. A rank with no cap would change the order of a set that
    // goes through whole, so collecting for it would stall the refold behind the last
    // design to be scored and select nothing in return; that case streams, and the
    // warning above says the ranking is doing nothing.
    //
    // collect/flatMap rather than .take(n): once a channel has more than one consumer
    // Nextflow hands the subworkflow a DataflowStream read adapter, which has no take
    // operator, and the call dies with "Missing process or function take(...)". Whether
    // that happens depends on the order the consumers attach, so .take() works or not
    // depending on which other scoring paths the run enables. flat: false keeps each
    // design a triple rather than flattening it into the list.
    def ch_selected = Channel.empty()
    def ch_considered_rows = Channel.empty()

    if (max) {
        // Ranked but uncapped, so the report can say where the cut fell.
        def ch_ranked = ch_gated.passed
            .collect(flat: false)
            .map { entries -> rankRefoldDesigns(entries as List, sort, null) }
            .map { entries ->
                // A sort key no design carries is almost always a typo or a metric
                // named for a different workflow: BindCraft, for instance, writes
                // Average_i_pTM and 1_i_pTM, never a bare i_pTM. Every comparison then
                // ties, so the cap quietly takes the first N designs in whatever order
                // they finished, which looks exactly like a working ranking in the
                // report. Say so, and name the metrics that are actually there.
                if (sort && entries && !entries.any { entry -> refoldNumeric(((entry[2] ?: [:]) as Map).get(sort.key)) != null }) {
                    log.warn(
                        "--refold_sort_by names '${sort.key}', which no design has a numeric value for, " +
                        "so --refold_max ${max} is taking the first ${max} designs scored rather than the best. " +
                        refoldSortKeySuggestions(entries as List)
                    )
                }
                return entries
            }

        ch_selected = ch_ranked
            .flatMap { entries -> entries.take(max) }
            .map { meta, structure, _scores -> tuple(meta, structure) }

        ch_considered_rows = ch_ranked.flatMap { entries ->
            entries.withIndex().collect { entry, index ->
                refoldReportRow(
                    entry[0].id.toString(),
                    (entry[2] ?: [:]) as Map,
                    sort,
                    index < max,
                    index + 1,
                )
            }
        }
    }
    else {
        ch_selected = ch_gated.passed.map { meta, structure, _scores -> tuple(meta, structure) }
        ch_considered_rows = ch_gated.passed.map { meta, _structure, scores ->
            refoldReportRow(meta.id.toString(), (scores ?: [:]) as Map, sort, true, null)
        }
    }

    def ch_rejected_rows = ch_gated.rejected.map { meta, _structure, scores ->
        refoldReportRow(meta.id.toString(), (scores ?: [:]) as Map, sort, false, null)
    }

    def ch_report = ch_considered_rows
        .mix(ch_rejected_rows)
        .collectFile(
            name: 'refold_selection.tsv',
            storeDir: "${params.outdir}/${publish_dir}",
            seed: refoldReportHeader(),
            newLine: true,
            sort: false,
        )

    emit:
    selected = ch_selected
    report = ch_report
}
