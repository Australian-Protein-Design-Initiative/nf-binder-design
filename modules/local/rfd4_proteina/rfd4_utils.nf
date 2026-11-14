// Utility functions for the RFD4-Proteina workflow
//
// The contig chain rules here mirror bin/rfd4_proteina/rfd4_config.py `infer-chains`, which is
// what the tests exercise. Keep the two in step if the rules change.

/**
 * Split a contig into chain runs on /0 chainbreaks, dropping empty runs.
 */
List<List<String>> splitRfd4ContigRuns(String contig) {
    def runs = [[]]
    contig.split(',').each { raw ->
        def token = raw.trim()
        if (!token) {
            return
        }
        if (token == '/0') {
            runs << []
        }
        else {
            runs[-1] << token
        }
    }
    return runs.findAll { it }
}

/**
 * Chain id of a preserved segment token (e.g. "A17-131" -> "A"), or null for an expandable
 * segment (a bare length or length range, optionally P/R/D suffixed, or a "name:" prefixed one).
 */
String rfd4PreservedChain(String token) {
    if (token.contains(':') || token ==~ /^\d+(-\d+)?[PRD]?$/) {
        return null
    }
    def matcher = token =~ /^([A-Za-z][A-Za-z0-9_]*?)(\d+)(-(\d+))?$/
    if (!matcher.matches()) {
        throw new Exception("Could not parse contig token: '${token}'")
    }
    return matcher[0][1]
}

/**
 * Target and binder chain ids for a two-chain binder contig.
 *
 * A run containing a preserved segment takes that segment's chain id; a run with none is the
 * designed binder, which RFD4-Proteina gives the first unused letter.
 */
List<String> rfd4ChainsFromContig(String contig) {
    def runs = splitRfd4ContigRuns(contig)
    if (!runs) {
        throw new Exception('Empty contig')
    }

    def runChains = runs.collect { run ->
        run.collect { rfd4PreservedChain(it) }.find { it != null }
    }

    def targetChains = runChains.findAll { it != null }
    if (!targetChains) {
        throw new Exception(
            'Contig has no preserved segment, so there is no target chain. A contig made entirely ' +
            'of expandable segments is rejected by RFD4-Proteina.'
        )
    }

    def used = targetChains.toSet()
    def binderChains = []
    runChains.each { chain ->
        if (chain != null) {
            return
        }
        def fresh = ('A'..'Z').find { !used.contains(it) }
        used << fresh
        binderChains << fresh
    }

    if (!binderChains) {
        throw new Exception(
            'Contig has no designed chain: every run contains a preserved segment. Add an ' +
            'expandable segment, e.g. "50-120, /0, A17-131".'
        )
    }
    return [targetChains[0], binderChains[0]]
}

/**
 * A param's trimmed string value, or null when it is unset.
 *
 * Params default to boolean false here rather than null, and Groovy renders that as the string
 * "false", which is truthy. A bare `params.x?.toString()` therefore reports an unset param as the
 * value "false" instead of as absent.
 */
String rfd4UnsetToNull(value) {
    if (value == null || value instanceof Boolean && !value) {
        return null
    }
    def text = value.toString().trim()
    return (text && text.toLowerCase() != 'false' && text.toLowerCase() != 'null') ? text : null
}

/**
 * Target and binder chain ids for the run, in that order.
 *
 * Explicit --rfd4_target_chain / --rfd4_binder_chain win. Otherwise they are derived from
 * --contigs in params mode. In config mode the spec may be YAML, which cannot be parsed on
 * the Nextflow head node, so the defaults apply and the user overrides them if their target
 * chain is not A.
 */
List<String> resolveRfd4Chains(params) {
    // These params default to boolean false, and false.toString() is the non-empty - therefore
    // truthy - string "false". Taken at face value that makes every chain id literally "false",
    // so the derivation below never runs. Treat an unset param as unset.
    def explicitTarget = rfd4UnsetToNull(params.rfd4_target_chain)
    def explicitBinder = rfd4UnsetToNull(params.rfd4_binder_chain)
    if (explicitTarget && explicitBinder) {
        return [explicitTarget, explicitBinder]
    }

    def derived = null
    def contigs = rfd4UnsetToNull(params.contigs)
    if (contigs) {
        derived = rfd4ChainsFromContig(contigs)
    }
    else {
        log.info(
            'Assuming target chain A and binder chain B. Set --rfd4_target_chain / ' +
            '--rfd4_binder_chain if your condition spec uses other chain ids.'
        )
        derived = ['A', 'B']
    }

    return [explicitTarget ?: derived[0], explicitBinder ?: derived[1]]
}
