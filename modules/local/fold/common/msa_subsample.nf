/*
Parse --msa_subsample for fold.nf and fan out CF-random-style shallow MSA depths.
*/

def msaSubsampleDepthDefaults() {
    return '1:2,2:4,4:8,8:16,16:32,32:64,64:128'
}

boolean msaSubsampleIsEnabled(v) {
    if (v == null) {
        return false
    }
    if (v instanceof Boolean) {
        return v
    }
    def s = v.toString().trim().toLowerCase()
    if (s in ['false', '0', '', 'null', 'none']) {
        return false
    }
    return true
}

boolean msaSubsampleIncludeFull(v) {
    if (v == null) {
        return true
    }
    if (v instanceof Boolean) {
        return v
    }
    def s = v.toString().trim().toLowerCase()
    return !(s in ['false', '0', '', 'null', 'none'])
}

// Count '>' headers in an a3m (sequences, including query).
int msaSubsampleCountA3mSequences(path) {
    def n = 0
    path.withReader { r ->
        r.eachLine { line ->
            if (line.startsWith('>')) {
                n += 1
            }
        }
    }
    return n
}

// List of [max_seq, max_extra_seq] ints; empty when subsample is off.
List msaSubsampleParseDepths(v) {
    if (!msaSubsampleIsEnabled(v)) {
        return []
    }
    def raw
    if (v instanceof Boolean || v.toString().trim().toLowerCase() in ['true', '1', 'yes']) {
        raw = msaSubsampleDepthDefaults()
    }
    else {
        raw = v.toString().trim()
    }
    def depths = []
    raw.split(',').each { part ->
        def p = part.trim()
        if (!p) {
            return
        }
        def bits = p.split(':')
        if (bits.size() != 2) {
            throw new IllegalArgumentException(
                "Invalid --msa_subsample depth '${p}' (expected max_seq:max_extra_seq)"
            )
        }
        def maxSeq
        def maxExtra
        try {
            maxSeq = bits[0].trim() as int
            maxExtra = bits[1].trim() as int
        }
        catch (NumberFormatException e) {
            throw new IllegalArgumentException(
                "Invalid --msa_subsample depth '${p}' (non-integer max_seq/max_extra_seq)"
            )
        }
        if (maxSeq < 1 || maxExtra < 0) {
            throw new IllegalArgumentException(
                "Invalid --msa_subsample depth '${p}' (max_seq>=1, max_extra_seq>=0)"
            )
        }
        depths << [maxSeq, maxExtra]
    }
    if (!depths) {
        throw new IllegalArgumentException('--msa_subsample enabled but no depths parsed')
    }
    return depths
}

/*
Depth jobs for fan-out: null = full MSA (no subsample); otherwise [max, extra].
When subsample is off, returns [null] (single unsubsampled job).
When nSeq is set, drop shallow depths with max_seq >= nSeq (noop full-MSA
shuffles); always keep full when include_full is on.
*/
List msaSubsampleDepthJobs(msa_subsample, msa_subsample_include_full, Integer nSeq = null) {
    def depths = msaSubsampleParseDepths(msa_subsample)
    if (!depths) {
        return [null]
    }
    if (nSeq != null) {
        depths = depths.findAll { d -> d[0] < nSeq }
    }
    def jobs = []
    if (msaSubsampleIncludeFull(msa_subsample_include_full)) {
        jobs << null
    }
    jobs.addAll(depths)
    if (!jobs) {
        throw new IllegalArgumentException(
            "All --msa_subsample depths have max_seq >= n_seq (${nSeq}); " +
            "enable --msa_subsample_include_full or choose shallower depths"
        )
    }
    return jobs
}

// Stable non-negative seed for -resume (no random draw at parse time).
// toString() so this is Java's String hash. A GString's hashCode is 37 plus
// that, and Nextflow 26.10's interpolation is a String, so hashing the
// GString directly changed the seed between versions and broke -resume.
int msaSubsampleStableSeed(String id, int batchOrRun, int maxSeq, int maxExtra) {
    def h = "${id}|${batchOrRun}|${maxSeq}|${maxExtra}".toString().hashCode()
    return h == Integer.MIN_VALUE ? 0 : Math.abs(h)
}

String msaSubsampleDepthTag(maxSeq, maxExtra) {
    if (maxSeq == null) {
        return 'full'
    }
    return "${maxSeq}_${maxExtra}"
}
