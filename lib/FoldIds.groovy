/*
Shared id sanitisation + collision checks for fold / fold_pulldown. Returns
error lists so callers can raise them via error() (Groovy classes cannot).
*/
class FoldIds {

    /** Same rule fold_pulldown_msa.nf uses for filenames/meta.id: keep the
     * FASTA header id readable as a path component, replace everything else. */
    static String sanitize(id) {
        return id.toString().replaceAll(/[^a-zA-Z0-9_.-]/, "_")
    }

    /** First whitespace-delimited token of each '>' header line in a FASTA file. */
    static List<String> extractFastaIds(f) {
        def file = (f instanceof File) ? f : f.toFile()
        def ids = []
        file.eachLine { line ->
            def l = line.trim()
            if (l.startsWith('>')) {
                ids << l.substring(1).split(/\s+/)[0]
            }
        }
        return ids
    }

    private static List<String> duplicated(List<String> ids) {
        def counts = ids.countBy { it }
        return counts.findAll { _id, n -> n > 1 }.keySet().sort()
    }

    /**
     * fold_pulldown id validation. Checks (raw ids, before sanitisation):
     *  - duplicate ids within targets, or within binders
     *  - an id present in both targets and binders (same-named per-chain MSA
     *    files/dirs collide in every engine but AF2, which keys A/B by role)
     * and (sanitised ids, as actually used for filenames/meta.id):
     *  - sanitisation collisions within targets, or within binders
     *  - pair-id collisions: pair_id = "${target}_and_${binder}" is ambiguous
     *    whenever an id itself contains '_and_' (e.g. target "a_and_b" +
     *    binder "c" collides with target "a" + binder "b_and_c").
     * Returns a list of error strings (empty if everything is fine).
     */
    static List<String> validatePulldownIds(List<String> targetIds, List<String> binderIds) {
        def errors = []

        def dupTargets = duplicated(targetIds)
        if (dupTargets) {
            errors << "duplicate target id(s): ${dupTargets.join(', ')}"
        }
        def dupBinders = duplicated(binderIds)
        if (dupBinders) {
            errors << "duplicate binder id(s): ${dupBinders.join(', ')}"
        }

        def targetSet = targetIds as Set
        def binderSet = binderIds as Set
        def sharedRaw = (targetSet.intersect(binderSet)).sort()
        if (sharedRaw) {
            errors << (
                "id(s) present in both targets and binders: ${sharedRaw.join(', ')} " +
                "(per-chain MSA files/dirs are keyed by this id and collide across " +
                "every engine except AF2 - rename one side)"
            )
        }

        def sanTargets = targetIds.collect { sanitize(it) }
        def sanBinders = binderIds.collect { sanitize(it) }
        def dupSanTargets = duplicated(sanTargets) - dupTargets
        if (dupSanTargets) {
            errors << (
                "target id(s) collide after sanitisation (non [a-zA-Z0-9_.-] " +
                "characters replaced with '_'): ${dupSanTargets.join(', ')}"
            )
        }
        def dupSanBinders = duplicated(sanBinders) - dupBinders
        if (dupSanBinders) {
            errors << (
                "binder id(s) collide after sanitisation (non [a-zA-Z0-9_.-] " +
                "characters replaced with '_'): ${dupSanBinders.join(', ')}"
            )
        }

        def pairCounts = [:]
        sanTargets.each { t ->
            sanBinders.each { b ->
                def pid = "${t}_and_${b}"
                pairCounts[pid] = (pairCounts[pid] ?: 0) + 1
            }
        }
        def ambiguousPairs = pairCounts.findAll { _pid, n -> n > 1 }.keySet().sort()
        if (ambiguousPairs) {
            errors << (
                "pair id(s) are ambiguous - more than one (target, binder) combination " +
                "sanitises to the same '<target>_and_<binder>' id: ${ambiguousPairs.join(', ')} " +
                "(an id containing the literal '_and_' is the usual cause)"
            )
        }

        return errors
    }

    /** --method fold: fail fast on duplicate meta.id (= FASTA file baseName). */
    static List<String> validateFoldIds(List<String> ids) {
        def dup = duplicated(ids)
        if (!dup) {
            return []
        }
        return [
            "duplicate id(s) derived from input FASTA filenames: ${dup.join(', ')} " +
            "(two input files share a basename, e.g. 'x.fa' and 'x.fasta' - rename one)"
        ]
    }
}
