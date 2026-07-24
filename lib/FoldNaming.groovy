/*
Single source of truth for the flat <outdir>/fold/predictions/ file naming used
by fold.nf. Both a folding module's second publishDir (which copies the
structure into fold/predictions/) and the scoring subworkflow (which records the
published name in the fold_scores TSV's predictions_file column) build the name
from here, so the two can never drift.

The published name is <prefix><native-basename>, where <prefix> carries the tool
tag plus - when fold.nf fans --n_predictions across jobs or --msa_subsample
across depths (meta.fold_namespaced) - a batch/depth qualifier that keeps names
unique across tasks. AF2 always carries its run number.
*/

class FoldNaming {

    // "_msa<tag>" when a subsample depth job set meta.msa_depth_tag, else "".
    static String msaBit(Map meta) {
        meta.msa_depth_tag ? "_msa${meta.msa_depth_tag}" : ''
    }

    // Flat-gather filename prefix for boltz / rf3 / protenix (append the native
    // structure basename). Namespaced jobs get a batch qualifier.
    static String flatPrefix(String tool, Map meta) {
        def msa = msaBit(meta)
        return meta.fold_namespaced \
            ? "${tool}_batch${meta.fold_batch}${msa}_" \
            : "${tool}${msa}_"
    }

    // AF2 prefix (append the native structure basename). AF2 always encodes its
    // per-run index; depth jobs add the msa bit.
    static String af2Prefix(Map meta) {
        "af2_${meta.id}_run${meta.af2_run}${msaBit(meta)}_"
    }
}
