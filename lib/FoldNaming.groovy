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

    // The engine's own label for one predicted sample, as FOLD_PARSE_CONFIDENCE records
    // it in the `model` column of the score tables.
    //
    // It is what tells two samples of one design apart within a fold job, and the refold
    // needs the same string so its RMSD and interface rows join the confidence rows on
    // (id, tool, batch, model). An engine given a list of seeds returns several samples
    // from one job, so the batch alone no longer separates them.
    //
    // The argument is the model's path, not just its basename: Protenix is the one
    // engine whose sample label cannot be read from the basename alone.
    static String sampleLabel(String tool, path) {
        def full = path.toString()
        def bn = full.replaceFirst(/^.*\//, '')
        switch (tool) {
            case 'boltz':
                def m = (bn =~ /_model_(\d+)\./)
                return m ? m[0][1] : bn.replaceFirst(/\.[^.]+$/, '')
            case 'rf3':
            case 'af3':
                def m = (bn =~ /_(seed-\d+_sample-\d+)_/)
                return m ? m[0][1] : bn.replaceFirst(/\.[^.]+$/, '')
            case 'openfold3':
            case 'chai1':
            case 'esmfold2':
            case 'esmfold2_fast':
                def m = (bn =~ /_(seed_\d+_sample_\d+)_/)
                return m ? m[0][1] : bn.replaceFirst(/\.[^.]+$/, '')
            case 'protenix':
                def m = (bn =~ /_(sample_\d+)\./)
                return m ? "${protenixSeedBit(full)}${m[0][1]}".toString() : bn.replaceFirst(/\.[^.]+$/, '')
            default:
                // AF2 and anything else: the basename stem, which is unique within a run.
                return bn.replaceFirst(/\.[^.]+$/, '')
        }
    }

    // "seed_<N>_" read back out of a Protenix output path, else "".
    //
    // Every other engine puts the seed in the file name. Protenix puts it only in the
    // directory - <name>/seed_<N>/predictions/<name>_sample_<r>.cif - so several seeds
    // from one invocation share a basename, and a flat gather of them would collide.
    // Recovering it from the path is what lets --protenix_seeds take a list.
    static String protenixSeedBit(path) {
        def m = (path.toString() =~ /(?:^|\/)seed_(\d+)\//)
        return m ? "seed_${m[0][1]}_" : ''
    }

    // AF2 prefix (append the native structure basename). AF2 always encodes its
    // per-run index; depth jobs add the msa bit. meta.af2_tool distinguishes the
    // multimer engine ('af2') from the monomer chain-break one ('af2_mono') so the
    // two can be selected together without colliding in fold/predictions/.
    static String af2Prefix(Map meta) {
        "${meta.af2_tool ?: 'af2'}_${meta.id}_run${meta.af2_run}${msaBit(meta)}_"
    }

    // AlphaFold3's Input.sanitised_name(): the job name AF3 uses for its output
    // directory and file prefixes (make_af3_input.py writes the same value).
    static String af3Name(id) {
        id.toString().replace(' ', '_').replaceAll(/[^A-Za-z0-9_.-]/, '')
    }

    // OpenFold3 query key / output directory name (make_openfold3_input.py
    // sanitises with the same rule as AF3).
    static String openfold3Name(id) {
        af3Name(id)
    }

    // Which basenames in an engine's output are predicted model structures, as
    // opposed to MSAs, confidence JSONs, PAE arrays or logs. Each folding module's
    // flat predictions/ publishDir selects on this, and so does any consumer that
    // needs the structures back out of a predictions channel (see REFOLD).
    //
    // AF2 is the exception that takes options: which of its structures counts
    // depends on --af2_no_relax and --af2_keep_models, so the caller passes those
    // in rather than this class reading params.
    static boolean isModelStructure(String tool, String basename, Map opts = [:]) {
        def bn = basename.toString().replaceFirst(/^.*\//, '')
        switch (tool) {
            case 'af2':
            case 'af2_mono':
                if (opts.no_relax) {
                    return opts.keep_best ? (bn == 'ranked_0.cif') : (bn ==~ /unrelaxed_model_.*\.cif/)
                }
                return bn ==~ /relaxed_model_.*\.cif/
            case 'boltz':
                return bn ==~ /.*_model_\d+\.(cif|pdb)/
            case 'rf3':
                return bn ==~ /.*_sample-\d+_model\.cif/
            case 'af3':
                return bn ==~ /.*_seed-\d+_sample-\d+_model\.cif/
            case 'protenix':
                return bn ==~ /.*_sample_\d+\.cif/
            // OpenFold3, Chai-1 and both ESMFold2 variants share this naming.
            case 'openfold3':
            case 'chai1':
            case 'esmfold2':
            case 'esmfold2_fast':
                return bn ==~ /.*_seed_\d+_sample_\d+_model\.cif/
            default:
                throw new IllegalArgumentException("FoldNaming.isModelStructure: unknown tool '${tool}'")
        }
    }

    // The model structures in one task's predictions output, in a stable order.
    static List modelStructures(String tool, files, Map opts = [:]) {
        def list = (files instanceof List) ? files : [files]
        return list.findAll { isModelStructure(tool, it.name, opts) }.sort { it.name }
    }
}
