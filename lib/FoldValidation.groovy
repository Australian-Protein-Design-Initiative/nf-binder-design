/*
Shared fold / fold_pulldown parameter validation. Returns [errors, warnings]
lists so callers can raise them via error()/log.warn() (Groovy classes cannot).
*/
class FoldValidation {
    // 'af2'      - AlphaFold2-multimer.
    // 'af2_mono' - AF2 monomer weights on a concatenated complex, chains separated
    //              only by an --af2_chain_break_offset jump in residue_index.
    static final List VALID_METHODS = ['af2', 'af2_mono', 'boltz', 'rf3', 'protenix']
    static final List VALID_MSA_METHODS = ['jackhmmer_af2', 'mmseqs2_colabfold']

    static List parseMethods(methodsParam) {
        return methodsParam.toString().split(',').collect { it.trim().toLowerCase() }
    }

    /**
     * Validate shared fold-engine params.
     * @param params Nextflow params map
     * @param methods parsed method list
     * @param opts optional map: checkColabfoldDbs (default true), hasMultimer (Boolean or null),
     *             af2DbPath (Path/File/String or null for uniprot/ check when hasMultimer+af2)
     * @return [errors, warnings] as Lists of Strings
     */
    static List validate(params, List methods, Map opts = [:]) {
        def errors = []
        def warnings = []
        def checkColabfoldDbs = opts.containsKey('checkColabfoldDbs') ? opts.checkColabfoldDbs : true
        def hasMultimer = opts.containsKey('hasMultimer') ? opts.hasMultimer : null

        methods.each { m ->
            if (!(m in VALID_METHODS)) {
                errors << "unknown --methods entry '${m}' (valid: ${VALID_METHODS.join(', ')})"
            }
        }
        if (!(params.msa_method in VALID_MSA_METHODS)) {
            errors << "unknown --msa_method '${params.msa_method}' (valid: ${VALID_MSA_METHODS.join(', ')})"
        }

        if ('af2_mono' in methods) {
            if ((params.af2_chain_break_offset as int) < 33) {
                errors << (
                    "--af2_chain_break_offset must be > 32 (got '${params.af2_chain_break_offset}'): " +
                    "AF2 clips relative positions at 32, so a smaller jump does not read " +
                    "as a chain break and the chains are modelled as covalently joined."
                )
            }
            if (params.af2_keep_models != 'best') {
                warnings << (
                    "--methods af2_mono with --af2_keep_models=${params.af2_keep_models}: without an " +
                    "initial guess the monomer models frequently fail to dock the chains at all, and " +
                    "only ranking separates the good pose from the failures (ranking_confidence for " +
                    "monomer presets is mean pLDDT). Prefer --af2_keep_models best, or filter on " +
                    "interface pLDDT before using af2_mono poses in any consensus."
                )
            }
        }

        if ('af2' in methods || 'af2_mono' in methods) {
            if (!(params.af2_keep_models in ['all', 'best'])) {
                errors << "--af2_keep_models must be 'all' or 'best' (got '${params.af2_keep_models}')"
            }
            if (params.n_predictions && params.af2_keep_models == 'all' && (params.n_predictions as int) % 5 != 0) {
                def n_runs = ((params.n_predictions as int) + 4).intdiv(5)
                warnings << (
                    "AF2 --af2_keep_models=all emits 5 models per run; --n_predictions " +
                    "${params.n_predictions} is not a multiple of 5, so AF2 will produce ${n_runs * 5} " +
                    "models across ${n_runs} runs (nearest multiple of 5 >= n_predictions)."
                )
            }
        }

        ['boltz_batch_size', 'rf3_batch_size', 'protenix_batch_size'].each { pname ->
            def v = params[pname]
            if (v != null && !(v instanceof Boolean)) {
                def n = v as int
                if (n < 1) {
                    errors << "--${pname} must be >= 1 (got '${v}')"
                }
            }
        }
        if (params.n_predictions && (params.n_predictions as int) < 1) {
            errors << "--n_predictions must be >= 1 (got '${params.n_predictions}')"
        }

        def engens_clustering = params.engens_clustering.toString().split(',').collect { it.trim().toLowerCase() }
        engens_clustering.each { c ->
            if (!(c in ['gmm', 'km', 'hdbscan'])) {
                errors << "unknown --engens_clustering entry '${c}' (valid: gmm, km, hdbscan)"
            }
        }
        if (!(params.engens_dimred.toString().toLowerCase() in ['umap'])) {
            errors << "--engens_dimred must be 'umap' (got '${params.engens_dimred}')"
        }
        if (!(params.engens_gmm_ic.toString().toLowerCase() in ['aic', 'bic'])) {
            errors << "--engens_gmm_ic must be 'aic' or 'bic' (got '${params.engens_gmm_ic}')"
        }
        if ((params.engens_min_structures as int) < 2) {
            errors << "--engens_min_structures must be >= 2 (got '${params.engens_min_structures}')"
        }
        if ((params.engens_max_clusters as int) < 2) {
            errors << "--engens_max_clusters must be >= 2 (got '${params.engens_max_clusters}')"
        }
        def engens_featurizers = params.engens_featurizers.toString().split(',').collect { it.trim().toLowerCase() }.findAll { it }
        if (!engens_featurizers) {
            errors << "--engens_featurizers must list at least one of: default, 3di, pb"
        }
        engens_featurizers.each { f ->
            if (!(f in ['default', '3di', 'pb'])) {
                errors << "unknown --engens_featurizers entry '${f}' (valid: default, 3di, pb)"
            }
        }

        if (MsaSubsample.isEnabled(params.msa_subsample)) {
            def raw = (params.msa_subsample instanceof Boolean \
                || params.msa_subsample.toString().trim().toLowerCase() in ['true', '1', 'yes']) \
                ? MsaSubsample.MSA_DEPTH_DEFAULTS \
                : params.msa_subsample.toString().trim()
            raw.split(',').each { part ->
                def pdepth = part.trim()
                if (!pdepth) {
                    return
                }
                if (!(pdepth ==~ /^\d+:\d+$/)) {
                    errors << (
                        "invalid --msa_subsample depth '${pdepth}' " +
                        "(expected max_seq:max_extra_seq, or true for CF-random defaults)"
                    )
                }
                else {
                    def bits = pdepth.split(':')
                    if ((bits[0] as int) < 1) {
                        errors << "--msa_subsample max_seq must be >= 1 (got '${pdepth}')"
                    }
                }
            }
        }

        if (checkColabfoldDbs && params.msa_method == 'mmseqs2_colabfold' \
            && !params.use_remote_server && !(params.uniref30 && params.colabfold_envdb)) {
            errors << (
                "--msa_method mmseqs2_colabfold needs --uniref30 and --colabfold_envdb " +
                "(local ColabFold DBs) or --use_remote_server true."
            )
        }

        // NB: the two af2 gates below are deliberately 'af2' only, not af2_mono.
        // af2_mono uses the monomer weights and never pairs, so it needs neither the
        // uniprot/ all-seqs DB nor the native jackhmmer multimer MSA pipeline -
        // ColabFold MSAs are exactly what it wants. Do not widen these to af2_mono.
        if (hasMultimer == true) {
            if (MsaSubsample.isEnabled(params.msa_subsample)) {
                errors << "--msa_subsample is not supported for multimer inputs (monomer only)."
            }
            if (params.msa_method == 'mmseqs2_colabfold' && !params.use_msa_server) {
                warnings << (
                    "multimer input with --msa_method mmseqs2_colabfold - ColabFold a3m " +
                    "headers carry no taxonomy, so RF3/Protenix/Boltz will run UNPAIRED. Use " +
                    "--msa_method jackhmmer_af2, or --use_msa_server true (Boltz fetches + pairs itself)."
                )
            }
            if ('af2' in methods && params.msa_method != 'jackhmmer_af2') {
                errors << (
                    "AF2 multimer requires --msa_method jackhmmer_af2 (AF2's native " +
                    "multimer MSA pipeline). Drop af2 from --methods for ColabFold multimer runs."
                )
            }
            if ('af2' in methods && opts.af2DbPath) {
                def uniprot_dir = new File("${opts.af2DbPath}/uniprot")
                if (!uniprot_dir.exists()) {
                    errors << (
                        "AF2 multimer needs a uniprot/ DB under --af2_db_path, absent in " +
                        "'${opts.af2DbPath}'. Point --af2_db_path at the 2021 snapshot (e.g. " +
                        "/mnt/datasets/alphafold/alphafold_20211129), which has uniprot/ + pdb_seqres/."
                    )
                }
            }
        }

        return [errors, warnings]
    }

    /**
     * Count FASTA records and detect empty sequences. Returns
     * [n_chains, empty_records]. Accepts File or Path.
     */
    static List countFastaChains(f) {
        def file = (f instanceof File) ? f : f.toFile()
        int n_chains = 0
        int empty_records = 0
        boolean in_record = false
        int seq_len = 0
        file.eachLine { line ->
            def l = line.trim()
            if (l.startsWith('>')) {
                if (in_record && seq_len == 0) {
                    empty_records += 1
                }
                n_chains += 1
                in_record = true
                seq_len = 0
            }
            else if (in_record) {
                seq_len += l.length()
            }
        }
        if (in_record && seq_len == 0) {
            empty_records += 1
        }
        return [n_chains, empty_records]
    }
}
