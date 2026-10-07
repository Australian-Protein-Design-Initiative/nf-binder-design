/*
Shared fold / fold_pulldown parameter validation. Returns [errors, warnings]
lists so callers can raise them via error()/log.warn() (Groovy classes cannot).
*/
class FoldValidation {
    // 'af2'      - AlphaFold2-multimer.
    // 'af2_mono' - AF2 monomer weights on a concatenated complex, chains separated
    //              only by an --af2_chain_break_offset jump in residue_index.
    static final List VALID_METHODS = ['af2', 'af2_mono', 'boltz', 'rf3', 'protenix', 'af3', 'openfold3', 'esmfold2', 'esmfold2_fast']
    static final List VALID_MSA_METHODS = ['jackhmmer_af2', 'mmseqs2_colabfold']
    // Engines that use --templates; the others fold without them.
    static final List TEMPLATE_METHODS = ['af3', 'boltz', 'openfold3']

    static List parseMethods(methodsParam) {
        return methodsParam.toString().split(',').collect { it.trim().toLowerCase() }
    }

    /**
     * Validate shared fold-engine params.
     * @param params Nextflow params map
     * @param methods parsed method list
     * @param opts optional map: checkColabfoldDbs (default true), hasMultimer (Boolean or null),
     *             af2DbPath (Path/File/String or null for uniprot/ check when hasMultimer+af2),
     *             pulldown (default false): fold_pulldown pairs nothing across chains (the
     *             binder is query-only) and builds AF2's multimer MSAs itself, so the
     *             ColabFold-multimer error/warning below do not apply
     * @return [errors, warnings] as Lists of Strings
     */
    static List validate(params, List methods, Map opts = [:]) {
        def errors = []
        def warnings = []
        def checkColabfoldDbs = opts.containsKey('checkColabfoldDbs') ? opts.checkColabfoldDbs : true
        def hasMultimer = opts.containsKey('hasMultimer') ? opts.hasMultimer : null
        def pulldown = opts.containsKey('pulldown') ? opts.pulldown : false

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

        ['boltz_batch_size', 'rf3_batch_size', 'protenix_batch_size', 'af3_batch_size', 'openfold3_batch_size',
         'esmfold2_batch_size'].each { pname ->
            def v = params[pname]
            if (v != null && !(v instanceof Boolean)) {
                def n = v as int
                if (n < 1) {
                    errors << "--${pname} must be >= 1 (got '${v}')"
                }
            }
        }
        if ('af3' in methods) {
            errors.addAll(af3WeightsErrors(params))
            if (!(params.af3_flash_attention in ['auto', 'triton', 'cudnn', 'xla'])) {
                errors << "--af3_flash_attention must be auto, triton, cudnn or xla (got '${params.af3_flash_attention}')"
            }
            if (params.af3_seeds && params.af3_seeds.toString().contains(',')) {
                warnings << (
                    "--af3_seeds takes a single base seed; only '${params.af3_seeds.toString().split(',')[0].trim()}' " +
                    "is used (batches use seed, seed+1, ...)."
                )
            }
        }
        // Protenix names its structures without the seed, so several seeds in one job
        // overwrite each other in fold/predictions/ and mis-pair confidence files.
        if ('protenix' in methods && params.protenix_seeds && params.protenix_seeds.toString().contains(',')) {
            errors << (
                "--protenix_seeds takes a single seed (got '${params.protenix_seeds}'). For more " +
                "structures use --n_predictions / --protenix_batch_size; batches use seed, seed+1, ..."
            )
        }
        if ('openfold3' in methods && params.openfold3_seeds && params.openfold3_seeds.toString().contains(',')) {
            warnings << (
                "--openfold3_seeds takes a single base seed; only '${params.openfold3_seeds.toString().split(',')[0].trim()}' " +
                "is used (batches use seed, seed+1, ...)."
            )
        }
        if (('esmfold2' in methods) || ('esmfold2_fast' in methods)) {
            if (params.esmfold2_seeds && params.esmfold2_seeds.toString().contains(',')) {
                warnings << (
                    "--esmfold2_seeds takes a single base seed; only '${params.esmfold2_seeds.toString().split(',')[0].trim()}' " +
                    "is used (batches use seed, seed+1, ...)."
                )
            }
            if (!(params.esmfold2_kernel_backend in ['fused', 'cuequivariance', 'none'])) {
                errors << (
                    "--esmfold2_kernel_backend must be fused, cuequivariance or none " +
                    "(got '${params.esmfold2_kernel_backend}')"
                )
            }
        }
        if ('esmfold2_fast' in methods) {
            warnings << (
                "ESMFold2-Fast (--methods esmfold2_fast) has no MSA encoder, so it folds from sequence " +
                "alone and does not use the MSA. Use --methods esmfold2 for MSA-conditioned ESMFold2."
            )
        }

        if (params.templates) {
            def untemplated = methods.findAll { !(it in TEMPLATE_METHODS) }
            if (untemplated) {
                warnings << (
                    "--templates is only used by ${TEMPLATE_METHODS.join(', ')}; " +
                    "${untemplated.join(', ')} will fold without templates."
                )
            }
            ['template_min_identity', 'template_min_coverage'].each { p ->
                def v = params[p] as double
                if (v < 0 || v > 1) {
                    errors << "--${p} must be between 0 and 1 (got '${params[p]}')"
                }
            }
            if ((params.template_max_per_chain as int) < 1) {
                errors << "--template_max_per_chain must be >= 1 (got '${params.template_max_per_chain}')"
            }
        }
        else if (pulldown && params.binder_templates) {
            warnings << "--binder_templates has no effect without --templates."
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

        // af2_mono uses the monomer weights and never pairs, so ColabFold MSAs suit it
        // and the jackhmmer requirement below is 'af2' only. Under jackhmmer_af2 though,
        // af2_mono's per-chain MSAs come from AF2's multimer search, which needs uniprot/.
        if (hasMultimer == true) {
            if (MsaSubsample.isEnabled(params.msa_subsample)) {
                errors << "--msa_subsample is not supported for multimer inputs (monomer only)."
            }
            if (!pulldown && params.msa_method == 'mmseqs2_colabfold' && !params.use_msa_server) {
                warnings << (
                    "multimer input with --msa_method mmseqs2_colabfold - ColabFold a3m " +
                    "headers carry no taxonomy, so RF3/Protenix/Boltz/AF3/OpenFold3/ESMFold2 will run UNPAIRED. Use " +
                    "--msa_method jackhmmer_af2, or --use_msa_server true (Boltz fetches + pairs itself)."
                )
            }
            if (!pulldown && 'af2' in methods && params.msa_method != 'jackhmmer_af2') {
                errors << (
                    "AF2 multimer requires --msa_method jackhmmer_af2 (AF2's native " +
                    "multimer MSA pipeline). Drop af2 from --methods for ColabFold multimer runs."
                )
            }
            def af2_multimer_search = ('af2' in methods) \
                || (!pulldown && 'af2_mono' in methods && params.msa_method == 'jackhmmer_af2')
            if (af2_multimer_search && opts.af2DbPath) {
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
     * AF3 weights are user-supplied (never in the container), so check the model
     * dir up front rather than failing inside a GPU job. AF3's loader rejects a
     * directory with more than one matching model file.
     */
    static List af3WeightsErrors(params) {
        def help = (
            "AlphaFold3 weights are not distributed with the pipeline. After reading the terms " +
            "(https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md), run " +
            "models/download_af3_weights.sh, or point --af3_model_dir at a directory containing " +
            "af3.bin.zst (or af3.bin)."
        )
        def dir_param = params.af3_model_dir
        if (!dir_param) {
            return ["--methods af3 needs --af3_model_dir. ${help}"]
        }
        def dir_str = dir_param.toString()
        if (dir_str ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/ && !dir_str.startsWith('file://')) {
            return []
        }
        def dir = new File(dir_str.replaceFirst(/^file:\/\//, ''))
        if (!dir.isDirectory()) {
            return ["--af3_model_dir '${dir_str}' does not exist or is not a directory. ${help}"]
        }
        def models = dir.listFiles().findAll { f ->
            f.isFile() && (f.name ==~ /.*\.bin(\.zst)?$/ || f.name ==~ /.*\.bin\.zst\.\d+$/)
        }
        if (!models) {
            return ["--af3_model_dir '${dir_str}' contains no AlphaFold3 weights (*.bin or *.bin.zst). ${help}"]
        }
        if (models.size() > 1) {
            return [
                "--af3_model_dir '${dir_str}' contains ${models.size()} model files " +
                "(${models*.name.join(', ')}); AlphaFold3 requires exactly one."
            ]
        }
        return []
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
