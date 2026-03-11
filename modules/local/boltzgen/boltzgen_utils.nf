// Utility functions for BoltzGen workflows

process PARSE_BOLTZGEN_CONFIG {
    tag "parse_boltzgen_config"

    input:
    path config_yaml
    // config_dir is a path outside the `work` directory - while not ideal, 
    // this process is always intended to be executed locally or on a shared 
    // filesystem. We do this in a process{} rather than a bare Groovy 
    // cmd.execute() since we need pyyaml as a dependency to run 
    // parse_boltzgen_config.py, and cannot rely on this being available 
    // via the system (non-containerized) Python installation.
    // This script needs to capture the real (non-`work/`) paths of 
    // input files relative to the real path of config_yaml,
    // so we cannot remain sandboxed within the `work` directory for this step.
    val  config_dir

    output:
    stdout

    script:
    """
    python3 ${projectDir}/bin/boltzgen/parse_boltzgen_config.py ${config_yaml} --config-dir ${config_dir}
    """
}

def detectParams(params) {
    // If config_yaml and protocol are already set, return them
    def config_yaml = params.config_yaml
    def protocol = params.protocol
    
    // Validate run directory exists (needed for reading params.json)
    def run_dir = file(params.run)
    if (!run_dir.exists() || !run_dir.isDirectory()) {
        error "Input directory does not exist or is not a directory: ${params.run}"
    }

    // Try to read params.json for auto-detection
    def params_json_path = file("${params.run}/../../params.json")
    def params_json = null
    if (params_json_path.exists()) {
        params_json = new groovy.json.JsonSlurper().parse(params_json_path)
    }

    // Determine config_yaml - use existing value if set, otherwise try to read from params.json
    if (!config_yaml || config_yaml == false) {
        if (params_json && params_json.params && params_json.params.config_yaml) {
            config_yaml = params_json.params.config_yaml
            log.info("Auto-detected config_yaml from params.json: ${config_yaml}")
        }
    }
    if (!config_yaml || config_yaml == false) {
        error "config_yaml not specified and could not be auto-detected from ${params_json_path}. Please specify --config_yaml."
    }

    // Determine protocol - use existing value if set, otherwise try to read from params.json
    if (!protocol || protocol == false) {
        if (params_json && params_json.params && params_json.params.protocol) {
            protocol = params_json.params.protocol
            log.info("Auto-detected protocol from params.json: ${protocol}")
        }
    }
    if (!protocol || protocol == false) {
        error "Protocol not specified and could not be auto-detected from ${params_json_path}. Please specify --protocol."
    }

    return [config_yaml: config_yaml, protocol: protocol]
}

def buildFilteringArgs(alpha, filter_biased, metrics_override, additional_filters, size_buckets, refolding_rmsd_threshold) {
    def filtering_args_list = []
    
    if (alpha != false) {
        filtering_args_list.add("--alpha ${alpha}")
    }
    
    // filter_biased defaults to true in boltzgen, only pass if user wants to disable it
    if (filter_biased != false) {
        if (filter_biased.toString().toLowerCase() == 'false') {
            filtering_args_list.add("--filter_biased false")
        }
        // If true, don't pass it (use boltzgen default)
    }
    
    // Handle metrics_override - can be a list or a single string
    // Quote values to handle special characters like = in shell
    if (metrics_override) {
        def metrics_list = metrics_override instanceof List ? metrics_override : [metrics_override]
        def filtered_metrics = metrics_list.findAll { it && it.toString().trim() }
        if (filtered_metrics.size() > 0) {
            def metrics_override_str = filtered_metrics.collect { "'${it}'" }.join(' ')
            filtering_args_list.add("--metrics_override ${metrics_override_str}")
        }
    }
    
    // Handle additional_filters - can be a list or a single string
    // Quote values to handle special characters like < and > in shell
    if (additional_filters) {
        def filters_list = additional_filters instanceof List ? additional_filters : [additional_filters]
        def filtered_filters = filters_list.findAll { it && it.toString().trim() }
        if (filtered_filters.size() > 0) {
            def additional_filters_str = filtered_filters.collect { "'${it}'" }.join(' ')
            filtering_args_list.add("--additional_filters ${additional_filters_str}")
        }
    }
    
    // Handle size_buckets - can be a list or a single string
    // Quote values to handle special characters like - and : in shell
    if (size_buckets) {
        def buckets_list = size_buckets instanceof List ? size_buckets : [size_buckets]
        def filtered_buckets = buckets_list.findAll { it && it.toString().trim() }
        if (filtered_buckets.size() > 0) {
            def size_buckets_str = filtered_buckets.collect { "'${it}'" }.join(' ')
            filtering_args_list.add("--size_buckets ${size_buckets_str}")
        }
    }
    
    if (refolding_rmsd_threshold != false) {
        filtering_args_list.add("--refolding_rmsd_threshold ${refolding_rmsd_threshold}")
    }
    
    return filtering_args_list
}

