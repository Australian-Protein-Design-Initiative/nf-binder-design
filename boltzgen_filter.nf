#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters
params.help = false
params.config_yaml = false
params.run = 'results/boltzgen/merged'
params.outdir = 'results'
params.protocol = false
params.budget = 10
params.alpha = false
params.filter_biased = false
params.metrics_override = []
params.additional_filters = []
params.size_buckets = []
params.refolding_rmsd_threshold = false

include { BOLTZGEN_FILTERING } from './modules/boltzgen/boltzgen_filtering'

include { detectParams; buildFilteringArgs } from './modules/boltzgen/boltzgen_utils'

def paramsToMap(params) {
    def map = [:]
    params.each { key, value ->
        if (value instanceof Path || value instanceof File) {
            map[key] = value.toString()
        }
        else if (!(value instanceof Closure) && !(key in [
            'class',
            'launchDir',
            'projectDir',
            'workDir',
        ])) {
            map[key] = value
        }
    }
    return map
}

workflow {
    // Show help message
    if (params.help) {
        log.info(
            """
        ==================================================================
        🧬 BOLTZGEN FILTERING WORKFLOW 🧬
        ==================================================================

        Re-run filtering on existing BoltzGen merged results with custom parameters.

        If left unset, --config_yaml and --protocol are auto-detected from ${params.run}/../params.json.

        Optional arguments:
            --run                         Input directory containing merged results [default: ${params.run}]
            --outdir                      Output directory [default: ${params.outdir}]
            --config_yaml                 Path to BoltzGen YAML config file [auto-detected from --run/../params.json if not specified]
            --protocol                    Protocol type (protein-anything, peptide-anything, protein-small_molecule, nanobody-anything) [auto-detected from --run/../params.json if not specified]
            --budget                      Final diversity-optimized set size [default: ${params.budget}]
            --alpha                       Trade-off for sequence diversity selection: 0.0=quality-only, 1.0=diversity-only
            --filter_biased               Remove amino-acid composition outliers (default: true, use --filter_biased false to disable)
            --metrics_override            Per-metric inverse-importance weights for ranking. Format: metric_name=weight (e.g., 'plip_hbonds_refolded=4' 'delta_sasa_refolded=2')
            --additional_filters          Extra hard filters. Format: feature>threshold or feature<threshold (e.g., 'design_ALA>0.3' 'design_GLY<0.2')
            --size_buckets                Optional constraint for maximum number of designs in size ranges. Format: min-max:count (e.g., '10-20:5' '20-30:10')
            --refolding_rmsd_threshold     Threshold used for RMSD-based filters (lower is better)

        Example:
            nextflow run boltzgen_filter.nf --budget 20 --alpha 0.05

        """.stripIndent()
        )
        exit(0)
    }

    // Auto-detect config_yaml and protocol from params.json if not specified
    def loaded_params = detectParams(params)
    def config_yaml = loaded_params.config_yaml
    def protocol = loaded_params.protocol

    // Validate config_yaml exists
    ch_config_yaml = Channel.fromPath(config_yaml).first()

    // Parse config.yaml to extract input files
    def config_file = file(config_yaml)
    def config_dir = config_file.parent
    def parse_cmd = ["python3", "${projectDir}/bin/boltzgen/parse_boltzgen_config.py", config_file.toString(), "--config-dir", config_dir.toString()]
    def parse_output = parse_cmd.execute().text.trim()
    def input_file_paths = parse_output.split('\n').findAll { it.trim() }

    // Create channel of input files
    if (input_file_paths.isEmpty()) {
        ch_input_files = Channel.value([])
    } else {
        ch_input_files = Channel.from(input_file_paths)
            .map { file_path -> file(file_path.trim()) }
            .collect()
            .map { files -> files }
            .first()
    }

    // Set design_name from config_yaml basename
    def config_file_obj = new File(config_yaml)
    def config_basename = config_file_obj.getName()
    def design_name = config_basename.replaceFirst(/\.(yaml|yml)$/, '')

    // Build filtering arguments string
    def filtering_args = buildFilteringArgs(
        params.alpha,
        params.filter_biased,
        params.metrics_override,
        params.additional_filters,
        params.size_buckets,
        params.refolding_rmsd_threshold
    )

    // Create channels for constant values
    ch_merged_dir = Channel.fromPath("${params.run}", type: 'dir').first()
    ch_design_name = Channel.value(design_name)
    ch_protocol = Channel.value(protocol)
    ch_filtering_arg_list = Channel.value(filtering_args)

    // Run filtering
    BOLTZGEN_FILTERING(
        ch_merged_dir,
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        Channel.value(params.budget),
        ch_filtering_arg_list,
    )

    ///////////////////////////////////////////////////////////////////////////
    workflow.onComplete = {
        // Write the pipeline parameters to a JSON file
        def filter_params_json = [:]

        filter_params_json['params'] = paramsToMap(params)

        filter_params_json['workflow'] = [
            name: workflow.manifest.name,
            version: workflow.manifest.version,
            revision: workflow.revision ?: null,
            commit: workflow.commitId ?: null,
            runName: workflow.runName,
            start: workflow.start.format('yyyy-MM-dd HH:mm:ss'),
            complete: workflow.complete.format('yyyy-MM-dd HH:mm:ss'),
            duration: workflow.duration,
            success: workflow.success,
        ]

        def output_file = "${params.outdir}/boltzgen/filtered/final_ranked_designs/config/filter_params.json"
        def json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(filter_params_json))

        // Create the directory for the output file
        new File(output_file).parentFile.mkdirs()
        new File(output_file).text = json_string

        log.info("Filtering parameters saved to: ${output_file}")
    }
}

