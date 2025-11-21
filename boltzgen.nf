#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOLTZGEN_DESIGN } from './modules/boltzgen/boltzgen_design'
include { BOLTZGEN_INVERSE_FOLDING } from './modules/boltzgen/boltzgen_inverse_folding'
include { BOLTZGEN_FOLDING } from './modules/boltzgen/boltzgen_folding'
include { BOLTZGEN_DESIGN_FOLDING } from './modules/boltzgen/boltzgen_design_folding'
include { BOLTZGEN_AFFINITY } from './modules/boltzgen/boltzgen_affinity'
include { BOLTZGEN_MERGE } from './modules/boltzgen/boltzgen_merge'
include { BOLTZGEN_ANALYSIS } from './modules/boltzgen/boltzgen_analysis'
include { BOLTZGEN_FILTERING } from './modules/boltzgen/boltzgen_filtering'

// Default parameters
params.config_yaml = false
params.outdir = 'results'
params.design_name = false
params.protocol = 'protein-anything'
params.num_designs = 100
params.batch_size = 10
params.budget = 10
params.devices = false
params.num_workers = false
params.inverse_fold_num_sequences = false

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
    if (params.config_yaml == false) {
        log.info(
            """
        ==================================================================
        🧬 BOLTZGEN PIPELINE 🧬
        ==================================================================

        Required arguments:
            --config_yaml              Path to BoltzGen YAML config file

        Optional arguments:
            --outdir                      Output directory [default: ${params.outdir}]
            --design_name                 Name of the design, used for output file prefixes [default: auto-derived from config_yaml basename]
            --protocol                    Protocol type (protein-anything, peptide-anything, protein-small_molecule, nanobody-anything) [default: ${params.protocol}]
            --num_designs                 Total number of designs to generate [default: ${params.num_designs}]
            --batch_size                  Number of designs per batch [default: ${params.batch_size}]
            --budget                      Final diversity-optimized set size [default: ${params.budget}]
            --devices                     Number of GPU devices [default: unspecified]
            --num_workers                 Number of DataLoader workers [default: unspecified]
            --inverse_fold_num_sequences  Number of sequences per design in inverse folding step [default: unspecified]

        """.stripIndent()
        )
        exit(1)
    }

    def design_name = params.design_name
    
    // Set design_name from config_yaml basename if not explicitly set
    if (!params.design_name) {
        def config_file = new File(params.config_yaml)
        def config_basename = config_file.getName()
        // Remove .yaml or .yml extension
        design_name = config_basename.replaceFirst(/\.(yaml|yml)$/, '')
    }

    // Validate config_yaml exists
    ch_config_yaml = Channel.fromPath(params.config_yaml).first()

    // Parse config.yaml to extract input files
    def config_file = file(params.config_yaml)
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

    // Generate batch start indices - create separate channels to avoid double consumption
    def batch_indices = (0..params.num_designs - 1).findAll { it % params.batch_size == 0 }
    
    ch_batch_n_designs = Channel.from(batch_indices)
        .map { start_idx ->
            Math.min(params.batch_size, params.num_designs - start_idx)
        }
    
    ch_batch_start_idx = Channel.from(batch_indices)

    // Create channels for constant values
    ch_design_name = Channel.value(design_name)
    ch_protocol = Channel.value(params.protocol)
    ch_devices = Channel.value(params.devices)
    ch_num_workers = Channel.value(params.num_workers)
    ch_inverse_fold_num_sequences = Channel.value(params.inverse_fold_num_sequences)

    // Phase 1: Design (Parallel)
    BOLTZGEN_DESIGN(
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        ch_batch_n_designs,
        ch_batch_start_idx,
        ch_devices,
        ch_num_workers,
    )

    // Phase 2: Inverse Folding (Parallel)
    // Extract start_index from batch_dir path (batch_dir is batch_{start_index}/)
    ch_batch_with_start = BOLTZGEN_DESIGN.out.batch_dir.map { batch_dir ->
        def dir_path = batch_dir.toString()
        def dir_name = new File(dir_path).name
        def start_idx = dir_name.replaceAll(/batch_/, '').toInteger()
        return [batch_dir, start_idx]
    }

    BOLTZGEN_INVERSE_FOLDING(
        ch_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        ch_batch_with_start.map { _batch_dir, start_idx -> start_idx },
        ch_devices,
        ch_num_workers,
        ch_inverse_fold_num_sequences,
    )

    // Phase 3: Folding (Parallel)
    ch_folding_batch_with_start = BOLTZGEN_INVERSE_FOLDING.out.batch_dir.map { batch_dir ->
        def dir_path = batch_dir.toString()
        def dir_name = new File(dir_path).name
        def start_idx = dir_name.replaceAll(/batch_/, '').toInteger()
        return [batch_dir, start_idx]
    }

    BOLTZGEN_FOLDING(
        ch_folding_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        ch_folding_batch_with_start.map { _batch_dir, start_idx -> start_idx },
        ch_devices,
        ch_num_workers,
    )

    // Phase 4: Design Folding (Parallel, if applicable)
    if (params.protocol in ['protein-anything', 'protein-small_molecule']) {
        ch_design_folding_batch_with_start = BOLTZGEN_FOLDING.out.batch_dir.map { batch_dir ->
            def dir_path = batch_dir.toString()
            def dir_name = new File(dir_path).name
            def start_idx = dir_name.replaceAll(/batch_/, '').toInteger()
            return [batch_dir, start_idx]
        }

        BOLTZGEN_DESIGN_FOLDING(
            ch_design_folding_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
            ch_config_yaml,
            ch_input_files,
            ch_design_name,
            ch_protocol,
            ch_design_folding_batch_with_start.map { _batch_dir, start_idx -> start_idx },
            ch_devices,
            ch_num_workers,
        )
        
        if (params.protocol == 'protein-small_molecule') {
            ch_affinity_batch_with_start = BOLTZGEN_DESIGN_FOLDING.out.batch_dir.map { batch_dir ->
                def dir_path = batch_dir.toString()
                def dir_name = new File(dir_path).name
                def start_idx = dir_name.replaceAll(/batch_/, '').toInteger()
                return [batch_dir, start_idx]
            }

            BOLTZGEN_AFFINITY(
                ch_affinity_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
                ch_config_yaml,
                ch_input_files,
                ch_design_name,
                ch_protocol,
                ch_affinity_batch_with_start.map { _batch_dir, start_idx -> start_idx },
                ch_devices,
                ch_num_workers,
            )
            ch_design_folded = BOLTZGEN_AFFINITY.out.batch_dir
        }
        else {
            ch_design_folded = BOLTZGEN_DESIGN_FOLDING.out.batch_dir
        }
    }
    else {
        ch_design_folded = BOLTZGEN_FOLDING.out.batch_dir
    }

    // Phase 5: Merge
    BOLTZGEN_MERGE(
        ch_design_folded.collect()
    )

    // Phase 6: Analysis (Single process)
    ch_merged_dir = BOLTZGEN_MERGE.out.merged_dir.first()

    BOLTZGEN_ANALYSIS(
        ch_merged_dir,
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
    )

    // Phase 7: Filtering (Single process)
    ch_analysis_merged_dir = BOLTZGEN_ANALYSIS.out.merged_dir.first()

    BOLTZGEN_FILTERING(
        ch_analysis_merged_dir,
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        Channel.value(params.budget),
    )

    ///////////////////////////////////////////////////////////////////////////
    workflow.onComplete = {
        // Write the pipeline parameters to a JSON file
        def params_json = [:]

        params_json['params'] = paramsToMap(params)

        params_json['workflow'] = [
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

        def output_file = "${params.outdir}/params.json"
        def json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params_json))

        new File(params.outdir).mkdirs()
        new File(output_file).text = json_string

        log.info("Pipeline parameters saved to: ${output_file}")
    }
}
