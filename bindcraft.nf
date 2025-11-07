#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIM_TO_CONTIGS } from './modules/trim_to_contigs.nf'
include { BINDCRAFT_CREATE_SETTINGS } from './modules/bindcraft_create_settings.nf'
include { BINDCRAFT } from './modules/bindcraft'
include { BINDCRAFT_REPORTING } from './modules/bindcraft_reporting.nf'

params.outdir = 'results'
params.bindcraft_advanced_settings_preset = 'default_4stage_multimer'
params.bindcraft_filters_preset = 'default_filters'
// params.bindcraft_advanced_settings_preset = 'default_4stage_multimer_flexible_hardtarget'
params.design_name = 'bindcraft_design'
params.input_pdb = false
params.hotspot_res = false
params.hotspot_subsample = 1.0
params.target_chains = false
params.binder_length_range = '60-150'
params.contigs = false
params.bindcraft_n_traj = 10
params.bindcraft_batch_size = 1
params.bindcraft_compress_html = true
params.bindcraft_compress_pdb = true
params.require_gpu = true
params.gpu_devices = ''
params.gpu_allocation_detect_process_regex = '(python.*/app/dl_binder_design/af2_initial_guess/predict\\.py|python.*/app/BindCraft/bindcraft\\.py|boltz predict|python.*/app/RFdiffusion/scripts/run_inference\\.py)'

// Function to validate hotspot_res parameter format
def validateHotspotRes(hotspot_res) {
    if (hotspot_res == null) {
        return true  // No validation needed if parameter is not provided
    }

    // Check for empty or whitespace-only string
    if (hotspot_res.trim().isEmpty()) {
        error "Invalid hotspot_res format: '${hotspot_res}'. Hotspot residues cannot be empty or whitespace only."
    }

    // Split by comma and trim whitespace
    def residues = hotspot_res.split(',').collect { it.trim() }

    // Check for empty elements (e.g., "A12," or ",A12" or "A12,,B15")
    if (residues.any { it.isEmpty() }) {
        error "Invalid hotspot_res format: '${hotspot_res}'. Hotspots must be specified like: 'A12,A15,A99' (comma separated {chainID}{resnum}, no empty elements)."
    }

    // Check each residue follows the format {chainID}{resnum}
    def pattern = ~/^[A-Z]\d+$/

    for (def residue : residues) {
        if (!(residue =~ pattern)) {
            error "Invalid hotspot residue format: '${residue}'. Hotspots must be specified like: 'A12,A15,A99' (comma separated {chainID}{resnum})."
        }
    }

    return residues
}

// Function to validate binder_length_range parameter format
def validateBinderLengthRange(binder_length_range) {
    if (binder_length_range == null) {
        return true  // No validation needed if parameter is not provided
    }

    // Check for empty or whitespace-only string
    if (binder_length_range.trim().isEmpty()) {
        error "Invalid binder_length_range format: '${binder_length_range}'. Binder length range cannot be empty or whitespace only."
    }

    // Check format: two integers separated by dash
    def pattern = ~/^\d+-\d+$/
    if (!(binder_length_range =~ pattern)) {
        error "Invalid binder_length_range format: '${binder_length_range}'. Binder length range must be specified like: '60-150' (two integers separated by dash)."
    }

    // Split by dash and convert to integers
    def parts = binder_length_range.split('-')
    def min_length = Integer.parseInt(parts[0])
    def max_length = Integer.parseInt(parts[1])

    // Check that first value is less than or equal to second value
    if (min_length > max_length) {
        error "Invalid binder_length_range: '${binder_length_range}'. First value (${min_length}) must be less than or equal to second value (${max_length})."
    }

    return true
}

// Function to validate target_chains parameter format
def validateTargetChains(target_chains) {
    if (target_chains == null) {
        return true  // No validation needed if parameter is not provided
    }

    // Check for empty or whitespace-only string
    if (target_chains.trim().isEmpty()) {
        error "Invalid target_chains format: '${target_chains}'. Target chains cannot be empty or whitespace only."
    }

    // Split by comma and trim whitespace
    def chains = target_chains.split(',').collect { it.trim() }

    // Check for empty elements (e.g., "A," or ",A" or "A,,B")
    if (chains.any { it.isEmpty() }) {
        error "Invalid target_chains format: '${target_chains}'. Target chains must be specified like: 'A' or 'A,C' (comma separated single capital letters, no empty elements)."
    }

    // Check each chain ID is a single capital letter
    def pattern = ~/^[A-Z]$/

    for (def chain : chains) {
        if (!(chain =~ pattern)) {
            error "Invalid target chain format: '${chain}'. Target chains must be specified like: 'A' or 'A,C' (comma separated single capital letters)."
        }
    }

    return chains
}

// Function to extract chain IDs from contigs string
def extractChainsFromContigs(contigs) {
    if (!contigs) {
        return []
    }

    // Find all capital letters in the contigs string
    def chainPattern = ~/[A-Z]/
    def matches = contigs.findAll(chainPattern)
    return matches.unique()
}

// Function to validate chain consistency between hotspot_res and target_chains
def validateChainConsistency(hotspot_residues, target_chains) {
    if (hotspot_residues == true || target_chains == true) {
        return true  // No validation needed if either parameter is not provided
    }

    // Extract chain IDs from hotspot residues
    def hotspot_chains = hotspot_residues.collect { residue ->
        residue[0]  // First character is the chain ID
    }.unique()

    // Check if all hotspot chains are in target chains
    def missing_chains = hotspot_chains.findAll { chain ->
        !target_chains.contains(chain)
    }

    if (missing_chains.size() > 0) {
        error "Chain consistency validation failed. Hotspot residues reference chains: ${hotspot_chains}, target chains specified: ${target_chains}, missing chains in target_chains: ${missing_chains}. All chains used in hotspot_res must be listed in target_chains."
    }

    return true
}

def paramsToMap(params) {
    def map = [:]
    params.each { key, value ->
        if (value instanceof Path || value instanceof File) {
            map[key] = value.toString()
        } else if (!(value instanceof Closure) && !(key in [
            'class', 'launchDir', 'projectDir', 'workDir'])) {
            map[key] = value
            }
    }
    return map
}

if (!params.input_pdb) {
    log.info"""
    ==================================================================
    🧬 BINDCRAFT NEXTFLOW WRAPPER 🧬
    ==================================================================

    Required arguments:
        --input_pdb                        The input PDB file

    Optional arguments:
        --outdir                  Output directory [default: ${params.outdir}]
        --design_name             Name of the design, used for output file prefixes [default: ${params.design_name}]
        --target_chains           Target chain(s) for binder design, "A" or "A,B" (mutually exclusive with --contigs) [default: ${params.target_chains}]
        --hotspot_res             Hotspot residues, eg "A473,A995,A411,A421" - you must include the chain ID in every hotspot
        --contigs                 Contigs to trim input PDB to, eg "[F2-23/F84-175/F205-267/0 G91-171/G209-263/0]" (mutually exclusive with --target_chains, automatically extracts target chains)
        --binder_length_range     Dash-separated min and max length for binders [default: ${params.binder_length_range}]
        --hotspot_subsample       Fraction of hotspot residues to randomly subsample (0.0-1.0) [default: ${params.hotspot_subsample}]
        --bindcraft_n_traj        Total number of designs attempts (trajectories) to generate [default: ${params.bindcraft_n_traj}]
        --bindcraft_batch_size    Number of designs to generate per batch [default: ${params.bindcraft_batch_size}]
        --bindcraft_advanced_settings_preset
                                  Preset for advanced settings [default: ${params.bindcraft_advanced_settings_preset}]
        --bindcraft_filters_preset
                                  Preset for filters [default: ${params.bindcraft_filters_preset}]
        --bindcraft_compress_html
                                  Compress batch output *.html) file with gzip [default: ${params.bindcraft_compress_html}]
        --bindcraft_compress_pdb
                                  Compress batch output *.pdb with gzip [default: ${params.bindcraft_compress_pdb}]

        --require_gpu           Fail tasks that go too slow without a GPU if no GPU is detected [default: ${params.require_gpu}]
        --gpu_devices           GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
        --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]

    """.stripIndent()
    exit 1
}

workflow {
    // Validate parameters
    def hotspot_residues = validateHotspotRes(params.hotspot_res)
    validateBinderLengthRange(params.binder_length_range)

    // Handle mutual exclusivity of contigs and target_chains
    def target_chains
    if (params.contigs && params.target_chains) {
        error 'Cannot specify both --contigs and --target_chains. Use --contigs to automatically extract target chains, or use --target_chains to specify them manually.'
    } else if (params.contigs) {
        // Extract target chains from contigs
        def extracted_chains = extractChainsFromContigs(params.contigs)
        if (extracted_chains.size() == 0) {
            error "No chain IDs found in contigs string: '${params.contigs}'"
        } else {
            target_chains = extracted_chains.join(',')
            log.info "Extracted target chains from contigs: ${target_chains}"
        }
    } else if (params.target_chains) {
        // Use manually specified target_chains
        target_chains = validateTargetChains(params.target_chains)
        if (target_chains != true) {
            target_chains = target_chains.join(',')
        }
    } else {
        // Neither contigs nor target_chains specified
        error 'Must specify either --contigs or --target_chains. Use --contigs to automatically extract target chains from contigs string, or use --target_chains to manually specify target chains (e.g., "A" or "A,B").'
    }

    // Check chain consistency if both parameters are valid
    if (hotspot_residues != true && target_chains != true) {
        validateChainConsistency(hotspot_residues, target_chains)
    }

    // TODO: Allow multiple --contigs definitions (comma-separated?), treated as alternative targets

    // TODO: Consider allowing multiple input PDBs and change the batching structure to accommodate this.
    //       Probably in this case we would disallow use of contigs (for simplicity) and require pre-trimmed PDBs,
    //       with a common --target_chains setting
    // (we don't want a complex configuration requiring mapping contig strings to input PDBs etc)

    // TODO: Consider 'samplesheet' support for complex combinations of targets, hotspots etc.

    ch_input_pdb = Channel.fromPath(params.input_pdb).first()
    def design_indices = 0..(params.bindcraft_n_traj - 1)
    def batches = design_indices.collate(params.bindcraft_batch_size)

    // Create batch info channel
    ch_batch_info = Channel.from(batches.withIndex())
        .map { batch, index -> [index, batch.size()] }

    if (params.contigs) {
        TRIM_TO_CONTIGS(
            ch_input_pdb,
            params.contigs
        )
        ch_input_pdb = TRIM_TO_CONTIGS.out.pdb
    }

    BINDCRAFT_CREATE_SETTINGS(
        ch_batch_info,
        ch_input_pdb,
        params.hotspot_res,
        target_chains,
        params.binder_length_range,
        params.design_name,
        params.hotspot_subsample
    )

    BINDCRAFT(
        ch_input_pdb,
        BINDCRAFT_CREATE_SETTINGS.out.settings_json,
        params.bindcraft_advanced_settings_preset,
        params.bindcraft_filters_preset,
        BINDCRAFT_CREATE_SETTINGS.out.batch_id,
        params.bindcraft_compress_html,
        params.bindcraft_compress_pdb
    )

    // TODO: Sort by Average_i_pTM
    
    // Merge CSV outputs from each batch into master files
    ch_final_stats_merged = BINDCRAFT.out.final_stats_csv
        .collectFile(name: 'final_design_stats.csv',
                     storeDir: "${params.outdir}/bindcraft",
                     keepHeader: true,
                     skip: 1)    

    ch_trajectory_stats_merged = BINDCRAFT.out.trajectory_stats_csv
        .collectFile(name: 'trajectory_stats.csv',
                     storeDir: "${params.outdir}/bindcraft",
                     keepHeader: true,
                     skip: 1)

    ch_mpnn_design_stats_merged = BINDCRAFT.out.mpnn_design_stats_csv
        .collectFile(name: 'mpnn_design_stats.csv',
                     storeDir: "${params.outdir}/bindcraft",
                     keepHeader: true,
                     skip: 1)

    // Collect and sum failure_csv rows - each file has a header and single data row with numeric values to sum
    ch_failure_csv_merged = BINDCRAFT.out.failure_csv
        .collect()
        .map { files ->
            def header = ''
            def sums = [:]
            def columnOrder = []

            files.eachWithIndex { file, index ->
                def lines = file.readLines()
                if (lines.size() >= 2) {
                    if (index == 0) {
                        // Get header from first file
                        header = lines[0]
                        columnOrder = header.split(',')
                    // Initialize sums map
                    def values = lines[1].split(',')
                    columnOrder.eachWithIndex { col, i ->
                        def value = values[i]
                        try {
                            sums[col] = Integer.parseInt(value)
                             } catch (NumberFormatException e) {
                            sums[col] = value
                        }
                    }
                    } else {
                    // Sum values from subsequent files
                    def values = lines[1].split(',')
                    columnOrder.eachWithIndex { col, i ->
                        def value = values[i]
                        if (sums[col] instanceof Number) {
                            try {
                                sums[col] += Integer.parseInt(value)
                                 } catch (NumberFormatException e) {
                            // Skip non-numeric values
                            }
                        }
                    }
                    }
                }
            }

            // Create summed CSV content
            def summedValues = columnOrder.collect { col ->
                sums[col] instanceof Number ? sums[col].toString() : sums[col]
            }.join(',')
            return "${header}\n${summedValues}"
        }
        .collectFile(name: 'failure_csv.csv', storeDir: "${params.outdir}/bindcraft")

    // Collect per-batch directories into a single list for reporting
    ch_batch_dirs_list = BINDCRAFT.out.batch_dir.collect()

    //ch_batch_dirs_list.view()

    // Generate BindCraft report
    BINDCRAFT_REPORTING(
        ch_batch_dirs_list,
        ch_failure_csv_merged,
        ch_final_stats_merged,
        ch_mpnn_design_stats_merged,
        ch_trajectory_stats_merged
    )

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
            success: workflow.success
        ]

        def output_file = "${params.outdir}/params.json"
        def json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params_json))

        new File(params.outdir).mkdirs()
        new File(output_file).text = json_string

        log.info "Pipeline parameters saved to: ${output_file}"
    }
}
