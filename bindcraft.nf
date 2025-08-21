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
params.target_chains = 'A'
params.binder_length_range = '60-150'
params.contigs = false
params.bindcraft_n_traj = 10
params.bindcraft_batch_size = 1
params.bindcraft_compress_html = true
params.bindcraft_compress_pdb = true
params.require_gpu = true
params.gpu_devices = ''
params.gpu_allocation_detect_process_regex = '(python.*/app/dl_binder_design/af2_initial_guess/predict\\.py|python.*/app/BindCraft/bindcraft\\.py|boltz predict|python.*/app/RFdiffusion/scripts/run_inference\\.py)'

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
        --target_chains           Target chain(s) for binder design [default: ${params.target_chains}]
        --hotspot_res             Hotspot residues, eg "A473,A995,A411,A421" - you must include the chain ID in every hotspot
        --contigs                 Contigs to trim input PDB to, eg "[F2-23/F84-175/F205-267/0 G91-171/G209-263/0]"
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
    // TODO: In theory, if --contigs is specified then --target_chains can be inferred from the contigs
    //       Do that and make --target_chains and --contigs mutually exclusive
    // TODO: Allow multiple --contigs definitions (comma-separated?), treated as alternative targets

    // TODO: Consider allowing multiple input PDBs and change the batching structure to accommodate this.
    //       Probably in this case we would disallow use of contigs (for simiplicity) and require pre-trimmed PDBs,
    //       with a common --target_chains setting
    // (we don't want a complex configuration requiring mapping contig strings to input PDBs etc)

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
        params.target_chains,
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

workflow.onComplete {
    // Write the pipeline parameters to a JSON file
    def params_json = [:]

    params_json['params'] = paramsToMap(params)

    params_json['workflow'] = [
        name: workflow.manifest.name,
        version: workflow.manifest.version,
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