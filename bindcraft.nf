#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIM_TO_CONTIGS } from './modules/trim_to_contigs.nf'
include { BINDCRAFT_CREATE_SETTINGS } from './modules/bindcraft_create_settings.nf'
include { BINDCRAFT } from './modules/bindcraft'

params.outdir = 'results'
params.bindcraft_advanced_settings_preset = 'default_4stage_multimer'
// params.bindcraft_advanced_settings_preset = 'default_4stage_multimer_flexible_hardtarget'
params.design_name = 'bindcraft_design'
params.input_pdb = false
params.hotspot_res = false
params.target_chains = 'A'
params.binder_length_range = '60-150'
params.bindcraft_n_designs = 10
params.bindcraft_batch_size = 1
params.require_gpu = true
params.gpu_devices = ''
params.gpu_allocation_detect_process_regex = '(python.*/app/dl_binder_design/af2_initial_guess/predict\\.py|python.*/app/BindCraft/bindcraft\\.py|boltz predict|python.*/app/RFdiffusion/scripts/run_inference\\.py)'

if (!params.input_pdb || !params.hotspot_res) {
    log.info"""
    ==================================================================
    🧬 BINDCRAFT NEXTFLOW WRAPPER 🧬
    ==================================================================

    Required arguments:
        --input_pdb                        The input PDB file
        --hotspot_res                      Hotspot residues, eg "A473,A995,A411,A421"

    Optional arguments:
        --outdir                Output directory [default: ${params.outdir}]
        --design_name           Name of the design, used for output file prefixes [default: ${params.design_name}]
        --target_chains         Target chain(s) for binder design [default: ${params.target_chains}]
        --contigs               Contigs to trim input PDB to, eg "[F2-23/F84-175/F205-267/0 G91-171/G209-263/0]"
        --binder_length_range   Dash-separated min and max length for binders [default: ${params.binder_length_range}]
        --bindcraft_n_designs   Total number of designs to generate [default: ${params.bindcraft_n_designs}]
        --bindcraft_batch_size  Number of designs to generate per batch [default: ${params.bindcraft_batch_size}]
        --bindcraft_advanced_settings_preset
                                Preset for advanced settings [default: ${params.bindcraft_advanced_settings_preset}]

        --require_gpu           Fail tasks that go too slow without a GPU if no GPU is detected [default: ${params.require_gpu}]
        --gpu_devices           GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
        --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]

    """.stripIndent()
    exit 1
}

workflow {
    // TODO: In theory, if --contigs is specified then --target_chains can be inferred from the contigs
    //       Do that and make --target_chains and --contigs mutually exclusive
    // TODO: Allow multiple --contigs definitions (comma-separated?)
    // TODO: Allow a hotspot subsampling parameter (eg "--hotspot_res_subsample 0.5") - the settings or hotspots chosen need to
    //       be properly recorded along with the results

    // TODO: Consider allowing multiple input PDBs and change the batching structure to accommodate this.
    //       Probably in this case we would disallow use of contigs (for simiplicity) and require pre-trimmed PDBs,
    //       with a common --target_chains setting
    // (we don't want a complex configuration requiring mapping contig strings to input PDBs etc)

    ch_input_pdb = Channel.fromPath(params.input_pdb).first()
    def design_indices = 0..(params.bindcraft_n_designs - 1)
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
        params.design_name
    )

    BINDCRAFT(
        ch_input_pdb,
        BINDCRAFT_CREATE_SETTINGS.out.settings_json,
        params.bindcraft_advanced_settings_preset,
        BINDCRAFT_CREATE_SETTINGS.out.batch_id
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
        .collectFile(name: 'failure_stats_summed.csv', storeDir: "${params.outdir}/bindcraft")
}
