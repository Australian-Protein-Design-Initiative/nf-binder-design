#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
RFDiffusion3-based binder design workflow

Usage via main.nf:
  nextflow run main.nf --method rfd3 --rfd3_config config.json
  nextflow run main.nf --method rfd3 --input_pdb target.pdb --contigs "A17-131,/0,50-120" --hotspot_res "A56,A115,A123"

*/

// Default parameters
params.input_pdb = false
params.outdir = 'results'

params.design_name = 'rfd3'
params.contigs = ''
params.hotspot_res = false

// RFDiffusion3 params
params.rfd3_config = false
params.rfd3_n_designs = 1
params.rfd3_batch_size = 1
params.rfd3_step_scale = 3
params.rfd3_gamma_0 = 0.2
params.rfd3_allow_loopy = false
params.rfd3_extra_args = ''

params.rfd3_filters = false

// RosettaFold3 params
//params.rf3_ckpt_path = '/weights/rf3_foundry_01_24_latest_remapped.ckpt'
params.rf3_ckpt_path = '/models/foundry/rf3_foundry_01_24_latest_remapped.ckpt'

// MPNN params - legacy (pmpnn_*) for compatibility with rfd workflow
params.pmpnn_seqs_per_struct = 1
params.pmpnn_temperature = false
params.pmpnn_augment_eps = false
params.pmpnn_omit_aas = false
params.pmpnn_weights = false

// Boltz refolding params
params.refold_with = '' // comma-separated list of methods, e.g., 'boltz'
params.refold_max = false
params.refold_filter_sort = 'pair_pae_min'
params.refold_use_msa_server = false
params.refold_create_target_msa = false
params.refold_target_templates = false
params.refold_target_fasta = false
params.uniref30 = false
params.colabfold_envdb = false
params.output_rmsd_aligned = false

// MPNN params - new (mpnn_*)
params.mpnn_model_type = 'protein_mpnn'
params.mpnn_legacy_weights = true
params.mpnn_designed_chains = 'A'
params.mpnn_batch_size = false
params.mpnn_temperature = 0.1
params.mpnn_structure_noise = 0
params.mpnn_omit = 'CX'
params.mpnn_checkpoint_path = false

params.require_gpu = true
params.gpu_devices = ''
params.gpu_allocation_detect_process_regex = '(boltz predict|rfd3 design|rf3 fold)'

include { UNIQUE_ID } from '../modules/local/common/unique_id'
include { RFDIFFUSION3 } from '../modules/local/rfd3/rfdiffusion3'
include { GENERATE_RFD3_CONFIG } from '../modules/local/rfd3/generate_rfd3_config'
include { MPNN } from '../modules/local/rfd3/mpnn'
include { ROSETTAFOLD3 } from '../modules/local/rfd3/rosettafold3'
include { RFD3_RMSD } from '../modules/local/rfd3/rfd3_rmsd'
include { COMBINE_RFD3_SCORES } from '../modules/local/rfd3/combine_rfd3_scores'
include { FILTER_DESIGNS as RFD3_FILTER_DESIGNS } from '../modules/local/rfd3/filter_designs'
include { buildMpnnArgs; normaliseContigToV3 } from '../modules/local/rfd3/rfd3_utils'
include { BOLTZ_REFOLD_SCORING_RFD3 } from '../subworkflows/local/boltz_refold_scoring_rfd3'

workflow RFD3 {

    main:

    if (params.rfd3_config && params.input_pdb) {
        throw new Exception('--rfd3_config and --input_pdb are mutually exclusive. Use --rfd3_config (input path is read from the config) or --input_pdb with --contigs.')
    }

    // Show help message
    if (params.input_pdb == false && params.rfd3_config == false) {
        log.info(
            """
        ==================================================================
        PROTEIN BINDER DESIGN PIPELINE - RFDiffusion3
        ==================================================================

        Required (one of):
            --rfd3_config         Path to RFDiffusion3 JSON/YAML config file
                                  (input path is read from the config, e.g. "designname": {"input": "path/to/target.pdb"})
            --input_pdb           Input PDB/CIF file (params mode only)

        Config mode (--rfd3_config):
            Input file path is taken from the config; do not use --input_pdb.

        Params mode (--input_pdb, generates config automatically):
            --contigs             Contig string for rfd3 (v3 format: "A17-131,/0,50-120")
                                  v1 format is also accepted and auto-translated [default: ${params.contigs}]
            --hotspot_res         Hotspot residues, eg "A56,A115,A123" [default: ${params.hotspot_res}]

        RFDiffusion3 options:
            --outdir              Output directory [default: ${params.outdir}]
            --design_name         Name of the design [default: ${params.design_name}]
            --rfd3_n_designs      Total number of designs [default: ${params.rfd3_n_designs}]
            --rfd3_batch_size     Designs per batch (diffusion_batch_size) [default: ${params.rfd3_batch_size}]
            --rfd3_step_scale     inference_sampler.step_scale [default: ${params.rfd3_step_scale}]
            --rfd3_gamma_0        inference_sampler.gamma_0 [default: ${params.rfd3_gamma_0}]
            --rfd3_allow_loopy   Allow loopy designs (disable is_non_loopy) [default: ${params.rfd3_allow_loopy}]
            --rfd3_extra_args     Additional CLI arguments for rfd3 [default: ${params.rfd3_extra_args}]
            --rfd3_filters        Semicolon-separated filters for RFD3 backbones (binder = chain B), e.g. "rg<25" [default: disabled]

        MPNN options (new names / legacy names):
            --mpnn_model_type / (n/a)                  Model type [default: ${params.mpnn_model_type}]
            --mpnn_legacy_weights / (n/a)              Use legacy weights [default: ${params.mpnn_legacy_weights}]
            --mpnn_designed_chains / (n/a)             Chains to redesign [default: ${params.mpnn_designed_chains}]
            --mpnn_batch_size / --pmpnn_seqs_per_struct  Sequences per structure [default: ${params.pmpnn_seqs_per_struct}]
            --mpnn_temperature / --pmpnn_temperature   Sampling temperature [default: ${params.mpnn_temperature}]
            --mpnn_structure_noise / --pmpnn_augment_eps  Structure noise [default: ${params.mpnn_structure_noise}]
            --mpnn_omit / --pmpnn_omit_aas             Omit residue types (1-letter eg "CX") [default: ${params.mpnn_omit}]
            --mpnn_checkpoint_path / --pmpnn_weights   Custom weights path [default: ${params.mpnn_checkpoint_path}]

        Refolding options:
            --refold_with         Comma-separated list of refolding methods (valid options: 'boltz') [default: ${params.refold_with}]
            --refold_max          Maximum designs to refold [default: ${params.refold_max}]
            --refold_filter_sort  Metric to sort by before refolding. Use '-' prefix for descending. [default: ${params.refold_filter_sort}]
            --refold_create_target_msa  Create target MSA for refolding structure prediction [default: ${params.refold_create_target_msa}]

        Other options:
            --require_gpu         Fail tasks without a GPU [default: ${params.require_gpu}]
            --gpu_devices         GPU devices to use (comma-separated or 'all') [default: ${params.gpu_devices}]

        """.stripIndent()
        )
        exit(1)
    }
    def refold_methods = params.refold_with ? params.refold_with.toString().split(',').collect{it.trim()} : []
    println("Refold methods: ${refold_methods}")

    def ch_input_pdb
    if (params.rfd3_config) {
        def config_file = file(params.rfd3_config)
        def config_dir = config_file.parent.toString()
        def parse_cmd = ["python3", "${projectDir}/bin/rfd3/stage_rfd3_config.py", "parse-inputs", config_file.toString(), "--config-dir", config_dir]
        def parse_output = parse_cmd.execute().text.trim()
        def input_paths = parse_output.split('\n').findAll { it.trim() }
        if (input_paths.isEmpty()) {
            throw new Exception("No 'input' path found in --rfd3_config ${params.rfd3_config}")
        }
        ch_input_pdb = Channel.fromPath(input_paths[0]).first()
    } else {
        ch_input_pdb = Channel.fromPath(params.input_pdb).first()
    }

    // Generate unique ID for this run
    UNIQUE_ID()
    ch_unique_id = UNIQUE_ID.out.id_file.map { it.text.trim() }
    ch_design_name = ch_unique_id.map { uid -> "${params.design_name}_${uid}" }

    // Calculate number of batches
    def n_batches = Math.ceil(params.rfd3_n_designs / params.rfd3_batch_size).toInteger()
    if (n_batches < 1) { n_batches = 1 }

    if (params.rfd3_config) {
        // Config mode: user provides their own JSON/YAML config
        ch_config = Channel.fromPath(params.rfd3_config).first()

        ch_rfd3_startnum = Channel.of(0..n_batches - 1)

        RFDIFFUSION3(
            ch_config,
            ch_input_pdb,
            ch_design_name,
            params.rfd3_batch_size,
            1,  // n_batches=1 per task, we parallelise via Nextflow
            ch_rfd3_startnum,
            params.rfd3_step_scale,
            params.rfd3_gamma_0,
            params.rfd3_extra_args,
        )
    }
    else {
        // Params mode: generate config from --contigs, --hotspot_res, etc.
        if (!params.contigs) {
            throw new Exception('Either --rfd3_config or --contigs must be provided')
        }

        GENERATE_RFD3_CONFIG(
            params.design_name,
            ch_input_pdb,
            normaliseContigToV3(params.contigs),
            params.hotspot_res ?: false,
            false,  // partial_t - not used in de novo design
            params.rfd3_allow_loopy,
        )

        ch_rfd3_startnum = Channel.of(0..n_batches - 1)

        RFDIFFUSION3(
            GENERATE_RFD3_CONFIG.out.config_json,
            ch_input_pdb,
            ch_design_name,
            params.rfd3_batch_size,
            1,  // n_batches=1 per task, we parallelise via Nextflow
            ch_rfd3_startnum,
            params.rfd3_step_scale,
            params.rfd3_gamma_0,
            params.rfd3_extra_args,
        )
    }

    // Build MPNN CLI args from resolved params
    def mpnn_args = buildMpnnArgs(params)

    ch_rfd3_cifs = RFDIFFUSION3.out.cifs.flatten()

    if (params.rfd3_filters) {
        RFD3_FILTER_DESIGNS(
            ch_rfd3_cifs,
            params.rfd3_filters,
            'B',
            'filter',
        )
        ch_backbones = RFD3_FILTER_DESIGNS.out.accepted
    }
    else {
        ch_backbones = ch_rfd3_cifs
    }

    MPNN(
        ch_backbones,
        mpnn_args,
    )

    ch_mpnn_with_meta = MPNN.out.cifs.flatten().map { c -> tuple([id: c.name], c) }

    // Run RosettaFold3 structure prediction on each MPNN-designed structure (run unique id as value channel)
    ROSETTAFOLD3(ch_mpnn_with_meta, ch_unique_id)

    ch_rmsd_input = ch_mpnn_with_meta.join(ROSETTAFOLD3.out.refolded_cif)
    RFD3_RMSD(ch_rmsd_input)

    ch_boltz_complex = Channel.empty()
    ch_boltz_monomer = Channel.empty()

    if (refold_methods.contains('boltz')) {
        ch_boltz_input = ROSETTAFOLD3.out.refolded_cif
            .join(ROSETTAFOLD3.out.scores_with_meta)
            .map { meta, cif, score_tsv ->
                def lines = score_tsv.readLines()
                def header = lines[0].split('\t')
                def values = lines[1].split('\t')
                def scoreMap = [:]
                for (int i = 0; i < header.size(); i++) {
                    scoreMap[header[i]] = values[i]
                }
                return [meta, cif, scoreMap]
            }
            .toSortedList { a, b ->
                def sort_key = params.refold_filter_sort ?: 'pair_pae_min'
                def descending = false
                if (sort_key.startsWith('-')) {
                    descending = true
                    sort_key = sort_key.substring(1)
                }
                def valA = a[2][sort_key]
                def valB = b[2][sort_key]
                
                valA = valA && valA != 'None' && valA != '' ? valA.toDouble() : (descending ? -Double.MAX_VALUE : Double.MAX_VALUE)
                valB = valB && valB != 'None' && valB != '' ? valB.toDouble() : (descending ? -Double.MAX_VALUE : Double.MAX_VALUE)

                return descending ? valB <=> valA : valA <=> valB
            }
            .flatMap { list ->
                if (params.refold_max) {
                    return list.take(params.refold_max as int)
                }
                return list
            }
            .map { meta, cif, scoreMap -> tuple(meta, cif) }

        BOLTZ_REFOLD_SCORING_RFD3(
            ch_boltz_input,
            'B', // binder_chain
            'A', // target_chain
            params.refold_max,
            params.refold_create_target_msa,
            params.refold_use_msa_server,
            params.refold_target_fasta,
            params.refold_target_templates,
            params.colabfold_envdb,
            params.uniref30,
            params.outdir
        )
        ch_boltz_complex = BOLTZ_REFOLD_SCORING_RFD3.out.boltz_scores_complex
        ch_boltz_monomer = BOLTZ_REFOLD_SCORING_RFD3.out.boltz_scores_monomer
    }

    ch_rmsd_target_aligned_binder = RFD3_RMSD.out.rmsd_target_aligned_binder
        .map { meta, tsv -> tsv }
        .collectFile(name: 'rmsd_target_aligned_binder.tsv', storeDir: "${params.outdir}/rfd3/rosettafold3/rmsd", keepHeader: true, skip: 1)
    ch_rmsd_complex = RFD3_RMSD.out.rmsd_complex
        .map { meta, tsv -> tsv }
        .collectFile(name: 'rmsd_complex.tsv', storeDir: "${params.outdir}/rfd3/rosettafold3/rmsd", keepHeader: true, skip: 1)
    ch_rmsd_binder_aligned_binder = RFD3_RMSD.out.rmsd_binder_aligned_binder
        .map { meta, tsv -> tsv }
        .collectFile(name: 'rmsd_binder_aligned_binder.tsv', storeDir: "${params.outdir}/rfd3/rosettafold3/rmsd", keepHeader: true, skip: 1)
    ch_rmsd_target_aligned_target = RFD3_RMSD.out.rmsd_target_aligned_target
        .map { meta, tsv -> tsv }
        .collectFile(name: 'rmsd_target_aligned_target.tsv', storeDir: "${params.outdir}/rfd3/rosettafold3/rmsd", keepHeader: true, skip: 1)

    ch_rmsd_tuple = ch_rmsd_target_aligned_binder
        .combine(ch_rmsd_complex)
        .combine(ch_rmsd_binder_aligned_binder)
        .combine(ch_rmsd_target_aligned_target)
        .ifEmpty(Channel.of(tuple(
            file("${projectDir}/assets/dummy_files/empty"),
            file("${projectDir}/assets/dummy_files/empty"),
            file("${projectDir}/assets/dummy_files/empty"),
            file("${projectDir}/assets/dummy_files/empty"),
        )))

    ch_rf3_scores_merged = ROSETTAFOLD3.out.scores
        .collectFile(name: 'rf3_scores.tsv', keepHeader: true, skip: 1)
    ch_rfd3_scores_merged = RFDIFFUSION3.out.scores
        .collectFile(name: 'rfd3_scores.tsv', keepHeader: true, skip: 1)

    ch_combine_input = ch_rf3_scores_merged
        .combine(ch_rfd3_scores_merged)
        .combine(ch_rmsd_tuple)
        .combine(ch_boltz_complex.ifEmpty(file("${projectDir}/assets/dummy_files/empty")))
        .combine(ch_boltz_monomer.ifEmpty(file("${projectDir}/assets/dummy_files/empty")))
        .map { it ->
            def f = it.flatten()
            tuple(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7])
        }
    COMBINE_RFD3_SCORES(ch_combine_input)

    emit:
    rfd3_cifs = RFDIFFUSION3.out.cifs
    mpnn_cifs = MPNN.out.cifs
    mpnn_fastas = MPNN.out.fastas
    rf3_results = ROSETTAFOLD3.out.results
    combined_scores = COMBINE_RFD3_SCORES.out.combined_scores
    rfd3_rmsd_target_aligned_binder = ch_rmsd_target_aligned_binder
    rfd3_rmsd_complex = ch_rmsd_complex
    rfd3_rmsd_binder_aligned_binder = ch_rmsd_binder_aligned_binder
    rfd3_rmsd_target_aligned_target = ch_rmsd_target_aligned_target
}
