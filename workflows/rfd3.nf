#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
RFDiffusion3-based binder design workflow

Usage via main.nf:
  nextflow run main.nf --method rfd3 --rfd3_config config.json --input_pdb target.pdb
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

// RosettaFold3 params
//params.rf3_ckpt_path = '/weights/rf3_foundry_01_24_latest_remapped.ckpt'
params.rf3_ckpt_path = '/models/foundry/rf3_foundry_01_24_latest_remapped.ckpt'

// MPNN params - legacy (pmpnn_*) for compatibility with rfd workflow
params.pmpnn_seqs_per_struct = 1
params.pmpnn_temperature = false
params.pmpnn_augment_eps = false
params.pmpnn_omit_aas = false
params.pmpnn_weights = false

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
params.gpu_allocation_detect_process_regex = '(python.*/app/dl_binder_design/af2_initial_guess/predict\\.py|python.*/app/BindCraft/bindcraft\\.py|boltz predict|rfd3)'

include { UNIQUE_ID } from '../modules/local/common/unique_id'
include { RFDIFFUSION3 } from '../modules/local/rfd3/rfdiffusion3'
include { GENERATE_RFD3_CONFIG } from '../modules/local/rfd3/generate_rfd3_config'
include { MPNN } from '../modules/local/rfd3/mpnn'
include { ROSETTAFOLD3 } from '../modules/local/rfd3/rosettafold3'
include { buildMpnnArgs; normaliseContigToV3 } from '../modules/local/rfd3/rfd3_utils'

workflow RFD3 {

    main:

    // Show help message
    if (params.input_pdb == false && params.rfd3_config == false) {
        log.info(
            """
        ==================================================================
        PROTEIN BINDER DESIGN PIPELINE - RFDiffusion3
        ==================================================================

        Required arguments:
            --input_pdb           Input PDB/CIF file for the target

        Config mode (recommended for advanced use):
            --rfd3_config         Path to RFDiffusion3 JSON/YAML config file
                                  (see https://rosettacommons.github.io/foundry/models/rfd3/input.html)

        Params mode (generates config automatically):
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

        MPNN options (new names / legacy names):
            --mpnn_model_type / (n/a)                  Model type [default: ${params.mpnn_model_type}]
            --mpnn_legacy_weights / (n/a)              Use legacy weights [default: ${params.mpnn_legacy_weights}]
            --mpnn_designed_chains / (n/a)             Chains to redesign [default: ${params.mpnn_designed_chains}]
            --mpnn_batch_size / --pmpnn_seqs_per_struct  Sequences per structure [default: ${params.pmpnn_seqs_per_struct}]
            --mpnn_temperature / --pmpnn_temperature   Sampling temperature [default: ${params.mpnn_temperature}]
            --mpnn_structure_noise / --pmpnn_augment_eps  Structure noise [default: ${params.mpnn_structure_noise}]
            --mpnn_omit / --pmpnn_omit_aas             Omit residue types (1-letter eg "CX") [default: ${params.mpnn_omit}]
            --mpnn_checkpoint_path / --pmpnn_weights   Custom weights path [default: ${params.mpnn_checkpoint_path}]

        Other options:
            --require_gpu         Fail tasks without a GPU [default: ${params.require_gpu}]
            --gpu_devices         GPU devices to use (comma-separated or 'all') [default: ${params.gpu_devices}]

        """.stripIndent()
        )
        exit(1)
    }

    if (!params.input_pdb) {
        throw new Exception('--input_pdb must be provided')
    }

    // Generate unique ID for this run
    UNIQUE_ID()
    ch_unique_id = UNIQUE_ID.out.id_file.map { it.text.trim() }
    ch_design_name = ch_unique_id.map { uid -> "${params.design_name}_${uid}" }

    ch_input_pdb = Channel.fromPath(params.input_pdb).first()

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

    // Run MPNN inverse folding on each backbone structure
    ch_backbones = RFDIFFUSION3.out.cifs.flatten()

    MPNN(
        ch_backbones,
        mpnn_args,
    )

    // Run RosettaFold3 structure prediction on each MPNN-designed structure
    ROSETTAFOLD3(
        MPNN.out.cifs.flatten(),
    )

    emit:
    rfd3_cifs = RFDIFFUSION3.out.cifs
    mpnn_cifs = MPNN.out.cifs
    mpnn_fastas = MPNN.out.fastas
    rf3_results = ROSETTAFOLD3.out.results
}
