#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
eg

nextflow run main.nf \
  --input_pdb target.pdb \
  --rfd_n_designs=2 \
  -resume \
  -with-report report_$(date +%Y%m%d_%H%M%S).html \
  -with-trace trace_$(date +%Y%m%d_%H%M%S).txt

*/

// Default parameters
params.input_pdb = false
params.outdir = 'results'

params.design_name = 'design_ppi'
params.contigs = ''         // "[A371-508/A753-883/A946-1118/A1135-1153/0 70-100]"
params.hotspot_res = false  // "[A473,A995,A411,A421]"
params.rfd_batch_size = 1
params.rfd_n_designs = 2
params.rfd_model_path = false // "models/rfdiffusion/Complex_beta_ckpt.pt"
params.rfd_config = false // "base"
params.rfd_noise_scale = 0
params.rfd_extra_args = ''
params.rfd_compress_trajectories = true

params.pmpnn_relax_cycles = 0
params.pmpnn_seqs_per_struct = 1
params.pmpnn_weights = false
params.pmpnn_temperature = 0.000001
params.pmpnn_augment_eps = 0
params.pmpnn_omit_aas = 'CX'

params.max_rg = false
params.rfd_filters = false

params.af2ig_recycle = 3

params.require_gpu = true
params.gpu_device = 'all'

// Set to a path of existing RFDiffusion backbone models, skips running RFDiffusion
params.rfd_backbone_models = false

// Show help message
if (params.input_pdb == false && params.rfd_backbone_models == false) {
    log.info"""
    ==================================================================
    🧬 PROTEIN BINDER DESIGN PIPELINE 🧬
    ==================================================================

    Required arguments:
        --input_pdb           Input PDB file for the target

        __or__

        --rfd_backbone_models Existing RFDiffusion backbone models - skips running RFDiffusion and uses these instead (glob accepted, eg 'results/rfdiffusion/pdbs/*.pdb)

        Optional arguments:
            --outdir              Output directory [default: ${params.outdir}]
            --design_name         Name of the design, used for output file prefixes [default: ${params.design_name}]
            --contigs             Contig map for RFdiffusion [default: ${params.contigs}]
            --hotspot_res         Hotspot residues, eg "[A473,A995,A411,A421]" [default: ${params.hotspot_res}]
            --rfd_batch_size      Number of designs per batch [default: ${params.rfd_batch_size}]
            --rfd_model_path      Path to RFdiffusion model checkpoint file - leaving unset will allow RFDiffusion to choose based on other parameters [default: ${params.rfd_model_path}]
            --rfd_extra_args      Extra arguments for RFdiffusion [default: ${params.rfd_extra_args}]
            --rfd_config          'base', 'symmetry' or a path to a YAML file [default: ${params.rfd_config_name}]
            --rfd_compress_trajectories Compress trajectories with gzip [default: ${params.rfd_compress_trajectories}]
            --pmpnn_weights       Path to ProteinMPNN weights file (leave unset to use default weights) [default: ${params.pmpnn_weights}]
            --pmpnn_temperature   Temperature for ProteinMPNN [default: ${params.pmpnn_temperature}]
            --pmpnn_augment_eps   Variance of random noise to add to the atomic coordinates ProteinMPNN [default: ${params.pmpnn_augment_eps}]
            --pmpnn_relax_cycles  Number of relax cycles for ProteinMPNN [default: ${params.pmpnn_relax_cycles}]
            --pmpnn_seqs_per_struct Number of sequences per structure for ProteinMPNN [default: ${params.pmpnn_seqs_per_struct}]
            --pmpnn_omit_aas      A string of all residue types (one letter case-insensitive) that should not appear in the design [default: ${params.pmpnn_omit_aas}]
            --max_rg              Maximum radius of gyration for backbone filtering [default: disabled]
            --rfd_filters         Semicolon-separated list of filters for RFDiffusion backbones, eg "rg<25;compactness>0.8" [default: disabled]

        --af2ig_recycle       Number of recycle cycles for AF2 initial guess [default: ${params.af2ig_recycle}]
        --require_gpu         Fail tasks that go too slow without a GPU if no GPU is detected [default: ${params.require_gpu}]
        --gpu_device          GPU device to use [default: ${params.gpu_device}]

    """.stripIndent()
    exit 1
}

include { RFDIFFUSION } from './modules/rfdiffusion'
include { SILENT_FROM_PDBS } from './modules/silentfrompdbs'
include { DL_BINDER_DESIGN_PROTEINMPNN } from './modules/dl_binder_design'
include { AF2_INITIAL_GUESS } from './modules/af2_initial_guess'
include { COMBINE_SCORES } from './modules/combine_scores'
include { UNIQUE_ID } from './modules/unique_id'
include { FILTER_DESIGNS } from './modules/filter_designs'

workflow {
    // Generate unique ID for this run
    UNIQUE_ID()
    ch_unique_id = UNIQUE_ID.out.id_file.map { it.text.trim() }

    if (!params.input_pdb && !params.rfd_backbone_models) {
        throw new Exception('Either --input_pdb or --rfd_backbone_models must be provided')
    }

    // File inputs - converted to value channels with .first()
    // so these channels infinitely produce the file on demand
    if (params.rfd_config) {
        ch_rfd_config = Channel.fromPath(params.rfd_config).first()
    } else {
        ch_rfd_config = Channel.value(false)
    }
    if (params.input_pdb) {
        ch_input_pdb = Channel.fromPath(params.input_pdb).first()
    } else {
        ch_input_pdb = Channel.value(false)
    }
    if (params.rfd_model_path) {
        ch_rfd_model_path = Channel.fromPath(params.rfd_model_path).first()
    } else {
        ch_rfd_model_path = Channel.value(false)
    }

    // Create channel of start numbers for batches
    ch_rfd_startnum = Channel
        .of(0..params.rfd_n_designs - 1)
        .filter { v -> v % params.rfd_batch_size == 0 }
        //.view()

    if (params.rfd_backbone_models) {
        ch_rfd_backbone_models = Channel.fromPath("${params.rfd_backbone_models}")
    } else {
        // Run RFdiffusion in batches
        RFDIFFUSION(
            ch_rfd_config,
            ch_input_pdb,
            ch_rfd_model_path,
            params.contigs,
            params.hotspot_res,
            params.rfd_batch_size,
            ch_rfd_startnum,
            ch_unique_id,
            params.gpu_device
        )
        ch_rfd_backbone_models = RFDIFFUSION.out.pdbs.flatten()
    }

    if (params.rfd_filters) {
        FILTER_DESIGNS(
            ch_rfd_backbone_models,
            params.rfd_filters,
            'A',          // binder_chain is always "A" for RFDiffusion binder design
            'rfdiffusion' // step name
        )
        ch_filtered_backbones = FILTER_DESIGNS.out.accepted
    } else {
        ch_filtered_backbones = ch_rfd_backbone_models
    }

    // Create a channel that repeats each PDB params.pmpnn_seqs_per_struct times
    // and pairs it with an index from 0 to pmpnn_seqs_per_struct-1
    ch_pmpnn_inputs = ch_filtered_backbones
        | combine(Channel.of(0..(params.pmpnn_seqs_per_struct - 1)))

    // Run ProteinMPNN (dl_binder_design) on backbone-only PDBs
    DL_BINDER_DESIGN_PROTEINMPNN(
        ch_pmpnn_inputs.map { pdb, idx -> pdb },  // Extract just the PDB for the input path
        Channel.value(1),  // Always use seqs_per_struct=1
        params.pmpnn_relax_cycles,
        params.pmpnn_weights,
        params.pmpnn_temperature,
        params.pmpnn_augment_eps,
        ch_pmpnn_inputs.map { pdb, idx -> idx }  // Pass the index for unique naming
    )

    // Run AF2 initial guess to build/refine sidechains for compatible sequence
    AF2_INITIAL_GUESS(
        DL_BINDER_DESIGN_PROTEINMPNN.out.pdbs,
        params.gpu_device
    )

    // Combine all the score files into a single TSV file
    COMBINE_SCORES(
        AF2_INITIAL_GUESS.out.scores.collect(),
        AF2_INITIAL_GUESS.out.pdbs.collect()
    )
}
