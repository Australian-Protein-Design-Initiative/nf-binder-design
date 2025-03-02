#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
eg

nextflow run main.nf \
  --input_pdb 'designs/*.pdb' \
  --binder_chain 'A' \
  --target_contigs 'B1-100' \
  --rfd_batch_size=10 \
  --rfd_n_partial_per_binder=10 \
  -resume \
  -with-report report_$(date +%Y%m%d_%H%M%S).html \
  -with-trace trace_$(date +%Y%m%d_%H%M%S).txt

*/

// Default parameters
params.input_pdb = false
params.outdir = "results"

params.design_name = "fuzzed_ppi"
params.binder_chain = "A" // eg, our fixed target chain - usually A ?
params.target_contigs = 'auto' // 'auto' to detect from PDB, or eg "B10-110" for fixed target chain B, residues 10-110
params.rfd_batch_size = 10
params.rfd_n_partial_per_binder = 10
params.rfd_model_path = false // "models/rfdiffusion/Complex_beta_ckpt.pt"
params.rfd_config = "base"
params.rfd_noise_scale = 0
params.rfd_partial_T = 20
params.rfd_extra_args = ""

params.pmpnn_relax_cycles = 0
params.pmpnn_seqs_per_struct = 1
params.pmpnn_weights = false
params.pmpnn_temperature = 0.000001
params.pmpnn_augment_eps = 0

params.require_gpu = true

// Validate numeric parameters
def validate_numeric = { param_name, value ->
    if (!(value instanceof Number)) {
        error "Parameter $param_name must be a number, got: $value (${value.getClass().getName()})"
    }
}

validate_numeric('rfd_batch_size', params.rfd_batch_size)
validate_numeric('rfd_n_partial_per_binder', params.rfd_n_partial_per_binder)
validate_numeric('rfd_noise_scale', params.rfd_noise_scale)
validate_numeric('rfd_partial_T', params.rfd_partial_T)
validate_numeric('pmpnn_relax_cycles', params.pmpnn_relax_cycles)
validate_numeric('pmpnn_seqs_per_struct', params.pmpnn_seqs_per_struct)

// Show help message
if (params.input_pdb == false) {
    log.info"""
    ==================================================================
    🧬 PROTEIN BINDER DESIGN PIPELINE 🧬
    ==================================================================
    
    Required arguments:
        --input_pdb           Input PDBs file for the binders to diffuse (* glob accepted)

    Optional arguments:
        --design_name         Name of the design, used for output file prefixes [default: ${params.design_name}]
        --target_contigs      Contig map for target chain(s) - 'auto' to detect from PDB, or specify manually [default: ${params.target_contigs}]
        --binder_chain        Chain ID of the binder chain [default: ${params.binder_chain}]
        --rfd_batch_size      Number of designs per batch [default: ${params.rfd_batch_size}]
        --rfd_n_partial_per_binder Number of partial diffused designs per binder [default: ${params.rfd_n_partial_per_binder}]
        --rfd_model_path      Path to RFdiffusion model checkpoint file - you probaby don't want to set this manually [default: ${params.rfd_model_path}]
        --rfd_extra_args      Extra arguments for RFdiffusion [default: ${params.rfd_extra_args}]
        --rfd_config          'base', 'symmetry' or a path to a YAML file [default: ${params.rfd_config_name}]
        --rfd_partial_T       Number of timesteps to run partial diffusion for (lower = less diffusion) [default: ${params.rfd_partial_T}]
        --pmpnn_relax_cycles  Number of relax cycles for ProteinMPNN [default: ${params.pmpnn_relax_cycles}]
        --pmpnn_seqs_per_struct Number of sequences per structure for ProteinMPNN [default: ${params.pmpnn_seqs_per_struct}]
        --pmpnn_weights       Path to ProteinMPNN weights file (leave unset to use default weights) [default: ${params.pmpnn_weights}]
        --pmpnn_temperature   Temperature for ProteinMPNN [default: ${params.pmpnn_temperature}]
        --pmpnn_augment_eps   Variance of random noise to add to the atomic coordinates ProteinMPNN [default: ${params.pmpnn_augment_eps}]
        --outdir              Output directory [default: ${params.outdir}]

        --require_gpu       Fail tasks that go too slow without a GPU if no GPU is detected [default: ${params.require_gpu}]

    """.stripIndent()
    exit 1
}

include { RFDIFFUSION } from './modules/rfdiffusion'
include { RFDIFFUSION_PARTIAL } from './modules/rfdiffusion_partial'
include { SILENT_FROM_PDBS } from './modules/silentfrompdbs' 
include { DL_BINDER_DESIGN_PROTEINMPNN } from './modules/dl_binder_design'
include { AF2_INITIAL_GUESS } from './modules/af2_initial_guess'
include { GET_CONTIGS } from './modules/get_contigs'
include { RENUMBER_RESIDUES } from './modules/renumber_residues'
workflow {

    if (!params.input_pdb) {
        throw new Exception("--input_pdb must be provided")
    }

    // File inputs - converted to value channels with .first()
    // so these channels infinitely produce the file on demand
    ch_rfd_config = Channel.fromPath(params.rfd_config).first()
    ch_input_pdb = Channel.fromPath(params.input_pdb)
    if (params.rfd_model_path) {
        ch_rfd_model_path = Channel.fromPath(params.rfd_model_path).first()
    } else {
        ch_rfd_model_path = Channel.value(false)
    }

    // We renumber the residues in each chain to be 1 to chain_length
    // irrespective of missing non-sequential residue numbers etc
    // The --binder_chain is always set to chain A, with the non-diffusable chains
    // coming after it.
    ch_preprocessed_pdb = ch_input_pdb
        .map { pdb -> [pdb, params.binder_chain] } 
        | RENUMBER_RESIDUES

    // Get contigs string for each input PDB
    // NOTE/HACK: We hardcode 'A' here, since RENUMBER_RESIDUES always 
    // makes the binder chain 'A'
    ch_contigs = ch_preprocessed_pdb
        .map { pdb -> [pdb, 'A']} //params.binder_chain] }
        | GET_CONTIGS

    //ch_contigs.view()

    // Create batches with contigs
    ch_contigs
        | map { input_pdb, contigs ->
            def num_batches = (params.rfd_n_partial_per_binder / params.rfd_batch_size).toInteger()
            // For each batch, create a tuple of (pdb, contigs, start_num)
            (0..<num_batches).collect { batchNum ->
                tuple(input_pdb, contigs, batchNum * params.rfd_batch_size)
            }
        }
        | flatMap()  // Flatten the lists of tuples into individual tuples
        | set { ch_rfd_jobs }

    // Run RFdiffusion with partial diffusion in batches
    RFDIFFUSION_PARTIAL(
        ch_rfd_config,
        ch_rfd_jobs.map { input_pdb, contigs, start -> input_pdb },  // Extract PDB path
        ch_rfd_model_path,
        ch_rfd_jobs.map { input_pdb, contigs, start -> contigs },  // Extract contigs string
        params.rfd_batch_size,
        ch_rfd_jobs.map { input_pdb, contigs, start -> start }  // Extract start number
    )
    ch_rfd_backbone_models = RFDIFFUSION_PARTIAL.out.pdbs

    // Convert PDBs to silent file
    // SILENT_FROM_PDBS(
    //     RFDIFFUSION.out.pdbs.collect()
    // )

    // Run ProteinMPNN (dl_binder_design) on backbone-only PDBs
    DL_BINDER_DESIGN_PROTEINMPNN(
        ch_rfd_backbone_models,
        params.pmpnn_seqs_per_struct,
        params.pmpnn_relax_cycles,
        params.pmpnn_weights,
        params.pmpnn_temperature,
        params.pmpnn_augment_eps
    )


    // // Run AF2 initial guess to build/refine sidechains for compatible sequence
    AF2_INITIAL_GUESS(
        DL_BINDER_DESIGN_PROTEINMPNN.out.pdbs
    )

    // TODO: Use bin/af2_combine_scores.py to collate DL_BINDER_DESIGN_PROTEINMPNN.out.scores 
    // tables into a single TSV table eg many tables, multi-space sep, like:
    /*
    SCORE:     binder_aligned_rmsd pae_binder pae_interaction pae_target plddt_binder plddt_target plddt_total target_aligned_rmsd time description
    SCORE:       14.374   11.205   27.694    4.642   55.734   93.869   88.411   28.443  189.208        design_ppi_2_dldesign_0_af2pred
    SCORE:       13.619   11.015   27.882    4.714   55.609   93.577   88.143   30.422   37.691        design_ppi_2_dldesign_1_af2pred
    */
} 
