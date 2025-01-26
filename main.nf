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
params.outdir = "results"
params.design_name = "design_ppi"
params.contigs = "[A371-508/A753-883/A946-1118/A1135-1153/0 70-100]"
params.hotspot_res = "[A473,A995,A411,A421]"
params.rfd_batch_size = 1
params.rfd_n_designs = 2
params.rfd_model_path = false // "models/rfdiffusion/Complex_beta_ckpt.pt"
params.rfd_config = "base"
params.rfd_noise_scale = 0
params.rfd_extra_args = ""
params.pmpnn_relax_cycles = 0
params.pmpnn_seqs_per_struct = 1
params.require_gpu = true

// Set to a path of existing RFDiffusion backbone models, skips running RFDiffusion
params.rfd_backbone_models = false

// Show help message
if (params.input_pdb == false) {
    log.info"""
    ==================================================================
    🧬 PROTEIN BINDER DESIGN PIPELINE 🧬
    ==================================================================
    
    Required arguments:


    Optional arguments:
        --input_pdb           Input PDB file for the target
        --rfd_backbone_models Path to existing RFDiffusion backbone models - skips running RFDiffusion and uses these instead
        --design_name         Name of the design, used for output file prefixes [default: ${params.design_name}]
        --contigs             Contig map for RFdiffusion [default: ${params.contigs}]
        --hotspot_res         Hotspot residues [default: ${params.hotspot_res}]
        --rfd_batch_size      Number of designs per batch [default: ${params.rfd_batch_size}]
        --rfd_model_path      Path to RFdiffusion model checkpoint file - leaving unset will allow RFDiffusion to choose based on other parameters [default: ${params.rfd_model_path}]
        --rfd_extra_args      Extra arguments for RFdiffusion [default: ${params.rfd_extra_args}]
        --rfd_config          'base', 'symmetry' or a path to a YAML file [default: ${params.rfd_config_name}]
        --outdir              Output directory [default: ${params.outdir}]

        --require_gpu       Fail tasks that go too slow without a GPU if no GPU is detected [default: ${params.require_gpu}]

    """.stripIndent()
    exit 1
}

include { RFDIFFUSION } from './modules/rfdiffusion'
include { SILENT_FROM_PDBS } from './modules/silentfrompdbs' 
include { DL_BINDER_DESIGN_PROTEINMPNN } from './modules/dl_binder_design'
include { AF2_INITIAL_GUESS } from './modules/af2_initial_guess'

workflow {

    if (!params.input_pdb && !params.rfd_backbone_models) {
        throw new Exception("Either --input_pdb or --rfd_backbone_models must be provided")
    }

    // File inputs - converted to value channels with .first()
    // so these channels infinitely produce the file on demand
    ch_rfd_config = Channel.fromPath(params.rfd_config).first()
    ch_input_pdb = Channel.fromPath(params.input_pdb).first()
    if (params.rfd_model_path) {
        ch_rfd_model_path = Channel.fromPath(params.rfd_model_path).first()
    } else {
        ch_rfd_model_path = Channel.value(false)
    }

    // Create channel of start numbers for batches
    ch_rfd_startnum = Channel
        .of(0..params.rfd_n_designs-1)
        .filter { v -> v % params.rfd_batch_size == 0 }
        //.view()
    
    if (params.rfd_backbone_models) {
        ch_rfd_backbone_models = Channel.fromPath("${params.rfd_backbone_models}/*.pdb")
    } else {
        // Run RFdiffusion in batches
        RFDIFFUSION(
            ch_rfd_config,
            ch_input_pdb,
            ch_rfd_model_path,
            params.contigs,
            params.hotspot_res,
            params.rfd_batch_size,
            ch_rfd_startnum
        )
        ch_rfd_backbone_models = RFDIFFUSION.out.pdbs
    }

    // Convert PDBs to silent file
    // SILENT_FROM_PDBS(
    //     RFDIFFUSION.out.pdbs.collect()
    // )

    // Run ProteinMPNN (dl_binder_design) on backbone-only PDBs
    DL_BINDER_DESIGN_PROTEINMPNN(
        ch_rfd_backbone_models,
        params.pmpnn_relax_cycles,
        params.pmpnn_seqs_per_struct
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
