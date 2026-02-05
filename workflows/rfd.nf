#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
RFDiffusion-based binder design workflow

Usage via main.nf:
  nextflow run main.nf --method rfd --input_pdb target.pdb --rfd_n_designs=10

*/

// Default parameters
params.input_pdb = false
params.outdir = 'results'

params.design_name = 'design_ppi'
params.contigs = ''
params.hotspot_res = false
params.rfd_batch_size = 1
params.rfd_n_designs = 2
params.rfd_model_path = false
params.rfd_config = false
params.rfd_noise_scale = 0
params.rfd_extra_args = ''
params.rfd_compress_trajectories = true

params.pmpnn_relax_cycles = 3
params.pmpnn_seqs_per_struct = 1
params.pmpnn_weights = false
params.pmpnn_temperature = 0.000001
params.pmpnn_augment_eps = 0
params.pmpnn_omit_aas = 'CX'

params.max_rg = false
params.rfd_filters = false

params.af2ig_recycle = 3

params.refold_af2ig_filters = false
params.refold_max = false
params.refold_use_msa_server = false
params.refold_create_target_msa = false
params.refold_target_templates = false
params.refold_target_fasta = false
params.uniref30 = false
params.colabfold_envdb = false

params.output_rmsd_aligned = false

params.require_gpu = true
params.gpu_devices = ''
params.gpu_allocation_detect_process_regex = '(python.*/app/dl_binder_design/af2_initial_guess/predict\\.py|python.*/app/BindCraft/bindcraft\\.py|boltz predict|python.*/app/RFdiffusion/scripts/run_inference\\.py)'

params.rfd_backbone_models = false

include { RFDIFFUSION } from '../modules/local/rfd/rfdiffusion'
include { SILENT_FROM_PDBS } from '../modules/local/rfd/silentfrompdbs'
include { DL_BINDER_DESIGN_PROTEINMPNN } from '../modules/local/rfd/dl_binder_design'
include { AF2_INITIAL_GUESS } from '../modules/local/rfd/af2_initial_guess'
include { FILTER_DESIGNS } from '../modules/local/rfd/filter_designs'
include { AF2IG_SCORE_FILTER } from '../modules/local/rfd/af2ig_score_filter'

include { UNIQUE_ID } from '../modules/local/common/unique_id'

include { BOLTZ_REFOLD_SCORING } from '../subworkflows/local/boltz_refold_scoring'

workflow RFD {

    main:

    // Show help message
    if (params.input_pdb == false && params.rfd_backbone_models == false) {
        log.info(
            """
        ==================================================================
        PROTEIN BINDER DESIGN PIPELINE - RFDiffusion
        ==================================================================

        Required arguments:
            --input_pdb           Input PDB file for the target

            __or__

            --rfd_backbone_models Existing RFDiffusion backbone models - skips running RFDiffusion and uses these instead (glob accepted, eg 'results/rfdiffusion/pdbs/*.pdb)

            Optional arguments:
                --outdir              Output directory [default: ${params.outdir}]
                --design_name         Name of the design, used for output file prefixes [default: ${params.design_name}]
                --contigs             Contig map for RFdiffusion [default: ${params.contigs}]
                --hotspot_res         Hotspot residues, eg "A473,A995,A411,A421" - you must include the chain ID in every hotspot [default: ${params.hotspot_res}]
                --rfd_batch_size      Number of designs per batch [default: ${params.rfd_batch_size}]
                --rfd_model_path      Path to RFdiffusion model checkpoint file - leaving unset will allow RFDiffusion to choose based on other parameters [default: ${params.rfd_model_path}]
                --rfd_extra_args      Extra arguments for RFdiffusion [default: ${params.rfd_extra_args}]
                --rfd_config          'base', 'symmetry' or a path to a YAML file [default: ${params.rfd_config}]
                --rfd_compress_trajectories Compress trajectories with gzip [default: ${params.rfd_compress_trajectories}]
                --pmpnn_weights       Path to ProteinMPNN weights file (leave unset to use default weights) [default: ${params.pmpnn_weights}]
                --pmpnn_temperature   Temperature for ProteinMPNN [default: ${params.pmpnn_temperature}]
                --pmpnn_augment_eps   Variance of random noise to add to the atomic coordinates ProteinMPNN [default: ${params.pmpnn_augment_eps}]
                --pmpnn_relax_cycles  Number of relax cycles for ProteinMPNN [default: ${params.pmpnn_relax_cycles}]
                --pmpnn_seqs_per_struct Number of sequences per structure for ProteinMPNN [default: ${params.pmpnn_seqs_per_struct}]
                --pmpnn_omit_aas      A string of all residue types (one letter case-insensitive) that should not appear in the design [default: ${params.pmpnn_omit_aas}]
                --max_rg              Maximum radius of gyration for backbone filtering [default: disabled]
                --rfd_filters         Semicolon-separated list of filters for RFDiffusion backbones, eg "rg<25;compactness>0.8" [default: disabled]
                
                --refold_af2ig_filters       Semicolon-separated list of filters for AF2 initial guess designs, eg "pae_interaction<=10;plddt_binder>=80" [default: disabled]
                --af2ig_recycle       Number of recycle cycles for AF2 initial guess [default: ${params.af2ig_recycle}]
                    
                --refold_max          Maximum number of designs to refold with Boltz-2 [default: disabled]
                --refold_use_msa_server Use Boltz MSA server for target sequences [default: ${params.refold_use_msa_server}]
                --refold_create_target_msa Create MSA for target sequences [default: ${params.refold_create_target_msa}]
                --refold_target_templates Templates directory with .cif files for Boltz-2 [default: ${params.refold_target_templates}]
                --refold_target_fasta  FASTA file with full-length target sequences (headers should match PDB basenames) [default: ${params.refold_target_fasta}]
                --uniref30            UniRef30 database path for MSA creation [default: ${params.uniref30}]
                --colabfold_envdb     ColabFold environment database path for MSA creation [default: ${params.colabfold_envdb}]
                --output_rmsd_aligned Output aligned PDB files from RMSD calculations [default: ${params.output_rmsd_aligned}]
                    
                --require_gpu         Fail tasks that go too slow without a GPU if no GPU is detected [default: ${params.require_gpu}]
                --gpu_devices         GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
                --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]

        """.stripIndent()
        )
        exit(1)
    }

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
    }
    else {
        ch_rfd_config = Channel.value(false)
    }
    if (params.input_pdb) {
        ch_input_pdb = Channel.fromPath(params.input_pdb).first()
    }
    else {
        ch_input_pdb = Channel.value(false)
    }
    if (params.rfd_model_path) {
        ch_rfd_model_path = Channel.fromPath(params.rfd_model_path).first()
    }
    else {
        ch_rfd_model_path = Channel.value(false)
    }

    def hotspot_res = params.hotspot_res
    // We ensure hotspot_res has [square_brackets] by trimming any existing brackets and adding them back
    if (params.hotspot_res) {
        // Remove any leading/trailing whitespace, '[' and ']' characters, then add brackets back
        hotspot_res = params.hotspot_res.trim().replaceAll(/^\[+/, '').replaceAll(/\]+\$/, '')
        hotspot_res = "[${params.hotspot_res}]"
    }

    // Create channel of start numbers for batches
    ch_rfd_startnum = Channel.of(0..params.rfd_n_designs - 1)
        .filter { v -> v % params.rfd_batch_size == 0 }

    if (params.rfd_backbone_models) {
        ch_rfd_backbone_models = Channel.fromPath("${params.rfd_backbone_models}")
    }
    else {
        // Run RFdiffusion in batches
        RFDIFFUSION(
            ch_rfd_config,
            ch_input_pdb,
            ch_rfd_model_path,
            params.contigs,
            hotspot_res,
            params.rfd_batch_size,
            ch_rfd_startnum,
            ch_unique_id,
        )
        ch_rfd_backbone_models = RFDIFFUSION.out.pdbs.flatten()
    }

    if (params.rfd_filters) {
        FILTER_DESIGNS(
            ch_rfd_backbone_models,
            params.rfd_filters,
            'A',
            'rfdiffusion',
        )
        ch_filtered_backbones = FILTER_DESIGNS.out.accepted
    }
    else {
        ch_filtered_backbones = ch_rfd_backbone_models
    }

    // Create a channel that repeats each PDB params.pmpnn_seqs_per_struct times
    // and pairs it with an index from 0 to pmpnn_seqs_per_struct-1
    ch_pmpnn_inputs = ch_filtered_backbones
        | combine(Channel.of(0..(params.pmpnn_seqs_per_struct - 1)))

    // Run ProteinMPNN (dl_binder_design) on backbone-only PDBs
    DL_BINDER_DESIGN_PROTEINMPNN(
        ch_pmpnn_inputs.map { pdb, idx -> pdb },
        Channel.value(1),
        params.pmpnn_relax_cycles,
        params.pmpnn_weights,
        params.pmpnn_temperature,
        params.pmpnn_augment_eps,
        ch_pmpnn_inputs.map { pdb, idx -> idx },
    )

    // Run AF2 initial guess to build/refine sidechains for compatible sequence
    AF2_INITIAL_GUESS(
        DL_BINDER_DESIGN_PROTEINMPNN.out.pdbs
    )

    af2ig_scores = AF2_INITIAL_GUESS.out.pdbs_with_scores
        .map { pdbs, scores -> scores }
        .collectFile(
            name: 'af2ig_scores.tsv',
            storeDir: "${params.outdir}/af2_initial_guess",
            keepHeader: true,
            skip: 1,
        )

    // Run Boltz-2 refolding and scoring (or just AF2IG scoring if no refolding filters)
    BOLTZ_REFOLD_SCORING(
        AF2_INITIAL_GUESS.out.pdbs_with_scores,
        AF2_INITIAL_GUESS.out.pdbs,
        AF2_INITIAL_GUESS.out.scores,
        'A',  // binder_chain
        'B',  // target_chain
        params.refold_af2ig_filters,
        params.refold_max,
        params.refold_create_target_msa,
        params.refold_use_msa_server,
        params.refold_target_fasta,
        params.refold_target_templates,
        params.colabfold_envdb,
        params.uniref30,
        params.outdir,
    )

    emit:
    scores = BOLTZ_REFOLD_SCORING.out.scores
}
