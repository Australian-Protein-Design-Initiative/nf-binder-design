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
  --rfd_partial_T=5,10,20,50 \
  -resume \
  -with-report report_$(date +%Y%m%d_%H%M%S).html \
  -with-trace trace_$(date +%Y%m%d_%H%M%S).txt

*/

// Default parameters
params.input_pdb = false
params.outdir = 'results'

params.design_name = 'fuzzed_ppi'
params.binder_chain = 'A' // eg, our fixed target chain - usually A ?
params.target_contigs = 'auto' // 'auto' to detect from PDB, or eg "B10-110" for fixed target chain B, residues 10-110
params.hotspot_res = false
params.rfd_batch_size = 10
params.rfd_n_partial_per_binder = 10
params.rfd_model_path = false // "models/rfdiffusion/Complex_beta_ckpt.pt"
params.rfd_config = 'base'
// noise_scale_ca and noise_scale_frame
params.rfd_noise_scale = 0
params.rfd_partial_T = 20 // Can be a single value or comma-separated list like "5,10,20"
params.rfd_extra_args = ''
params.skip_renumber = false
params.rfd_compress_trajectories = true

params.pmpnn_relax_cycles = 3  // or 5, or maybe even more
params.pmpnn_seqs_per_struct = 1
params.pmpnn_weights = false
params.pmpnn_temperature = 0.000001
params.pmpnn_augment_eps = 0
params.pmpnn_omit_aas = 'CX'

params.max_rg = false
params.rfd_filters = false

params.refold_af2ig_filters = false
params.af2ig_recycle = 3

params.refold_max = false
params.refold_use_msa_server = false
params.refold_create_target_msa = false
params.refold_target_templates = false
params.refold_target_fasta = false
params.uniref30 = false
params.colabfold_envdb = false

params.require_gpu = true
params.gpu_devices = ''
params.gpu_allocation_detect_process_regex = '(python.*/app/dl_binder_design/af2_initial_guess/predict\\.py|python.*/app/BindCraft/bindcraft\\.py|boltz predict|python.*/app/RFdiffusion/scripts/run_inference\\.py)'

include { RFDIFFUSION } from './modules/rfdiffusion'
include { RFDIFFUSION_PARTIAL } from './modules/rfdiffusion_partial'
include { SILENT_FROM_PDBS } from './modules/silentfrompdbs'
include { DL_BINDER_DESIGN_PROTEINMPNN } from './modules/dl_binder_design'
include { AF2_INITIAL_GUESS } from './modules/af2_initial_guess'
include { GET_CONTIGS } from './modules/get_contigs'
include { RENUMBER_RESIDUES } from './modules/renumber_residues'
include { COMBINE_SCORES } from './modules/combine_scores'
include { UNIQUE_ID } from './modules/unique_id'
include { FILTER_DESIGNS } from './modules/filter_designs'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_AF2IG } from './modules/bindcraft_scoring.nf'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_BOLTZ_COMPLEX } from './modules/bindcraft_scoring.nf'
include { AF2IG_SCORE_FILTER } from './modules/af2ig_score_filter'
include { BOLTZ_COMPARE_COMPLEX } from './modules/boltz_compare_complex'
include { BOLTZ_COMPARE_BINDER_MONOMER } from './modules/boltz_compare_binder_monomer'
include { MMSEQS_COLABFOLDSEARCH } from './modules/mmseqs_colabfoldsearch'
include { PDB_TO_FASTA } from './modules/pdb_to_fasta'

// Validate numeric parameters
def validate_numeric = { param_name, value ->
    if (!(value instanceof Number)) {
        error "Parameter $param_name must be a number, got: $value (${value.getClass().getName()})"
    }
}

validate_numeric('rfd_batch_size', params.rfd_batch_size)
validate_numeric('rfd_n_partial_per_binder', params.rfd_n_partial_per_binder)
validate_numeric('rfd_noise_scale', params.rfd_noise_scale)
// Special handling for rfd_partial_T - can be a single value or comma-separated list
def validate_rfd_partial_T = { value ->
    if (value instanceof Number) {
        return [value] // Return as a list with one element
    } else if (value instanceof String) {
        try {
            // Try to parse as comma-separated list of numbers
            return value.split(',').collect { it.trim().toInteger() }
        } catch (Exception e) {
            error "Parameter rfd_partial_T must be a number or comma-separated list of numbers, got: $value"
        }
    } else {
        error "Parameter rfd_partial_T must be a number or comma-separated list of numbers, got: $value (${value.getClass().getName()})"
    }
}
def partial_T_values = validate_rfd_partial_T(params.rfd_partial_T)

validate_numeric('pmpnn_relax_cycles', params.pmpnn_relax_cycles)
validate_numeric('pmpnn_seqs_per_struct', params.pmpnn_seqs_per_struct)

workflow {
    // Show help message
    if (params.input_pdb == false) {
        log.info"""
        ==================================================================
        🧬 PROTEIN BINDER DESIGN PIPELINE 🧬
        ==================================================================

        Required arguments:
            --input_pdb           Input PDBs file for the binders to diffuse (* glob accepted)

        Optional arguments:
            --outdir              Output directory [default: ${params.outdir}]
            --design_name         Name of the design, used for output file prefixes [default: ${params.design_name}]
            --target_contigs      Contig map for target chain(s) - 'auto' to detect from PDB, or specify manually [default: ${params.target_contigs}]
            --binder_chain        Chain ID of the binder chain [default: ${params.binder_chain}]
            --hotspot_res         Hotspot residues, eg "A473,A995,A411,A421" - you must include the chain ID in every hotspot [default: ${params.hotspot_res}]
            --skip_renumber       Skip the residue renumbering step [default: ${params.skip_renumber}]
            --rfd_batch_size      Number of designs per batch [default: ${params.rfd_batch_size}]
            --rfd_n_partial_per_binder Number of partial diffused designs per binder [default: ${params.rfd_n_partial_per_binder}]
            --rfd_model_path      Path to RFdiffusion model checkpoint file - you probaby don't want to set this manually [default: ${params.rfd_model_path}]
            --rfd_extra_args      Extra arguments for RFdiffusion [default: ${params.rfd_extra_args}]
            --rfd_config          'base', 'symmetry' or a path to a YAML file [default: ${params.rfd_config_name}]
            --rfd_partial_T       Number of timesteps to run partial diffusion for (lower = less diffusion)
                                Can be a single value or comma-separated list like "5,10,20,50" [default: ${params.rfd_partial_T}]
            --rfd_compress_trajectories Compress trajectories with gzip [default: ${params.rfd_compress_trajectories}]

            --pmpnn_relax_cycles  Number of relax cycles for ProteinMPNN [default: ${params.pmpnn_relax_cycles}]
            --pmpnn_seqs_per_struct Number of sequences per structure for ProteinMPNN [default: ${params.pmpnn_seqs_per_struct}]
            --pmpnn_weights       Path to ProteinMPNN weights file (leave unset to use default weights) [default: ${params.pmpnn_weights}]
            --pmpnn_temperature   Temperature for ProteinMPNN [default: ${params.pmpnn_temperature}]
            --pmpnn_augment_eps   Variance of random noise to add to the atomic coordinates ProteinMPNN [default: ${params.pmpnn_augment_eps}]
            --pmpnn_omit_aas      A string of all residue types (one letter case-insensitive) that should not appear in the design [default: ${params.pmpnn_omit_aas}]
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
            
            --require_gpu         Fail tasks that go too slow without a GPU if no GPU is detected [default: ${params.require_gpu}]
            --gpu_devices         GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
            --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]

        """.stripIndent()
        exit 1
    }

    // Generate unique ID for this run
    UNIQUE_ID()
    ch_unique_id = UNIQUE_ID.out.id_file.map { it.text.trim() }

    if (!params.input_pdb) {
        throw new Exception('--input_pdb must be provided')
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

    def hotspot_res = params.hotspot_res
    // We ensure hotspot_res has [square_brackets] by trimming any existing brackets and adding them back
    if (params.hotspot_res) {
        // Remove any leading/trailing whitespace, '[' and ']' characters, then add brackets back
        hotspot_res = params.hotspot_res.trim().replaceAll(/^\[+/, '').replaceAll(/\]+\$/, '')
        hotspot_res = "[${params.hotspot_res}]"
    }

    // We renumber the residues in each chain to be 1 to chain_length
    // irrespective of missing non-sequential residue numbers etc
    // The --binder_chain is always set to chain A, with the non-diffusable chains
    // coming after it.
    if (params.skip_renumber) {
        // Skip renumbering and use input PDBs directly
        ch_preprocessed_pdb = ch_input_pdb
        log.info 'Skipping residue renumbering as --skip-renumber was set'
    } else {
        // Warn about renumbering when using hotspot residues
        if (params.hotspot_res) {
            log.warn 'WARNING: Target residues will be renumbered starting at 1 - do your chosen hotspots account for this?'
        } else {
            log.warn 'WARNING: No hotspots defined - binders will tend to drift from target'
        }

        // Apply renumbering as normal
        ch_preprocessed_pdb = ch_input_pdb
            .map { pdb -> [pdb, params.binder_chain] }
            | RENUMBER_RESIDUES
    }

    // Get contigs string for each input PDB
    // NOTE/HACK: We hardcode 'A' here, since RENUMBER_RESIDUES always
    // makes the binder chain 'A'
    def binder_chain_for_contigs = params.skip_renumber ? params.binder_chain : 'A'
    ch_contigs = ch_preprocessed_pdb
        .map { pdb -> [pdb, binder_chain_for_contigs] }
        | GET_CONTIGS

    //ch_contigs.view()

    // Create batches with contigs
    ch_contigs
        | map { input_pdb, contigs ->
            def num_batches = (params.rfd_n_partial_per_binder / params.rfd_batch_size).toInteger()
            // Create jobs for each combination of batch and partial_T
            def all_jobs = []
            partial_T_values.each { partial_T ->
                (0..<num_batches).each { batchNum ->
                    all_jobs << tuple(input_pdb, contigs, batchNum * params.rfd_batch_size, partial_T)
                }
            }
            return all_jobs
        }
        | flatMap()  // Flatten the lists of tuples into individual tuples
        | set { ch_rfd_jobs }

    // Run RFdiffusion with partial diffusion in batches
    RFDIFFUSION_PARTIAL(
        ch_rfd_config,
        ch_rfd_jobs.map { input_pdb, contigs, start, partial_T -> input_pdb },  // Extract PDB path
        ch_rfd_model_path,
        ch_rfd_jobs.map { input_pdb, contigs, start, partial_T -> contigs },  // Extract contigs string
        hotspot_res,
        params.rfd_batch_size,
        ch_rfd_jobs.map { input_pdb, contigs, start, partial_T -> start },  // Extract start number
        ch_unique_id,
        ch_rfd_jobs.map { input_pdb, contigs, start, partial_T -> partial_T }  // Extract partial_T
    )
    ch_rfd_backbone_models = RFDIFFUSION_PARTIAL.out.pdbs.flatten()

    if (params.rfd_filters) {
        FILTER_DESIGNS(
            ch_rfd_backbone_models,
            params.rfd_filters,
            binder_chain_for_contigs,
            'rfdiffusion'  // step name
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
        DL_BINDER_DESIGN_PROTEINMPNN.out.pdbs
    )

    af2ig_scores = AF2_INITIAL_GUESS.out.scores_tsv.collectFile(
        name: 'af2ig_scores.tsv',
        storeDir: "${params.outdir}/af2_initial_guess",
        keepHeader: true,
        skip: 1)

    // Optional Boltz-2 refolding of filtered designs
    if (params.refold_af2ig_filters) {
        // Filter designs by score thresholds
        AF2IG_SCORE_FILTER(
            AF2_INITIAL_GUESS.out.scores_tsv,
            AF2_INITIAL_GUESS.out.pdbs.collect(),
            params.refold_af2ig_filters
        )
        
        // Flatten and create metadata for each accepted PDB
        ch_filtered_with_meta = AF2IG_SCORE_FILTER.out.accepted
            .flatten()
            .map { pdb -> [[id: pdb.baseName], pdb] }
        
        // Limit to refold_max if specified
        ch_filtered_for_refold = params.refold_max ? 
            ch_filtered_with_meta.take(params.refold_max) : 
            ch_filtered_with_meta
        
        // Determine target chain (all non-binder chains, typically 'B' in partial.nf context)
        def target_chain = 'B'
        
        // Optionally create target MSAs
        if (params.refold_create_target_msa && !params.refold_use_msa_server) {
            if (params.refold_target_fasta) {
                // Use the refold_target_fasta file directly for MSA creation
                ch_target_fastas = ch_filtered_for_refold.map { meta, pdb ->
                    [meta, file(params.refold_target_fasta)]
                }
                ch_target_msas = MMSEQS_COLABFOLDSEARCH(
                    ch_target_fastas,
                    params.colabfold_envdb,
                    params.uniref30
                )
            } else {
                // Create target FASTA files from PDBs for MSA creation
                ch_target_fastas = ch_filtered_for_refold
                    .map { meta, pdb -> pdb }
                    | PDB_TO_FASTA(
                        target_chain  // target chain
                    )
                    .map { fasta -> 
                        def basename = fasta.baseName.replaceAll(/_B$/, '')
                        [[id: basename], fasta]
                    }
                ch_target_msas = MMSEQS_COLABFOLDSEARCH(
                    ch_target_fastas,
                    params.colabfold_envdb,
                    params.uniref30
                )
            }
        } else if (params.refold_create_target_msa && params.refold_use_msa_server) {
            ch_target_msas = ch_filtered_for_refold.map { meta, pdb -> 
                [meta, file("${projectDir}/assets/dummy_files/boltz_will_make_target_msa")] 
            }
        } else {
            ch_target_msas = ch_filtered_for_refold.map { meta, pdb -> 
                [meta, file("${projectDir}/assets/dummy_files/empty_target_msa")] 
            }
        }
        
        // Run Boltz complex refolding with RMSD analysis
        BOLTZ_COMPARE_COMPLEX(
            ch_filtered_for_refold,
            binder_chain_for_contigs,  // binder_chain
            target_chain,              // target_chain
            params.refold_create_target_msa,
            params.refold_use_msa_server,
            ch_target_msas.map { meta, target_msa -> target_msa },
            file("${projectDir}/assets/dummy_files/empty_binder_msa"),
            file(params.refold_target_templates ?: "${projectDir}/assets/dummy_files/empty_templates"),
            params.refold_target_fasta ? file(params.refold_target_fasta): "${projectDir}/assets/dummy_files/empty"
        )
        
        // Run Boltz binder monomer prediction with RMSD analysis
        BOLTZ_COMPARE_BINDER_MONOMER(
            ch_filtered_for_refold
                .join(BOLTZ_COMPARE_COMPLEX.out.pdb)
                .map { meta, af2ig_pdb, boltz_pdb ->
                    [meta, af2ig_pdb, boltz_pdb]
                },
            binder_chain_for_contigs  // binder_chain
        )
        
        // Aggregate RMSD outputs
        ch_target_aligned_rmsd = BOLTZ_COMPARE_COMPLEX.out.rmsd_target_aligned
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_target_aligned_binder.tsv',
                storeDir: "${params.outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1
            )
        
        ch_complex_rmsd = BOLTZ_COMPARE_COMPLEX.out.rmsd_complex
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_complex_vs_af2ig.tsv',
                storeDir: "${params.outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1
            )
        
        ch_monomer_vs_af2ig_rmsd = BOLTZ_COMPARE_BINDER_MONOMER.out.rmsd_monomer_vs_af2ig
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_monomer_vs_af2ig.tsv',
                storeDir: "${params.outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1
            )
        
        ch_monomer_vs_complex_rmsd = BOLTZ_COMPARE_BINDER_MONOMER.out.rmsd_monomer_vs_complex
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_monomer_vs_complex.tsv',
                storeDir: "${params.outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1
            )
        
        // Aggregate confidence outputs
        ch_complex_confidence = BOLTZ_COMPARE_COMPLEX.out.confidence_tsv
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'boltz_scores_complex.tsv',
                storeDir: "${params.outdir}/boltz_refold",
                keepHeader: true,
                skip: 1
            )
        
        ch_monomer_confidence = BOLTZ_COMPARE_BINDER_MONOMER.out.confidence_tsv
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'boltz_scores_binder_monomer.tsv',
                storeDir: "${params.outdir}/boltz_refold",
                keepHeader: true,
                skip: 1
            )

        BINDCRAFT_SCORING_BOLTZ_COMPLEX(
            BOLTZ_COMPARE_COMPLEX.out.pdb.map { meta, pdb -> pdb },
            binder_chain_for_contigs,
            'default_4stage_multimer'
        )

        extra_scores = BINDCRAFT_SCORING_BOLTZ_COMPLEX.out.scores.collectFile(
            name: 'boltz_complex_extra_scores.tsv',
            storeDir: "${params.outdir}/boltz_refold",
            keepHeader: true,
            skip: 1)

        // Combine all the score files into a single TSV file
        COMBINE_SCORES(
            AF2_INITIAL_GUESS.out.scores.collect(),
            extra_scores,
            ch_complex_confidence,
            ch_monomer_vs_complex_rmsd,
            AF2_INITIAL_GUESS.out.pdbs.collect()
        )
    } else {

        BINDCRAFT_SCORING_AF2IG(
            AF2_INITIAL_GUESS.out.pdbs,
            binder_chain_for_contigs,
            'default_4stage_multimer'
        )

        extra_scores = BINDCRAFT_SCORING_AF2IG.out.scores.collectFile(
            name: "${params.outdir}/af2_initial_guess/af2ig_extra_scores.tsv",
            //storeDir: "${params.outdir}",
            keepHeader: true,
            skip: 1)

        // Combine all the score files into a single TSV file
        COMBINE_SCORES(
            AF2_INITIAL_GUESS.out.scores.collect(),
            extra_scores,
            "${projectDir}/assets/dummy_files/empty",
            "${projectDir}/assets/dummy_files/empty",
            AF2_INITIAL_GUESS.out.pdbs.collect()
        )
    }
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