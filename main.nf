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
params.contigs = ''
// "[A371-508/A753-883/A946-1118/A1135-1153/0 70-100]"
params.hotspot_res = false
// "A473,A995,A411,A421"
params.rfd_batch_size = 1
params.rfd_n_designs = 2
params.rfd_model_path = false
// "models/rfdiffusion/Complex_beta_ckpt.pt"
params.rfd_config = false
// "base"
// noise_scale_ca and noise_scale_frame
params.rfd_noise_scale = 0
params.rfd_extra_args = ''
params.rfd_compress_trajectories = true

params.pmpnn_relax_cycles = 3
// or 5, or maybe even more
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

// Set to a path of existing RFDiffusion backbone models, skips running RFDiffusion
params.rfd_backbone_models = false

include { RFDIFFUSION } from './modules/rfdiffusion'
include { SILENT_FROM_PDBS } from './modules/silentfrompdbs'
include { DL_BINDER_DESIGN_PROTEINMPNN } from './modules/dl_binder_design'
include { AF2_INITIAL_GUESS } from './modules/af2_initial_guess'
include { COMBINE_SCORES } from './modules/combine_scores'
include { UNIQUE_ID } from './modules/unique_id'
include { FILTER_DESIGNS } from './modules/filter_designs'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_AF2IG } from './modules/bindcraft_scoring'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_BOLTZ_COMPLEX } from './modules/bindcraft_scoring'
include { AF2IG_SCORE_FILTER } from './modules/af2ig_score_filter'
include { PDB_TO_FASTA } from './modules/pdb_to_fasta'
include { BOLTZ_COMPARE_COMPLEX } from './modules/boltz_compare_complex'
include { BOLTZ_COMPARE_BINDER_MONOMER } from './modules/boltz_compare_binder_monomer'
include { MMSEQS_COLABFOLDSEARCH } from './modules/mmseqs_colabfoldsearch'
include { COMPRESS_GPU_STATS; paramsToMap } from './modules/utils.nf'

workflow {

    // Show help message
    if (params.input_pdb == false && params.rfd_backbone_models == false) {
        log.info(
            """
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
                --hotspot_res         Hotspot residues, eg "A473,A995,A411,A421" - you must include the chain ID in every hotspot [default: ${params.hotspot_res}]
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
                --enable_gpu_stats    Enable GPU utilisation monitoring [default: ${params.enable_gpu_stats}]
                --gpu_stats_interval GPU monitoring sampling interval in seconds [default: ${params.gpu_stats_interval}]

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
    //.view()

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

    // Optional Boltz-2 refolding of filtered designs
    if (params.refold_af2ig_filters) {
        // Filter designs by score thresholds
        AF2IG_SCORE_FILTER(
            AF2_INITIAL_GUESS.out.pdbs_with_scores.map { pdbs, scores -> scores },
            AF2_INITIAL_GUESS.out.pdbs_with_scores.map { pdbs, scores -> pdbs },
            params.refold_af2ig_filters,
        )

        // Flatten and create metadata for each accepted PDB
        ch_filtered_with_meta = AF2IG_SCORE_FILTER.out.accepted
            .flatten()
            .map { pdb -> [[id: pdb.baseName], pdb] }

        // Limit to refold_max if specified
        ch_filtered_for_refold = params.refold_max
            ? ch_filtered_with_meta.take(params.refold_max)
            : ch_filtered_with_meta

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
                    params.uniref30,
                )
            }
            else {
                // Create target FASTA files from PDBs for MSA creation
                ch_target_fastas = ch_filtered_for_refold.map { meta, pdb -> pdb }
                    | PDB_TO_FASTA(
                        'B'
                    ).map { fasta ->
                        def basename = fasta.baseName.replaceAll(/_B$/, '')
                        [[id: basename], fasta]
                    }
                ch_target_msas = MMSEQS_COLABFOLDSEARCH(
                    ch_target_fastas,
                    params.colabfold_envdb,
                    params.uniref30,
                )
            }
        }
        else if (params.refold_create_target_msa && params.refold_use_msa_server) {
            ch_target_msas = ch_filtered_for_refold.map { meta, pdb ->
                [meta, file("${projectDir}/assets/dummy_files/boltz_will_make_target_msa")]
            }
        }
        else {
            ch_target_msas = ch_filtered_for_refold.map { meta, pdb ->
                [meta, file("${projectDir}/assets/dummy_files/empty_target_msa")]
            }
        }

        // Run Boltz complex refolding with RMSD analysis
        BOLTZ_COMPARE_COMPLEX(
            ch_filtered_for_refold,
            'A',
            'B',
            params.refold_create_target_msa,
            params.refold_use_msa_server,
            ch_target_msas.map { meta, target_msa -> target_msa },
            file("${projectDir}/assets/dummy_files/empty_binder_msa"),
            file(params.refold_target_templates ?: "${projectDir}/assets/dummy_files/empty_templates"),
            params.refold_target_fasta ? file(params.refold_target_fasta) : "${projectDir}/assets/dummy_files/empty",
        )

        // Run Boltz binder monomer prediction with RMSD analysis
        BOLTZ_COMPARE_BINDER_MONOMER(
            ch_filtered_for_refold.join(BOLTZ_COMPARE_COMPLEX.out.pdb).map { meta, af2ig_pdb, boltz_pdb ->
                [meta, af2ig_pdb, boltz_pdb]
            },
            'A',
        )

        // Aggregate RMSD outputs
        ch_target_aligned_rmsd = BOLTZ_COMPARE_COMPLEX.out.rmsd_target_aligned
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_target_aligned_binder.tsv',
                storeDir: "${params.outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1,
            )

        ch_complex_rmsd = BOLTZ_COMPARE_COMPLEX.out.rmsd_complex
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_complex_vs_af2ig.tsv',
                storeDir: "${params.outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1,
            )

        ch_monomer_vs_af2ig_rmsd = BOLTZ_COMPARE_BINDER_MONOMER.out.rmsd_monomer_vs_af2ig
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_monomer_vs_af2ig.tsv',
                storeDir: "${params.outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1,
            )

        ch_monomer_vs_complex_rmsd = BOLTZ_COMPARE_BINDER_MONOMER.out.rmsd_monomer_vs_complex
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_monomer_vs_complex.tsv',
                storeDir: "${params.outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1,
            )

        // Aggregate confidence outputs
        ch_complex_confidence = BOLTZ_COMPARE_COMPLEX.out.confidence_tsv
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'boltz_scores_complex.tsv',
                storeDir: "${params.outdir}/boltz_refold",
                keepHeader: true,
                skip: 1,
            )

        ch_monomer_confidence = BOLTZ_COMPARE_BINDER_MONOMER.out.confidence_tsv
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'boltz_scores_binder_monomer.tsv',
                storeDir: "${params.outdir}/boltz_refold",
                keepHeader: true,
                skip: 1,
            )

        // We BindCraft-score just the Boltz-2 refolded complexes 
        // (rather than the AF2 initial guess pdbs)
        BINDCRAFT_SCORING_BOLTZ_COMPLEX(
            BOLTZ_COMPARE_COMPLEX.out.pdb.map { meta, pdb -> pdb },
            'A',
            'default_4stage_multimer',
        )

        extra_scores = BINDCRAFT_SCORING_BOLTZ_COMPLEX.out.scores.collectFile(
            name: 'boltz_complex_extra_scores.tsv',
            storeDir: "${params.outdir}/boltz_refold",
            keepHeader: true,
            skip: 1,
        )

        // Combine all the score files into a single TSV file
        COMBINE_SCORES(
            AF2_INITIAL_GUESS.out.scores.collect(),
            extra_scores.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            ch_complex_confidence.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            ch_monomer_vs_complex_rmsd.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            ch_target_aligned_rmsd.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            AF2_INITIAL_GUESS.out.pdbs.collect(),
        )

        // Merge GPU stats including Boltz processes
        ch_gpu_stats_all = DL_BINDER_DESIGN_PROTEINMPNN.out.gpu_stats
            .mix(AF2_INITIAL_GUESS.out.gpu_stats)
            .mix(BOLTZ_COMPARE_COMPLEX.out.gpu_stats)
            .mix(BOLTZ_COMPARE_BINDER_MONOMER.out.gpu_stats)
        if (!params.rfd_backbone_models) {
            ch_gpu_stats_all = ch_gpu_stats_all.mix(RFDIFFUSION.out.gpu_stats)
        }
        ch_gpu_stats_merged = ch_gpu_stats_all
            .collectFile(name: 'gpu_stats.csv', keepHeader: true, skip: 1)
        
        COMPRESS_GPU_STATS(ch_gpu_stats_merged)
    }
    else {

        BINDCRAFT_SCORING_AF2IG(
            AF2_INITIAL_GUESS.out.pdbs,
            'A',
            'default_4stage_multimer',
        )

        extra_scores = BINDCRAFT_SCORING_AF2IG.out.scores.collectFile(
            name: "${params.outdir}/af2_initial_guess/af2ig_extra_scores.tsv",
            keepHeader: true,
            skip: 1,
        )

        // Combine all the score files into a single TSV file
        COMBINE_SCORES(
            AF2_INITIAL_GUESS.out.scores.collect(),
            extra_scores,
            "${projectDir}/assets/dummy_files/empty",
            "${projectDir}/assets/dummy_files/empty",
            "${projectDir}/assets/dummy_files/empty",
            AF2_INITIAL_GUESS.out.pdbs.collect(),
        )

        // Merge GPU stats (no Boltz processes in this branch)
        ch_gpu_stats_all = DL_BINDER_DESIGN_PROTEINMPNN.out.gpu_stats
            .mix(AF2_INITIAL_GUESS.out.gpu_stats)
        if (!params.rfd_backbone_models) {
            ch_gpu_stats_all = ch_gpu_stats_all.mix(RFDIFFUSION.out.gpu_stats)
        }
        ch_gpu_stats_merged = ch_gpu_stats_all
            .collectFile(name: 'gpu_stats.csv', keepHeader: true, skip: 1)
        
        COMPRESS_GPU_STATS(ch_gpu_stats_merged)
    }

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
            success: workflow.success,
        ]

        def output_file = "${params.outdir}/params.json"
        def json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params_json))

        new File(params.outdir).mkdirs()
        new File(output_file).text = json_string

        log.info("Pipeline parameters saved to: ${output_file}")
    }
}
