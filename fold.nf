#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Standalone multi-method structure folding dispatcher: predicts structures for
FASTA inputs with any combination of --methods af2,boltz,rf3,protenix, sharing
a single MSA-generation stage (FOLD_MSA) with a selectable --msa_method. This
deliberately overlaps workflows/boltz_pulldown.nf; it is a standalone
prototype for now and is expected to be folded into a renamed/expanded
boltz_pulldown inside main.nf later (not in this work). See
plans/fold-nf-multi-method-folding.md for the full design.

*** Phase 1 of that plan: MONOMER INPUTS ONLY. *** A multi-record FASTA
(multimer) fails fast with a clear error - see §4 of the plan for the
deferred paired-MSA multimer strategy (Phase 2). Protenix (Phase 3) is
monomer-only for the same reason.

Usage:
  nextflow run fold.nf --input 'input/*.fasta' --outdir results \
      --methods af2,boltz,rf3,protenix --msa_method jackhmmer_af2 -profile slurm,m3
*/

// Default parameters
params.help = false
params.input = false
params.outdir = 'results'
params.methods = 'af2'
params.msa_method = 'jackhmmer_af2'

// --- AF2 (monomer-verified route: jackhmmer MSA -> ALPHAFOLD2 predict) ---
params.alphafold2_db_path = '/mnt/datasets/alphafold/alphafold_20240229'
params.alphafold2_model_preset = 'monomer_ptm' // monomer|monomer_ptm|monomer_casp14|multimer
params.alphafold2_db_preset = 'full_dbs'
params.alphafold2_max_template_date = '2024-01-01'
params.alphafold2_random_seed = false
params.alphafold2_num_predictions_per_model = 1 // multimer only
params.alphafold2_models_to_relax = 'best' // all|best|none

// --- Boltz ---
params.use_msa_server = false // Boltz's own MMseqs2 server, independent of --msa_method
params.templates = false
params.boltz_recycling = false // --recycling_steps override (default: boltz's own default, currently 3)
params.boltz_diffusion_samples = false // --diffusion_samples override (default: 1)
params.boltz_sampling_steps = false // --sampling_steps override (default: 200)

// --- RF3 (same defaults as workflows/rfd3.nf) ---
params.rf3_ckpt_path = '/models/foundry/rf3_foundry_01_24_latest_remapped.ckpt'
params.rf3_num_steps = 50 // default for rf3 cli is 200, however 50 is faster with no difference in quality
params.rf3_n_recycles = 10
params.rf3_diffusion_batch_size = 5
params.rf3_early_stopping_plddt_threshold = 0.5 // exits early if mean pLDDT < 0.5 after the first recycle

// --- Protenix (defaults confirmed against `protenix pred --help` in
// protenix:v2.0.0-weights on 2026-07-17 - identical to that build's own CLI
// defaults) ---
params.protenix_seeds = '101'
params.protenix_cycle = 10 // pairformer cycles (protenix's own --cycle default)
params.protenix_step = 200 // diffusion steps (protenix's own --step default)
params.protenix_sample = 5 // number of diffusion samples (protenix's own --sample default)
params.protenix_model_name = 'protenix_base_default_v1.0.0' // checkpoint baked into the container at /models/protenix/checkpoint
params.protenix_use_msa = true // false forces protenix to ignore our a3m and predict MSA-free

// --- ColabFold MSA (--msa_method mmseqs2_colabfold) ---
params.use_remote_server = false // query the ColabFold MMseqs2 API instead of a local DB search
params.uniref30 = false
params.colabfold_envdb = false

// require_gpu / gpu_devices / gpu_allocation_detect_process_regex are already
// defined with pipeline-wide defaults in nextflow.config

include { FOLD_MSA } from './subworkflows/local/fold_msa'
include { ALPHAFOLD2 } from './subworkflows/local/alphafold2'
include { BOLTZ_FOLD } from './subworkflows/local/boltz_fold'
include { ROSETTAFOLD3_FOLD } from './subworkflows/local/rosettafold3_fold'
include { PROTENIX_FOLD } from './subworkflows/local/protenix_fold'

def VALID_METHODS = ['af2', 'boltz', 'rf3', 'protenix']
def VALID_MSA_METHODS = ['jackhmmer_af2', 'mmseqs2_colabfold']

workflow {
    if (params.help || params.input == false) {
        log.info(
            """
        ==================================================================
        🧬 FOLD WORKFLOW (multi-method structure prediction) 🧬
        ==================================================================

        Predict structures for one or more FASTA files with any combination of
        AlphaFold2, Boltz-2, RosettaFold3 and Protenix, sharing one
        MSA-generation stage.

        *** Phase 1: MONOMER INPUTS ONLY. *** A FASTA file with more than one
        record (multimer) will fail fast - multimer support is coming in a
        later phase (see plans/fold-nf-multi-method-folding.md §4).

        Required arguments:
            --input                            Single FASTA file, glob, or directory of FASTA files.
                                                Each file is one prediction unit (must be single-record
                                                in this phase).

        Optional arguments:
            --outdir                           Output directory [default: ${params.outdir}]
            --methods                          Comma-separated list of af2,boltz,rf3,protenix [default: ${params.methods}]
            --msa_method                       jackhmmer_af2|mmseqs2_colabfold [default: ${params.msa_method}]

            AlphaFold2 (--methods includes af2):
            --alphafold2_db_path                AlphaFold2 database directory [default: ${params.alphafold2_db_path}]
            --alphafold2_model_preset            monomer|monomer_ptm|monomer_casp14|multimer [default: ${params.alphafold2_model_preset}]
            --alphafold2_db_preset               full_dbs|reduced_dbs [default: ${params.alphafold2_db_preset}]
            --alphafold2_max_template_date        Maximum template release date [default: ${params.alphafold2_max_template_date}]
            --alphafold2_random_seed              Fix the data pipeline's random seed [default: unset, randomly generated]
            --alphafold2_models_to_relax          all|best|none [default: ${params.alphafold2_models_to_relax}]

            Boltz-2 (--methods includes boltz):
            --use_msa_server                   Use Boltz's own MMseqs2 MSA server instead of the shared --msa_method a3m [default: ${params.use_msa_server}]
            --templates                        Templates directory with .cif files [default: ${params.templates}]
            --boltz_recycling                  Boltz --recycling_steps override [default: boltz's own default]
            --boltz_diffusion_samples          Boltz --diffusion_samples override [default: boltz's own default]
            --boltz_sampling_steps              Boltz --sampling_steps override [default: boltz's own default]

            RosettaFold3 (--methods includes rf3):
            --rf3_ckpt_path                     RF3 checkpoint path [default: ${params.rf3_ckpt_path}]
            --rf3_num_steps                     [default: ${params.rf3_num_steps}]
            --rf3_n_recycles                    [default: ${params.rf3_n_recycles}]
            --rf3_diffusion_batch_size           [default: ${params.rf3_diffusion_batch_size}]
            --rf3_early_stopping_plddt_threshold [default: ${params.rf3_early_stopping_plddt_threshold}]

            Protenix (--methods includes protenix):
            --protenix_seeds                    Comma-separated seeds [default: ${params.protenix_seeds}]
            --protenix_cycle                    Pairformer cycles [default: ${params.protenix_cycle}]
            --protenix_step                     Diffusion steps [default: ${params.protenix_step}]
            --protenix_sample                   Number of diffusion samples [default: ${params.protenix_sample}]
            --protenix_model_name               Checkpoint name (baked into the container) [default: ${params.protenix_model_name}]
            --protenix_use_msa                  Feed our shared a3m to Protenix; false predicts MSA-free [default: ${params.protenix_use_msa}]

            ColabFold MSA (--msa_method mmseqs2_colabfold):
            --use_remote_server                 Query the ColabFold MMseqs2 API instead of a local DB search [default: ${params.use_remote_server}]
            --uniref30                          UniRef30 database path (local search only) [default: ${params.uniref30}]
            --colabfold_envdb                   ColabFold environment database path (local search only) [default: ${params.colabfold_envdb}]
                                                Note: no M3 default path exists for these local ColabFold DBs
                                                (see plans/fold-nf-multi-method-folding.md); use --use_remote_server
                                                true if you haven't set them up yourself.

            --gpu_devices                        GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
            --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]

        Example:
            nextflow run fold.nf --input 'input/*.fasta' --outdir results \\
                --methods af2,boltz,rf3,protenix --msa_method jackhmmer_af2 -profile slurm,m3

        """.stripIndent()
        )
        exit(1)
    }

    def methods = params.methods.toString().split(',').collect { it.trim().toLowerCase() }
    methods.each { m ->
        if (!(m in VALID_METHODS)) {
            error("fold.nf: unknown --methods entry '${m}' (valid: ${VALID_METHODS.join(', ')})")
        }
    }
    if (!(params.msa_method in VALID_MSA_METHODS)) {
        error("fold.nf: unknown --msa_method '${params.msa_method}' (valid: ${VALID_MSA_METHODS.join(', ')})")
    }
    if (params.msa_method == 'mmseqs2_colabfold' && !params.use_remote_server && !(params.uniref30 && params.colabfold_envdb)) {
        error(
            "fold.nf: --msa_method mmseqs2_colabfold needs --uniref30 and --colabfold_envdb " +
            "(local ColabFold DBs) or --use_remote_server true. No M3 default local ColabFold " +
            "DB path exists yet (see plans/fold-nf-multi-method-folding.md) - pass " +
            "--use_remote_server true to query the ColabFold API instead."
        )
    }

    // Accept a single file, a glob, or a directory of FASTA files (one FASTA
    // file = one prediction unit). Resolved
    // synchronously via the `file()` DSL function (which returns a List<Path>
    // for a glob pattern) rather than a lazy Channel.fromPath, so the
    // monomer-only guard below can fail the whole run up front, before any
    // process is scheduled. NOTE: Channel.fromPath(...).toList().getVal()
    // looks equivalent but deadlocks here - the dataflow scheduler isn't
    // pumped yet at this point in workflow-body execution, so avoid it for
    // this kind of eager, pre-channel validation.
    def p = params.input
    def resolved = file(p).isDirectory() ? file("${p}/*.{fasta,fa,faa}") : file(p)
    List input_paths = (resolved instanceof List) ? resolved : [resolved]

    if (!input_paths) {
        error("fold.nf: no FASTA files found for --input '${p}'")
    }

    // *** Phase 1: monomer only. *** Fail fast (rather than silently
    // mis-folding) on any multi-record FASTA - see plans/
    // fold-nf-multi-method-folding.md §4 for the deferred multimer strategy.
    input_paths.each { f ->
        def n_chains = 0
        f.eachLine { line -> if (line.startsWith('>')) { n_chains++ } }
        if (n_chains != 1) {
            error(
                "fold.nf: ${f} has ${n_chains} FASTA records - multimer folding not yet " +
                "supported in fold.nf (coming in Phase 2; see plans/fold-nf-multi-method-folding.md). " +
                'Split multi-chain inputs into single-record FASTA files for now.'
            )
        }
    }

    // meta.id is the FASTA stem, matching AF2's own per-target output
    // directory naming - the MSA -> predict wiring in
    // subworkflows/local/alphafold2.nf relies on that matching.
    ch_input = Channel.fromList(input_paths).map { f -> [[id: f.baseName, n_chains: 1], f] }

    FOLD_MSA(ch_input, methods, params.msa_method)

    if ('af2' in methods) {
        ALPHAFOLD2(FOLD_MSA.out.af2_msas)
    }
    if ('boltz' in methods) {
        BOLTZ_FOLD(FOLD_MSA.out.for_boltz)
    }
    if ('rf3' in methods) {
        ROSETTAFOLD3_FOLD(FOLD_MSA.out.for_rf3)
    }
    if ('protenix' in methods) {
        PROTENIX_FOLD(FOLD_MSA.out.for_protenix)
    }

    ///////////////////////////////////////////////////////////////////////////
    workflow.onComplete = {
        def fold_params_json = [:]

        fold_params_json['params'] = ParamsHelper.paramsToMap(params)

        fold_params_json['workflow'] = [
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
        def json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(fold_params_json))

        new File(output_file).parentFile.mkdirs()
        new File(output_file).text = json_string

        log.info("Parameters saved to: ${output_file}")
    }
}
