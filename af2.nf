#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Standalone AlphaFold2 prediction workflow: folds one or more FASTA files with
the Monash custom AlphaFold2 CUDA-12 container. Each FASTA file is one
prediction unit - multiple records in a single file are multimer chains, not
separate targets (the opposite convention to boltz_pulldown.nf's per-record
splitFasta).

Usage:
  nextflow run af2.nf --input 'input/*.fasta' --outdir results -profile slurm,m3

The two prediction stages (CPU jackhmmer/hhblits MSA generation, then GPU
structure prediction) are separate, reusable processes under
modules/local/af2/ so a future subworkflow can `include` them as an
alternative to Boltz-2.
*/

// Default parameters
params.help = false
params.input = false
params.outdir = 'results'

params.alphafold2_db_path = '/mnt/datasets/alphafold/alphafold_20240229'
params.alphafold2_model_preset = 'monomer_ptm' // monomer|monomer_ptm|monomer_casp14|multimer
params.alphafold2_db_preset = 'full_dbs'
params.alphafold2_max_template_date = '2024-01-01'
params.alphafold2_random_seed = false // set an int to fix the data pipeline's random seed
params.alphafold2_num_predictions_per_model = 1 // multimer only
params.alphafold2_models_to_relax = 'best' // all|best|none

// require_gpu / gpu_devices / gpu_allocation_detect_process_regex are already
// defined with pipeline-wide defaults in nextflow.config

include { ALPHAFOLD2_JACKHMMER_MSA } from './modules/local/af2/alphafold2_jackhmmer_msa'
// The GPU predict stage is now a shared subworkflow (also used by fold.nf) -
// see subworkflows/local/alphafold2.nf. MSA generation stays here unchanged:
// af2.nf always uses the native jackhmmer/hhblits stage, whereas fold.nf's
// FOLD_MSA can also feed the predict stage a ColabFold-a3m bridge.
include { ALPHAFOLD2 } from './subworkflows/local/alphafold2'

workflow {
    if (params.help || params.input == false) {
        log.info(
            """
        ==================================================================
        🧬 ALPHAFOLD2 WORKFLOW 🧬
        ==================================================================

        Predict structures for one or more FASTA files with the Monash custom
        AlphaFold2 CUDA-12 container.

        Required arguments:
            --input                            Single FASTA file, glob, or directory of FASTA files.
                                                Each file is one prediction unit (multiple records in
                                                one file are treated as multimer chains).

        Optional arguments:
            --outdir                           Output directory [default: ${params.outdir}]
            --alphafold2_db_path                AlphaFold2 database directory [default: ${params.alphafold2_db_path}]
            --alphafold2_model_preset            monomer|monomer_ptm|monomer_casp14|multimer [default: ${params.alphafold2_model_preset}]
            --alphafold2_db_preset               full_dbs|reduced_dbs [default: ${params.alphafold2_db_preset}]
            --alphafold2_max_template_date        Maximum template release date [default: ${params.alphafold2_max_template_date}]
            --alphafold2_random_seed              Fix the data pipeline's random seed [default: unset, randomly generated]
            --alphafold2_num_predictions_per_model Predictions per model, multimer only [default: ${params.alphafold2_num_predictions_per_model}]
            --alphafold2_models_to_relax          all|best|none [default: ${params.alphafold2_models_to_relax}]

            --gpu_devices                        GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
            --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]

        Example:
            nextflow run af2.nf --input 'input/*.fasta' --outdir results -profile slurm,m3

        Note: monomer presets emit one structure per model (5 for the standard
        preset); reducing that count would require model-subset selection,
        which is not exposed as a native flag in this container build and is
        therefore not currently supported. --num_multimer_predictions_per_model
        only applies when --alphafold2_model_preset multimer.

        """.stripIndent()
        )
        exit(1)
    }

    // Accept a single file, a glob, or a directory of FASTA files. meta.id is
    // the FASTA stem, which AlphaFold also uses to name its per-target output
    // directory - the MSA -> predict wiring below relies on that matching.
    def p = params.input
    ch_input = (
        file(p).isDirectory()
            ? Channel.fromPath("${p}/*.{fasta,fa,faa}")
            : Channel.fromPath(p)
    )
        .map { f -> [[id: f.baseName], f] }

    ALPHAFOLD2_JACKHMMER_MSA(ch_input)
    ALPHAFOLD2(ALPHAFOLD2_JACKHMMER_MSA.out.msa)

    ///////////////////////////////////////////////////////////////////////////
    workflow.onComplete = {
        def af2_params_json = [:]

        af2_params_json['params'] = ParamsHelper.paramsToMap(params)

        af2_params_json['workflow'] = [
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
        def json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(af2_params_json))

        new File(output_file).parentFile.mkdirs()
        new File(output_file).text = json_string

        log.info("Parameters saved to: ${output_file}")
    }
}
