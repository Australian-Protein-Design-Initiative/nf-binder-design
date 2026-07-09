#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Germinal binder design workflow

Usage via main.nf:
  nextflow run main.nf --method germinal --germinal_config configs/pdl1_vhh_protenix.yaml

*/

params.outdir = 'results'
params.germinal_config = false
params.germinal_pdb_dir = false
params.germinal_experiment_name = 'germinal_design'
params.germinal_n_traj = 10
params.germinal_batch_size = 1
params.germinal_max_passing_designs = 10000
params.germinal_max_hallucinated_trajectories = 10000
params.germinal_gpu_allocation_detect_process_regex = '(run_germinal\\.py|protenix|chai)'
params.require_gpu = true
params.gpu_devices = ''

include { GERMINAL as GERMINAL_PROCESS } from '../modules/local/germinal/germinal'
include { GERMINAL_MERGE } from '../modules/local/germinal/germinal_merge'

def resolveGerminalPdbDir(config_path) {
    if (params.germinal_pdb_dir) {
        return file(params.germinal_pdb_dir, checkIfExists: true)
    }

    def inferred = file("${config_path.parent}/../pdbs").normalize()
    if (!inferred.exists()) {
        error "Could not infer --germinal_pdb_dir from config location. Expected ${inferred} or pass --germinal_pdb_dir explicitly."
    }
    return inferred
}

workflow GERMINAL {

    main:

    if (!params.germinal_config) {
        log.info"""
        ==================================================================
        GERMINAL NEXTFLOW WRAPPER
        ==================================================================

        Required arguments:
            --germinal_config              Germinal Hydra config YAML (combined or partial)

        Optional arguments:
            --outdir                       Output directory [default: ${params.outdir}]
            --germinal_pdb_dir             Directory with target PDB and optional nb.pdb scaffold
                                           [default: ../pdbs relative to config directory]
            --germinal_experiment_name     Experiment name for output subdirectory [default: ${params.germinal_experiment_name}]
            --germinal_n_traj              Total number of trajectory attempts [default: ${params.germinal_n_traj}]
            --germinal_batch_size          Trajectories per batch (parallel task) [default: ${params.germinal_batch_size}]
            --germinal_max_passing_designs Max accepted designs per batch [default: ${params.germinal_max_passing_designs}]
            --germinal_max_hallucinated_trajectories
                                           Max hallucinated trajectories per batch [default: ${params.germinal_max_hallucinated_trajectories}]

            --require_gpu                  Fail if no GPU detected [default: ${params.require_gpu}]
            --gpu_devices                  GPU devices for allocation [default: ${params.gpu_devices}]
            --germinal_gpu_allocation_detect_process_regex
                                           Regex for busy GPU process detection [default: ${params.germinal_gpu_allocation_detect_process_regex}]

        """.stripIndent()
        exit 1
    }

    if (params.germinal_n_traj < 1) {
        error "germinal_n_traj must be >= 1"
    }
    if (params.germinal_batch_size < 1) {
        error "germinal_batch_size must be >= 1"
    }

    ch_config = Channel.fromPath(params.germinal_config, checkIfExists: true).first()
    ch_pdb_dir = ch_config.map { config ->
        resolveGerminalPdbDir(config)
    }

    def design_indices = 0..(params.germinal_n_traj - 1)
    def batches = design_indices.collate(params.germinal_batch_size)

    ch_batch_info = Channel.from(batches.withIndex())
        .map { batch, index -> [index, batch.size()] }

    GERMINAL_PROCESS(
        ch_batch_info,
        ch_config,
        ch_pdb_dir,
        params.germinal_experiment_name,
        params.germinal_max_passing_designs,
        params.germinal_max_hallucinated_trajectories
    )

    GERMINAL_MERGE(
        GERMINAL_PROCESS.out.batch_dir.collect(),
    )

    emit:
    all_trajectories = GERMINAL_MERGE.out.all_trajectories
    accepted_designs = GERMINAL_MERGE.out.accepted_designs
    failure_counts = GERMINAL_MERGE.out.failure_counts
    accepted_pdbs = GERMINAL_PROCESS.out.accepted_pdbs
}
