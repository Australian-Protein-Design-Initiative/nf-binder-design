#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Standalone EnGens conformational clustering of an existing structure ensemble.

Use this when you already have a folder (or glob) of .cif / .pdb files and only
want the EnGens report + representative conformations - not fold.nf prediction.
For post-prediction clustering inside fold.nf, leave --skip_engens unset.

Usage:
  nextflow run engens.nf --input path/to/structures/ --outdir results -profile slurm,m3
  nextflow run engens.nf --input 'predictions/*.cif' --id UL119_domain --outdir results
*/

params.help = false
params.input = false
params.outdir = 'results'
params.id = false // target id for publishDir / report; default = input dir basename or 'engens'

params.engens_dimred = 'umap'
params.engens_clustering = 'hdbscan'
params.engens_min_structures = 3
params.engens_max_clusters = 10
params.engens_gmm_ic = 'aic'
params.engens_seed = false

include { ENGENS } from './modules/local/engens/engens_cluster'

workflow {
    if (params.help || params.input == false) {
        log.info(
            """
        ==================================================================
        ENGENS WORKFLOW (conformational clustering of structure ensembles)
        ==================================================================

        Cluster an existing folder/glob of .cif / .pdb structures with EnGens
        (UMAP + HDBSCAN by default) and write results/engens/<id>/clusters.html.

        Required arguments:
            --input                            Directory of structures, or a glob of .cif/.pdb files.

        Optional arguments:
            --outdir                           Output directory [default: ${params.outdir}]
            --id                               Target id for results/engens/<id>/ [default: input
                                                directory basename, or 'engens' for a bare glob]
            --engens_clustering                hdbscan (default), gmm, km, or comma-separated
                                                combinations [default: ${params.engens_clustering}]
            --engens_dimred                    Dimensionality reduction (umap)
                                                [default: ${params.engens_dimred}]
            --engens_min_structures            Minimum usable structures before clustering
                                                [default: ${params.engens_min_structures}]
            --engens_max_clusters              Upper bound for auto cluster-count search
                                                [default: ${params.engens_max_clusters}]
            --engens_gmm_ic                    GMM information criterion: aic|bic
                                                [default: ${params.engens_gmm_ic}]
            --engens_seed                      Optional RNG seed for UMAP / clustering [default: unset]

        Example:
            nextflow run engens.nf --input results/fold/predictions/ --id UL119_domain \\
                --outdir results -profile slurm,m3

        """.stripIndent()
        )
        exit(params.help ? 0 : 1)
    }

    def engens_clustering = params.engens_clustering.toString().split(',').collect { it.trim().toLowerCase() }
    engens_clustering.each { c ->
        if (!(c in ['gmm', 'km', 'hdbscan'])) {
            error("engens.nf: unknown --engens_clustering entry '${c}' (valid: gmm, km, hdbscan)")
        }
    }
    if (!(params.engens_dimred.toString().toLowerCase() in ['umap'])) {
        error("engens.nf: --engens_dimred must be 'umap' (got '${params.engens_dimred}')")
    }
    if (!(params.engens_gmm_ic.toString().toLowerCase() in ['aic', 'bic'])) {
        error("engens.nf: --engens_gmm_ic must be 'aic' or 'bic' (got '${params.engens_gmm_ic}')")
    }
    if ((params.engens_min_structures as int) < 2) {
        error("engens.nf: --engens_min_structures must be >= 2 (got '${params.engens_min_structures}')")
    }
    if ((params.engens_max_clusters as int) < 2) {
        error("engens.nf: --engens_max_clusters must be >= 2 (got '${params.engens_max_clusters}')")
    }

    def input_path = file(params.input)
    def resolved
    def default_id
    if (input_path instanceof List) {
        // CLI glob already expanded by Nextflow's file()
        resolved = input_path
        default_id = 'engens'
    } else if (input_path.isDirectory()) {
        resolved = file("${params.input}/*.{cif,pdb,CIF,PDB}")
        default_id = input_path.name
    } else if (input_path.exists()) {
        resolved = [input_path]
        default_id = input_path.baseName
    } else {
        // Glob pattern that has not been expanded yet
        resolved = file(params.input)
        default_id = 'engens'
    }

    List structure_paths = (resolved instanceof List) ? resolved : [resolved]
    structure_paths = structure_paths.findAll { f ->
        def name = f.getName().toString().toLowerCase()
        (name.endsWith('.cif') || name.endsWith('.pdb')) && f.exists()
    }

    if (!structure_paths) {
        error("engens.nf: no .cif/.pdb files found for --input '${params.input}'")
    }

    // Duplicate basenames collide in the ENGENS work dir staging loop.
    def basenames = structure_paths.collect { it.getName().toString() }
    def dupes = basenames.countBy { it }.findAll { _name, n -> n > 1 }.keySet()
    if (dupes) {
        error(
            "engens.nf: duplicate structure basenames under --input (first: ${dupes.take(5).join(', ')}). " +
            "Flatten/rename files so each basename is unique."
        )
    }

    def target_id = params.id ? params.id.toString() : default_id.toString()
    log.info("engens.nf: clustering ${structure_paths.size()} structures as id='${target_id}'")

    ch_structures = Channel.value([[id: target_id], structure_paths])
    ch_qmd = Channel.value(file("${projectDir}/assets/engens/engens-analysis.qmd"))
    ENGENS(ch_qmd, ch_structures)
}
