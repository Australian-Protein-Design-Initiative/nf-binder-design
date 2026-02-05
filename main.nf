#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-binder-design
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Protein binder design pipeline with multiple methods
----------------------------------------------------------------------------------------
*/

// Method parameter for workflow selection
params.method = false
params.outdir = 'results'

// Conditional includes based on --method parameter
if (params.method == "rfd") {
    include { RFD } from './workflows/rfd'
} else if (params.method == "rfd_partial") {
    include { RFD_PARTIAL } from './workflows/rfd_partial'
} else if (params.method == "bindcraft") {
    include { BINDCRAFT } from './workflows/bindcraft'
} else if (params.method == "boltzgen") {
    include { BOLTZGEN } from './workflows/boltzgen'
} else if (params.method == "boltz_pulldown") {
    include { BOLTZ_PULLDOWN } from './workflows/boltz_pulldown'
}

def paramsToMap(params) {
    def map = [:]
    params.each { key, value ->
        if (value instanceof Path || value instanceof File) {
            map[key] = value.toString()
        }
        else if (!(value instanceof Closure) && !(key in [
            'class',
            'launchDir',
            'projectDir',
            'workDir',
        ])) {
            map[key] = value
        }
    }
    return map
}

workflow {

    main:

    // Show help if no method specified
    if (params.method == false) {
        log.info(
            """
        ==================================================================
        PROTEIN BINDER DESIGN PIPELINE
        ==================================================================

        Usage: nextflow run main.nf --method <method> [options]

        Available methods:
            rfd             RFDiffusion-based binder design
            rfd_partial     RFDiffusion partial diffusion for binder optimization
            bindcraft       BindCraft binder design
            boltzgen        BoltzGen binder design
            boltz_pulldown  Boltz pulldown predictions

        Example:
            nextflow run main.nf --method rfd --input_pdb target.pdb --rfd_n_designs 10

        For method-specific help, run with --method <method> and no other arguments.

        """.stripIndent()
        )
        exit(1)
    }

    // Dispatch to appropriate workflow
    if (params.method == "rfd") {
        RFD()
    } else if (params.method == "rfd_partial") {
        RFD_PARTIAL()
    } else if (params.method == "bindcraft") {
        BINDCRAFT()
    } else if (params.method == "boltzgen") {
        BOLTZGEN()
    } else if (params.method == "boltz_pulldown") {
        BOLTZ_PULLDOWN()
    } else {
        log.error("Unknown method: ${params.method}")
        log.info("Available methods: rfd, rfd_partial, bindcraft, boltzgen, boltz_pulldown")
        exit(1)
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
