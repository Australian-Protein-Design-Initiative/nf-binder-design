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
} else if (params.method == "germinal") {
    include { GERMINAL } from './workflows/germinal'
} else if (params.method == "boltzgen") {
    include { BOLTZGEN } from './workflows/boltzgen'
} else if (params.method == "boltz_pulldown") {
    include { BOLTZ_PULLDOWN } from './workflows/boltz_pulldown'
} else if (params.method == "fold") {
    include { FOLD } from './workflows/fold'
} else if (params.method == "fold_pulldown") {
    include { FOLD_PULLDOWN } from './workflows/fold_pulldown'
} else if (params.method == "rfd3") {
    include { RFD3 } from './workflows/rfd3'
} else if (params.method == "foldseek") {
    include { FOLDSEEK } from './workflows/foldseek'
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
            rfd3            RFDiffusion3-based binder design
            bindcraft       BindCraft binder design
            germinal        Germinal antibody/nanobody design
            boltzgen        BoltzGen binder design
            boltz_pulldown  Boltz pulldown predictions
            fold            Multi-method structure folding (AF2/Boltz/RF3/Protenix)
            fold_pulldown   Multi-method target x binder pulldown
            foldseek        FoldSeek structural similarity search

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
    } else if (params.method == "germinal") {
        GERMINAL()
    } else if (params.method == "boltzgen") {
        BOLTZGEN()
    } else if (params.method == "boltz_pulldown") {
        BOLTZ_PULLDOWN()
    } else if (params.method == "fold") {
        FOLD()
    } else if (params.method == "fold_pulldown") {
        FOLD_PULLDOWN()
    } else if (params.method == "rfd3") {
        RFD3()
    } else if (params.method == "foldseek") {
        FOLDSEEK()
    } else {
        log.error("Unknown method: ${params.method}")
        log.info("Available methods: rfd, rfd_partial, rfd3, bindcraft, germinal, boltzgen, boltz_pulldown, fold, fold_pulldown, foldseek")
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

        writeGpuTrace()
    }
}

// Aggregate the per-task GPU records written by nfbd_record_gpu_trace
// (bin/gpu_lock.sh) into one file alongside Nextflow's own trace.
//
// This exists because Nextflow's trace cannot answer "which GPU ran this?".
// Its observers run in the head process, while the device is chosen inside the
// container, so the only place that knows is the task itself.
//
// Rows are collected from the work directory rather than accumulated in memory,
// which means a -resume run still reports the GPU its cached tasks ran on
// originally, rather than silently dropping them.
def writeGpuTrace() {
    // Wrapped whole: an exception escaping workflow.onComplete is reported as
    // "Failed to invoke `workflow.onComplete` event handler" and turns a
    // successful run into one that looks failed. A diagnostic file is never
    // worth that, so every failure degrades to a warning.
    try {
        writeGpuTraceUnsafe()
    }
    catch (Exception e) {
        log.warn("Could not write the GPU trace: ${e}")
    }
}

def writeGpuTraceUnsafe() {
    def trace_dir = new File(params.gpu_trace_dir ?: "${workflow.workDir}/.gpu_trace")
    if (!trace_dir.isDirectory()) {
        return
    }

    def columns = [
        'timestamp',
        'hash',
        'process',
        'hostname',
        'n_gpus',
        'gpu_index',
        'gpu_uuid',
        'gpu_name',
        'driver_version',
        'memory_total_mib',
        'cuda_visible_devices',
    ]

    // A row that cannot be read is skipped rather than fatal: the file may have
    // been removed by a `nextflow clean`, or live on a filesystem that went
    // away, and one unreadable task must not cost the other several hundred.
    //
    // readLines().find is used rather than text.trim(): trim() strips the
    // trailing tab that delimits an empty final field, which would silently
    // shorten the row.
    //
    // Rows are also checked for width. Anything else in the directory -- a
    // leftover from an older column layout, a file a user dropped there -- would
    // otherwise be copied into the output and silently break every reader of it.
    def rows = (trace_dir.listFiles() ?: [] as File[])
        .findAll { it.isFile() && it.name.endsWith('.tsv') }
        .collect { f ->
            try {
                f.readLines().find { line -> line }
            }
            catch (Exception e) {
                log.warn("Skipping unreadable GPU trace record ${f.name}: ${e.message}")
                null
            }
        }
        .findAll { it && it.split('\t', -1).size() == columns.size() }
        .sort()

    if (!rows) {
        return
    }

    def header = columns.join('\t')

    def out = new File(params.gpu_trace_file ?: "${params.outdir}/logs/gpu_trace.txt")
    out.parentFile?.mkdirs()
    out.text = ([header] + rows).join('\n') + '\n'

    log.info("GPU trace saved to: ${out} (${rows.size()} tasks)")
}
