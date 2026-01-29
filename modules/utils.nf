nextflow.enable.dsl = 2

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

process COMPRESS_GPU_STATS {
    executor 'local'
    publishDir "${params.outdir}/logs", mode: 'copy'
    
    input:
    path csv_file
    
    output:
    path "gpu_stats.csv.gz", emit: gz_file
    
    script:
    """
    gzip -9 -c ${csv_file} >gpu_stats.csv.gz
    """
}
