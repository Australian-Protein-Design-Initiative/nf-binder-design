#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_DOWNLOAD_DB process

Downloads a FoldSeek database (e.g. CATH50, PDB, AFDB50) using `foldseek databases`.
The database is stored in the Nextflow work directory and passed downstream so
only one copy is downloaded per pipeline run (cached via Nextflow work cache).

When params.foldseek_databases_path is set, the downloaded database is published
there for reuse across runs. On subsequent runs, if the database already exists
at that path, the download is skipped entirely.

Available databases: see `foldseek databases --help`.
*/

process FOLDSEEK_DOWNLOAD_DB {
    tag "${db_name}"
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'
    cache 'lenient'

    // Publish downloaded DB to foldseek_databases_path when configured.
    // saveAs returns null when path is not set, so nothing is published.
    publishDir params.foldseek_databases_path ?: '', mode: 'copy', overwrite: false,
        saveAs: { params.foldseek_databases_path ? "${db_name}" : null }

    input:
    val db_name   // e.g. 'CATH50'

    output:
    path "db", emit: db
    val db_name, emit: db_name_out

    script:
    def args = task.ext.args ?: ''
    """
    # Some db names contain a "/" (e.g. Alphafold/UniProt50), so the output
    # prefix db/${db_name} needs its parent dir created first, not just "db".
    # `dirname` on a name without a slash returns ".", so this is a no-op
    # (mkdir -p db/.) for names like CATH50 or PDB.
    mkdir -p "db/\$(dirname "${db_name}")"
    foldseek databases ${db_name} db/${db_name} ./tmp ${args}
    """
}
