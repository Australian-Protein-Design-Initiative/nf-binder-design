#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_DOWNLOAD_CATH_NAMES process

Downloads the CATH names file (cath-names.txt) from the CATH FTP server.
This file maps CATH hierarchy codes (e.g. 1.10.8.10) to human-readable
descriptions at Class, Architecture, Topology, and Homologous Superfamily levels.

When params.foldseek_databases_path is set, the file is published there for
reuse across runs. On subsequent runs, if the file already exists at that path,
the download is skipped entirely.
*/

process FOLDSEEK_DOWNLOAD_CATH_NAMES {
    tag "cath-names"
    container 'quay.io/biocontainers/foldseek:10.941cd33--h5021889_1'
    cache 'lenient'

    publishDir params.foldseek_databases_path ?: '', mode: 'copy', overwrite: false,
        saveAs: { params.foldseek_databases_path ? "cath-names.txt" : null }

    input:
    val dummy   // trigger only when needed

    output:
    path "cath-names.txt", emit: names_file

    script:
    """
    wget -q -O cath-names.txt \
        "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-names.txt"
    """
}
