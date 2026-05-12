#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
FOLDSEEK_ANNOTATE_CATH process

Annotates FoldSeek results with CATH hierarchy descriptions.
Extracts CATH codes from FoldSeek target names (e.g. af_A0A096MJB5_3_128_1.25.40.10)
and adds columns for Class, Architecture, Topology, and Homologous Superfamily.

Requires the CATH names file (cath-names.txt) which can be downloaded via
FOLDSEEK_DOWNLOAD_CATH_NAMES or provided via params.foldseek_cath_names_path.
*/

process FOLDSEEK_ANNOTATE_CATH {
    tag "${meta.id}"
    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5'
    publishDir "${params.outdir}/foldseek/${db_name}", mode: 'copy',
        pattern: '{foldseek_results_annotated.tsv,foldseek_results_annotated.tsv.gz}'

    input:
    tuple val(meta), path(input_tsv)
    path cath_names
    val db_name

    output:
    path "foldseek_results_annotated.tsv", emit: tsv_annotated, optional: true
    path "foldseek_results_annotated.tsv.gz", emit: tsv_annotated_gz, optional: true

    script:
    def gzip_output = params.foldseek_gzip_output ?: false
    def script_dir = workflow.projectDir ?: '.'
    """
    python3 "${script_dir}/bin/foldseek/annotate_cath.py" \
        --input "${input_tsv}" \
        --cath-names "${cath_names}" \
        ${gzip_output ? '--gzip --output foldseek_results_annotated.tsv.gz' : '--output foldseek_results_annotated.tsv'}
    """
}
