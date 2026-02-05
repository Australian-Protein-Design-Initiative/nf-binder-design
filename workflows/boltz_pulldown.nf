#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
Boltz Pulldown workflow for pairwise target-binder predictions

Usage via main.nf:
  nextflow run main.nf --method boltz_pulldown --targets targets.fasta --binders binders.fasta

*/

// Default parameters
params.targets = false
params.binders = false
params.outdir = "results"
params.use_msa_server = false
params.create_binder_msa = false
params.create_target_msa = false
params.templates = false

params.gpu_devices = ''
params.gpu_allocation_detect_process_regex = "(python.*/app/dl_binder_design/af2_initial_guess/predict\\.py|python.*/app/BindCraft/bindcraft\\.py|boltz predict|python.*/app/RFdiffusion/scripts/run_inference\\.py)"

include { BOLTZ } from '../modules/local/common/boltz'
include { MMSEQS_COLABFOLDSEARCH } from '../modules/local/common/mmseqs_colabfoldsearch'
include { BOLTZ_PULLDOWN_REPORTING } from '../modules/local/common/boltz_pulldown_reporting'

// Helper function to sanitize strings for filenames
def sanitize(name) {
    return name.replaceAll(/[^a-zA-Z0-9_.-]/, "_")
}

process CREATE_BOLTZ_YAML {
    tag "${target_meta.id}_and_${binder_meta.id}"

    container "ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5"

    input:
    tuple val(target_meta), path(target_msa), val(binder_meta), path(binder_msa)
    path templates

    output:
    tuple val(meta), path(yaml), path(target_msa), path(binder_msa)

    script:
    def id = "${sanitize(target_meta.id)}_and_${sanitize(binder_meta.id)}"
    def use_msa_server_flag = params.use_msa_server ? "--use_msa_server" : ""
    def templates_flag = params.templates ? "--templates '${templates}'" : ""

    /* we use 'flag files' empty_target_msa and boltz_will_make_target_msa
       to indicate wheather we want no MSA (empty), or allow Boltz to generate one
       via a remote mmseqs2 server (boltz_will_make_target_msa).
       If another (.a3m) file is provided, we use that as the MSA - this will generally
       come from a mmseqs2 search task in the pipeline (eg MMSEQS_COLABFOLDSEARCH).

       An alternative approach here might be to have mutliple variants of
       CREATE_BOLTZ_YAML (CREATE_BOLTZ_YAML_NO_MSA, CREATE_BOLTZ_YAML_MSA, 
       CREATE_BOLTZ_YAML_BINDER_MSA_ONLY, CREATE_BOLTZ_YAML_TARGET_MSA_ONLY) and a
       conditional in the main pipeline that selects the appropriate variant.
    */
    def target_msa_flag = ""
    if (target_msa.name == "empty_target_msa") {
        target_msa_flag = "--target_msa empty"
    }
    else if (target_msa.name == "boltz_will_make_target_msa") {
        target_msa_flag = ""
    }
    else {
        target_msa_flag = "--target_msa '${target_msa}'"
    }

    def binder_msa_flag = ""
    if (binder_msa.name == "empty_binder_msa") {
        binder_msa_flag = "--binder_msa empty"
    }
    else if (binder_msa.name == "boltz_will_make_binder_msa") {
        binder_msa_flag = ""
    }
    else {
        binder_msa_flag = "--binder_msa '${binder_msa}'"
    }

    meta = [id: id, target: sanitize(target_meta.id), binder: sanitize(binder_meta.id)]
    yaml = "${id}.yml"
    """
    ${projectDir}/bin/create_boltz_yaml.py \
        --target_id '${target_meta.id}' \
        --target_seq '${target_meta.seq}' \
        --binder_id '${binder_meta.id}' \
        --binder_seq '${binder_meta.seq}' \
        ${target_msa_flag} \
        ${binder_msa_flag} \
        --output_yaml '${yaml}' \
        ${use_msa_server_flag} \
        ${templates_flag}
    """
}


process PARSE_BOLTZ_CONFIDENCE_JSON {
    tag "${meta.id}"
    container "ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.5"

    input:
    tuple val(meta), path(json_file), path(ipsae_tsv)

    output:
    stdout

    script:
    """
    python3 ${projectDir}/bin/parse_boltz_confidence.py \\
        --json "${json_file}" \\
        --id "${meta.id}" \\
        --target "${meta.target}" \\
        --binder "${meta.binder}" \\
        --merge-ipsae "${ipsae_tsv}"
    """
}

workflow BOLTZ_PULLDOWN {

    main:

    if (params.targets == false || params.binders == false) {
        log.info(
            """
        ==================================================================
        BOLTZ PULLDOWN PIPELINE
        ==================================================================
        
        Required arguments:
            --targets             FASTA file of target sequences
            --binders             FASTA file of binder sequences

        Optional arguments:
            --outdir              Output directory [default: ${params.outdir}]
            --create_binder_msa   Create MSA for binder (a3m format) [default: ${params.create_binder_msa}]
            --create_target_msa   Create MSA for target (a3m format) [default: ${params.create_target_msa}]
            --templates           Templates directory with .cif files [default: ${params.templates}]
            --use_msa_server      Use BOLTZ MSA server [default: ${params.use_msa_server}]
            --uniref30            UniRef30 database path [default: ${params.uniref30}]
            --colabfold_envdb     ColabFold environment database path [default: ${params.colabfold_envdb}]

            --gpu_devices         GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
            --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]
        """.stripIndent()
        )
        exit(1)
    }

    ch_targets_meta = Channel.fromPath(params.targets)
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }

    ch_binders_meta = Channel.fromPath(params.binders)
        .splitFasta(record: [id: true, seqString: true])
        .map { record -> [id: record.id, seq: record.seqString] }

    ch_targets_fasta_paths = Channel.fromPath(params.targets)
        .splitFasta(file: true)

    ch_binders_fasta_paths = Channel.fromPath(params.binders)
        .splitFasta(file: true)

    // [meta, fasta_path]
    ch_targets_fasta = ch_targets_meta.merge(ch_targets_fasta_paths)
    ch_binders_fasta = ch_binders_meta.merge(ch_binders_fasta_paths)

    if (params.create_target_msa && !params.use_msa_server) {
        ch_target_msas = MMSEQS_COLABFOLDSEARCH(ch_targets_fasta, params.colabfold_envdb, params.uniref30)
    }
    else if (params.create_target_msa && params.use_msa_server) {
        ch_target_msas = ch_targets_fasta.map { [it[0], file("${projectDir}/assets/dummy_files/boltz_will_make_target_msa")] }
    }
    else {
        // [meta, null]
        ch_target_msas = ch_targets_fasta.map { [it[0], file("${projectDir}/assets/dummy_files/empty_target_msa")] }
    }

    if (params.create_binder_msa && !params.use_msa_server) {
        ch_binder_msas = MMSEQS_COLABFOLDSEARCH(ch_binders_fasta, params.colabfold_envdb, params.uniref30)
    }
    else if (params.create_binder_msa && params.use_msa_server) {
        ch_binder_msas = ch_binders_fasta.map { [it[0], file("${projectDir}/assets/dummy_files/boltz_will_make_binder_msa")] }
    }
    else {
        // [meta, null]
        ch_binder_msas = ch_binders_fasta.map { [it[0], file("${projectDir}/assets/dummy_files/empty_binder_msa")] }
    }

    // All combinations of: [meta_target, a3m_path_target, meta_binder, a3m_path_binder]
    ch_pairs = ch_target_msas.combine(ch_binder_msas)

    ch_templates = params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates")

    CREATE_BOLTZ_YAML(ch_pairs, ch_templates)

    BOLTZ(CREATE_BOLTZ_YAML.out, ch_templates, "boltz_pulldown")

    PARSE_BOLTZ_CONFIDENCE_JSON(BOLTZ.out.confidence_json.join(BOLTZ.out.ipsae_tsv))

    ch_tsv_output = PARSE_BOLTZ_CONFIDENCE_JSON.out.collectFile(
        name: "boltz_pulldown.tsv",
        storeDir: "${params.outdir}/boltz_pulldown",
        keepHeader: true,
        skip: 1,
    )

    BOLTZ_PULLDOWN_REPORTING(
        file("${projectDir}/assets/boltz_pulldown_reporting.qmd"),
        ch_tsv_output,
    )

    emit:
    tsv = ch_tsv_output
}
