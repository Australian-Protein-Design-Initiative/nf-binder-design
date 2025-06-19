#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters
params.targets        = false
params.binders        = false
params.outdir         = "results"
params.use_msa_server = false
params.create_binder_msa = false
params.create_target_msa = false
params.templates = false

params.gpu_device = 'all'

include { BOLTZ } from './modules/boltz'
include { MMSEQS_COLABFOLDSEARCH } from './modules/mmseqs_colabfoldsearch'
include { BOLTZ_PULLDOWN_REPORTING } from './modules/boltz_pulldown_reporting.nf'

// Helper function to sanitize strings for filenames
def sanitize(name) {
    return name.replaceAll(/[^a-zA-Z0-9_.-]/, "_")
}

process CREATE_BOLTZ_YAML {
    tag "${target_meta.id}_and_${binder_meta.id}"

    container "ghcr.io/australian-protein-design-initiative/containers/mdanalysis:2.8.0"

    input:
    tuple val(target_meta), path(target_msa), val(binder_meta), path(binder_msa)
    path(templates)

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
    } else if (target_msa.name == "boltz_will_make_target_msa") {
        target_msa_flag = ""
    } else {
        target_msa_flag = "--target_msa '${target_msa}'"
    }

    def binder_msa_flag = ""
    if (binder_msa.name == "empty_binder_msa") {
        binder_msa_flag = "--binder_msa empty"
    } else if (binder_msa.name == "boltz_will_make_binder_msa") {
        binder_msa_flag = ""
    } else {
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
    container "ghcr.io/australian-protein-design-initiative/containers/mdanalysis:2.8.0"

    input:
    tuple val(meta), path(json_file)

    output:
    stdout

    script:
    """
    #!/usr/bin/env python
    
    #
    # This takes a boltz confidence*.json file and flattens it into a TSV
    # Additionally, we add the unique id, target and binder name columns
    # Nested chains_ptm and pair_chains_iptm become chains_0, pair_chains_iptm_0_0 etc

    import json
    import sys
    import pandas as pd

    with open("${json_file}", 'r') as f:
        data = json.load(f)

    # Flatten the nested dictionaries from the JSON data
    df_flat = pd.json_normalize(data, sep='_')

    # Add the metadata columns to the beginning of the DataFrame
    df_flat.insert(0, 'id', '${meta.id}')
    df_flat.insert(1, 'target', '${meta.target}')
    df_flat.insert(2, 'binder', '${meta.binder}')

    # Write the flattened data to stdout as a TSV
    df_flat.to_csv(sys.stdout, sep='\t', index=False, lineterminator='\\n')
    """
}

workflow {

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

            --gpu_device          GPU device to use [default: ${params.gpu_device}]
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
    } else if (params.create_target_msa && params.use_msa_server) {
        ch_target_msas = ch_targets_fasta.map { [it[0], file("${projectDir}/assets/dummy_files/boltz_will_make_target_msa")] }
    } else
    {
        // [meta, null]
        ch_target_msas = ch_targets_fasta.map { [it[0], file("${projectDir}/assets/dummy_files/empty_target_msa")] }
    }

    if (params.create_binder_msa && !params.use_msa_server) {
        ch_binder_msas = MMSEQS_COLABFOLDSEARCH(ch_binders_fasta, params.colabfold_envdb, params.uniref30)
    } else if (params.create_binder_msa && params.use_msa_server) {
        ch_binder_msas = ch_binders_fasta.map { [it[0], file("${projectDir}/assets/dummy_files/boltz_will_make_binder_msa")] }
    } else
    {
        // [meta, null]
        ch_binder_msas = ch_binders_fasta.map { [it[0], file("${projectDir}/assets/dummy_files/empty_binder_msa")] }
    }

    // All combinations of: [meta_target, a3m_path_target, meta_binder, a3m_path_binder]
    ch_pairs = ch_target_msas.combine(ch_binder_msas)

    // ch_pairs.view()

    CREATE_BOLTZ_YAML(ch_pairs, file(params.templates))

    BOLTZ(CREATE_BOLTZ_YAML.out, file(params.templates), params.gpu_device)

    // BOLTZ.out.confidence_json.view { meta, json ->
    //     "Finshed: ${meta.target} + ${meta.binder}"
    // }

    PARSE_BOLTZ_CONFIDENCE_JSON(BOLTZ.out.confidence_json)

    ch_tsv_output = PARSE_BOLTZ_CONFIDENCE_JSON.out
        .collectFile(name: "boltz_pulldown.tsv", 
                     storeDir: "${params.outdir}/boltz_pulldown", 
                     keepHeader: true, 
                     skip: 1)

    BOLTZ_PULLDOWN_REPORTING(file("${projectDir}/assets/boltz_pulldown_reporting.qmd"), 
                             ch_tsv_output)

    // TODO: Re-sort the table on iptm (index 5). 
    //       (An alternative would be to sort on confidence (index 3))
    // csvtk -t sort -k iptm:r ${params.outdir}/boltz_pulldown/boltz_pulldown.tsv >${params.outdir}/boltz_pulldown/boltz_pulldown.sorted.tsv

    /* TODO / IDEAS:
         - Allow PDB format templates (currently on mmCIF is supported) - use the ghcr.io/australian-protein-design-initiative/containers/cif-tools:latest
           container to convert PDBs to CIF files if needed
           - cif-tools pdb2cif: In order to use pdb2cif, PDBs need a HEADER line added at the top. Remove any EXPDTA lines.

	- Create symlinks that categorize output structures:
          - Target directories with the best structure for each binder (to that target).
          - best_by_binder directory, with the best structure for each binder.
          - Include links to these structures in report (possibly in existing tables)
          - https://github.com/jmbuhr/quarto-molstar viewer in report ?

	- Reporting (boltz_pulldown_report.qmd) 
          - Put a copy of the Quarto doc in the results folder
          - What is the gap between the best target for a binder, and the next best target (or distribution of ipTM scores per binder)
            Is it possible in many cases the 'best' target is so close to the second/third/fourth that we can't really say ?
            Replicates with different seeds might help here.
       
       - Option to run replicates with different random seed (or different MSA subsamples), and incorporate that into the reporting (eg, we get a distribution of scores for each target/binder pair, so we may be able to determine if targetA+binderA is significantly (p <= 0.05) better than targetB+binderA within prediction variation))

       - Support providing a directory of input PDBs (for benchmarking or binder validation)
         - This could be a single directory --input-pdbs and options --target-chain --binder-chain
         - Or, a --binder-pdbs and --target-pdbs (with a convenience script to split out chains for preparation)
         - We can extract sequences from the PDBs, and also do RMSD, TM score, DockQ scores from these

       - Calculate pDockQ score for each complex.
          - GPL version from AlphaPulldown: https://github.com/KosinskiLab/AlphaPulldown/blob/main/alphapulldown/analysis_pipeline/calculate_mpdockq.py
          - Apache version they pinched: https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py

          - For benchmarking against known structures, we can use real DockQ: 
              - https://github.com/bjornwallner/DockQ
          - There's a bioconda package with container: 
            http://datacache.galaxyproject.org/singularity/d/o/dockq%3A2.1.3--py312h031d066_0
            or 
            quay.io/biocontainers/dockq:2.1.3--py312h031d066_0

       - Calcuculate the iptm_ptm score as-per AlphaPulldown 
         - eg https://github.com/KosinskiLab/AlphaPulldown/blob/308f3ad38760d4e579aa5085316c3a443612fd0c/alphapulldown/folding_backend/alphafold_backend.py#L451
         - Is one of the existing Boltz scores almost equivalent to this ?
          
       - Additional interface analysis as-per Alphapulldown, or similar (eg Arpeggio, PISA ?) ?
         - https://github.com/KosinskiLab/AlphaPulldown/blob/main/alphapulldown/analysis_pipeline/pdb_analyser.py - needs pyrosetta
    */
}
