#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters
params.help = false
params.targets_fasta  = false
params.binders_fasta  = false
params.target_structures = false
params.binder_structures = false
params.complex_structures = false
params.target_chains = false
params.binder_chains = false
params.outdir         = "results"
params.use_msa_server = false
params.create_binder_msa = false
params.create_target_msa = false
params.templates = false
params.colabfold_envdb = false
params.uniref30 = false

params.gpu_devices = ''
params.gpu_allocation_detect_process_regex = "(python.*/app/dl_binder_design/af2_initial_guess/predict\\.py|python.*/app/BindCraft/bindcraft\\.py|boltz predict|python.*/app/RFdiffusion/scripts/run_inference\\.py)"

include { PDB_TO_FASTA as TARGET_PDB_TO_FASTA } from './modules/pdb_to_fasta.nf'
include { PDB_TO_FASTA as BINDER_PDB_TO_FASTA } from './modules/pdb_to_fasta.nf'
include { CREATE_BOLTZ_YAML as CREATE_BOLTZ_YAML_COMPLEX } from './modules/create_boltz_yaml.nf'
include { CREATE_BOLTZ_YAML_MONOMER as CREATE_BOLTZ_YAML_TARGET_MONOMER } from './modules/create_boltz_yaml.nf'
include { CREATE_BOLTZ_YAML_MONOMER as CREATE_BOLTZ_YAML_BINDER_MONOMER } from './modules/create_boltz_yaml.nf'
include { BOLTZ as BOLTZ_COMPLEX } from './modules/boltz'
include { BOLTZ as BOLTZ_TARGET_MONOMER } from './modules/boltz'
include { BOLTZ as BOLTZ_BINDER_MONOMER } from './modules/boltz'
include { MMSEQS_COLABFOLDSEARCH } from './modules/mmseqs_colabfoldsearch'
include { BOLTZ_PULLDOWN_REPORTING } from './modules/boltz_pulldown_reporting.nf'
include { RMSD4ALL as RMSD4ALL_COMPLEXES } from './modules/rmsd4all.nf'
include { RMSD4ALL as RMSD4ALL_TARGET_MONOMERS } from './modules/rmsd4all.nf'
include { RMSD4ALL as RMSD4ALL_BINDER_MONOMERS } from './modules/rmsd4all.nf'


process PARSE_BOLTZ_CONFIDENCE_JSON {
    tag "${meta.id}"
    container "ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.4"

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

// Function to display help message
def showHelp() {
    log.info(
        """
        ==================================================================
        BOLTZ PULLDOWN PIPELINE
        ==================================================================
        
        Required arguments (choose one for each):
            --targets_fasta       FASTA file of target sequences
            --target_structures   Directory containing target PDB files (requires --target_chains)

            --binders_fasta       FASTA file of binder sequences  
            --binder_structures   Directory containing binder PDB files (requires --binder_chains)

            --complex_structures  Directory containing binder+target complex PDB files 
                                 (can be combined with exactly one of target_* or binder_*; the other side is
                                  extracted from these complexes. For the side coming from complexes, chain
                                  defaults are A for target and B for binder if not specified.)

        Chain specification (required with structure inputs):
            --target_chains       Chain IDs to extract from target PDBs (e.g., "A" or "A,B")
            --binder_chains       Chain IDs to extract from binder PDBs (e.g., "A" or "A,B")

        Optional arguments:
            --outdir              Output directory [default: ${params.outdir}]
            --create_binder_msa   Create MSA for binder (a3m format) [default: ${params.create_binder_msa}]
            --create_target_msa   Create MSA for target (a3m format) [default: ${params.create_target_msa}]
            --templates           Templates directory with .pdb and .cif files [default: ${params.templates}]
            --use_msa_server      Use BOLTZ MSA server [default: ${params.use_msa_server}]
            --uniref30            UniRef30 database path [default: ${params.uniref30}]
            --colabfold_envdb     ColabFold environment database path [default: ${params.colabfold_envdb}]

            --gpu_devices         GPU devices to use (comma-separated list or 'all') [default: ${params.gpu_devices}]
            --gpu_allocation_detect_process_regex  Regex pattern to detect busy GPU processes [default: ${params.gpu_allocation_detect_process_regex}]
            --help                Show this help message and exit
        """.stripIndent()
    )
}

workflow {

    // Check for help flag
    if (params.help) {
        showHelp()
        exit(0)
    }

    // Validation logic
    def has_targets_fasta = params.targets_fasta != false
    def has_target_structures = params.target_structures != false
    def has_binders_fasta = params.binders_fasta != false
    def has_binder_structures = params.binder_structures != false
    def has_complex_structures = params.complex_structures != false
    
    def has_any_target_input = has_targets_fasta || has_target_structures
    def has_any_binder_input = has_binders_fasta || has_binder_structures
    
    // Check that we have at least one target input method
    if (!has_targets_fasta && !has_target_structures && !has_complex_structures) {
        log.error("ERROR: Must provide either --targets_fasta OR --target_structures OR --complex_structures")
        showHelp()
        exit(1)
    }
    
    // Check that we have at least one binder input method
    if (!has_binders_fasta && !has_binder_structures && !has_complex_structures) {
        log.error("ERROR: Must provide either --binders_fasta OR --binder_structures OR --complex_structures")
        showHelp()
        exit(1)
    }
    
    // Check mutual exclusivity for targets
    if (has_targets_fasta && has_target_structures) {
        log.error("ERROR: Cannot specify both --targets_fasta and --target_structures. Choose one.")
        exit(1)
    }
    
    // Check mutual exclusivity for binders
    if (has_binders_fasta && has_binder_structures) {
        log.error("ERROR: Cannot specify both --binders_fasta and --binder_structures. Choose one.")
        exit(1)
    }
    
    // Mixed complex mode validation: if complex_structures is provided, exactly one of target_* or binder_* must also be provided
    if (has_complex_structures) {
        if (has_any_target_input && has_any_binder_input) {
            log.error("ERROR: When using --complex_structures, specify exactly one of target_* OR binder_* (not both).")
            exit(1)
        }
        if (!has_any_target_input && !has_any_binder_input) {
            log.error("ERROR: When using --complex_structures, you must also specify either target_* OR binder_* (exactly one).")
            exit(1)
        }
    }
    
    // Check that chain parameters are provided when using structures
    if (has_target_structures && params.target_chains == false) {
        log.error("ERROR: --target_chains is required when using --target_structures")
        exit(1)
    }
    
    if (has_binder_structures && params.binder_chains == false) {
        log.error("ERROR: --binder_chains is required when using --binder_structures")
        exit(1)
    }
    
    // Determine which side (target or binder) comes from complex, and set chains/defaults accordingly
    def use_complex_for_target = has_complex_structures && !has_any_target_input
    def use_complex_for_binder = has_complex_structures && !has_any_binder_input
    def complex_target_chains = null
    def complex_binder_chains = null
    if (use_complex_for_target) {
        complex_target_chains = (params.target_chains != false) ? params.target_chains : "A"
        if (params.target_chains == false) {
            log.info("Using default target chain: A")
        }
    }
    if (use_complex_for_binder) {
        complex_binder_chains = (params.binder_chains != false) ? params.binder_chains : "B"
        if (params.binder_chains == false) {
            log.info("Using default binder chain: B")
        }
    }
    
    // Check that chain parameters are not provided when using FASTA
    if (has_targets_fasta && params.target_chains != false) {
        log.error("ERROR: --target_chains can only be used with --target_structures or --complex_structures, not --targets_fasta")
        exit(1)
    }
    
    if (has_binders_fasta && params.binder_chains != false) {
        log.error("ERROR: --binder_chains can only be used with --binder_structures or --complex_structures, not --binders_fasta")
        exit(1)
    }

    showHelp()

    // TODO: Current issues with WIP implementation below:
    //       - Providing PDBs are templates fails - it seems despite the docs, boltz wants
    //         .cif files ?
    //       - When providing target_structures, we should be able to omit --target_chains 
    //         (and vice versa for binders), in which case all chains are used (fixing this will likely
    //         require changes to the pdb_to_fasta -> create_boltz_yaml scripts and flow)
    //       - Should we simplify and all just --targets and --binders with mixed pdb and fasta inputs ?
    //         In this case the monomer vs input monomer RMSDs would only apply to the input structures.
    //       - Given that the processing of pdb or fasta -> boltz yaml is relatively low resource, maybe
    //         we should combine this into a single script, and possibly call it in a bash for loop
    //         in a single process, that then outputs a channel something like:
    //            [meta, config, various_, files_, config_needs]

    // TODO: RMSD comparison implementation - WIP: 
    //      - Add prediction of each monomer structure (target and binder) by default, and provide reporting on 
    //       RMSDs between original models (in a complex) and the monomer forms.
    //       The --skip-target-monomers and/or --skip-binder-monomers options turns off the monomer 
    //       prediction.
    //       We generally DON'T want to use the provided input structures as templates for these
    //       monomer predictions. We MAY want to use the --templates (other structures) in some cases
    //       to help guide the monomer predictions (eg in cases where the experimental target structure
    //       isn't solved)
    //       
    //       In here, if we have --target-structures and/or --binder-structures,
    //       and optionally --target-chains and --binder-chains, extract the relevant sequences from 
    //       each PDB file to feed into Boltz Pulldown.
    //       We want to pair up predicted models with original structures at the end and
    //       calculate the CA RMSD, possibly TM-score. Provide a table, and heatmap or similar
    //       in the Quarto report.
    //       
    //       The --complex-structures option can be used to specify existing binder+target complex 
    //       This is primarily intended for denovo binder re-prediction.
    //       
    //       External --templates can still be provided (eg experimental target structures of homologs 
    //       to help guide predictions). If we did want to ever bias the predictions like 'initial guess',
    //       we could always provide denovo complexes as separate binder and target monomer structures 
    //       to as --templates. Probably we never want to do this.
    //       
    //       Should the binder complex + binder monomers be in a separate pipeline for simplicity ?
    //       
    //       Consider a samplesheet mode, that has columns:
    //       structure_id, fasta_path, structure_path, chains, type (target or binder)
    //       where one of fasta_path or structure_path+chains is required, but both is an error.
    //       

    // Create target channels based on input type
    if (has_targets_fasta) {
        // FASTA input
        ch_targets_meta = Channel.fromPath(params.targets_fasta)
            .splitFasta(record: [id: true, seqString: true])
            .map { record -> [id: record.id, seq: record.seqString] }

        ch_targets_fasta_paths = Channel.fromPath(params.targets_fasta)
            .splitFasta(file: true)

        // [meta, fasta_path]
        ch_targets_fasta = ch_targets_meta.merge(ch_targets_fasta_paths)
    } else if (has_target_structures) {
        // Structure input - extract sequences from PDBs
        ch_target_pdbs = Channel.fromPath("${params.target_structures}/*.pdb")
        //ch_target_pdbs.view()
        
        // [[id, seq, type], fasta_file]
        ch_targets_fasta = TARGET_PDB_TO_FASTA(ch_target_pdbs, params.target_chains)
            .map { fasta_file ->
                [[id: fasta_file.text.split('\n')[0][1..-1], 
                  seq: fasta_file.text.split('\n')[1..-1].join(''), 
                  type: 'target'], 
                  fasta_file]
            }
        //ch_targets_fasta.view()
    } else if (has_complex_structures) {
        // Complex structure input for target - extract target sequences from PDBs
        ch_complex_pdbs = Channel.fromPath("${params.complex_structures}/*.pdb")
        
        // [[id, seq, type], fasta_file]
        ch_targets_fasta = TARGET_PDB_TO_FASTA(ch_complex_pdbs, complex_target_chains)
            .map { fasta_file ->
                [[id: fasta_file.text.split('\n')[0][1..-1], 
                  seq: fasta_file.text.split('\n')[1..-1].join(''), 
                  type: 'target'], 
                  fasta_file]
            }
        //ch_targets_fasta.view()
    }

    // Create binder channels based on input type
    if (has_binders_fasta) {
        // FASTA input
        ch_binders_meta = Channel.fromPath(params.binders_fasta)
            .splitFasta(record: [id: true, seqString: true])
            .map { record -> [id: record.id, seq: record.seqString] }

        ch_binders_fasta_paths = Channel.fromPath(params.binders_fasta)
            .splitFasta(file: true)

        // [meta, fasta_path]
        ch_binders_fasta = ch_binders_meta.merge(ch_binders_fasta_paths)
    } else if (has_binder_structures) {
        // Structure input - extract sequences from PDBs
        ch_binder_pdbs = Channel.fromPath("${params.binder_structures}/*.pdb")
        
        // [[id, seq, type], fasta_file]
        ch_binders_fasta = BINDER_PDB_TO_FASTA(ch_binder_pdbs, params.binder_chains)
            .map { fasta_file ->
                [[id: fasta_file.text.split('\n')[0][1..-1], 
                  seq: fasta_file.text.split('\n')[1..-1].join(''), 
                  type: 'binder'], 
                  fasta_file]
            }
    } else if (use_complex_for_binder) {
        // Complex structure input for binder - extract binder sequences from PDBs
        ch_complex_pdbs = Channel.fromPath("${params.complex_structures}/*.pdb")
        
        // [[id, seq, type], fasta_file]
        ch_binders_fasta = BINDER_PDB_TO_FASTA(ch_complex_pdbs, complex_binder_chains)
            .map { fasta_file ->
                [[id: fasta_file.text.split('\n')[0][1..-1], 
                  seq: fasta_file.text.split('\n')[1..-1].join(''), 
                  type: 'binder'], 
                  fasta_file]
            }
    }

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

    CREATE_BOLTZ_YAML_COMPLEX(ch_pairs, params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates"))

    BOLTZ_COMPLEX(
        CREATE_BOLTZ_YAML_COMPLEX.out, 
        params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates"), 
        tuple('boltz_pulldown', 'complex')
    )

    // Monomer predictions: run Boltz on targets and binders individually
    ch_targets_for_monomer = ch_target_msas.map { [it[0], it[1], 'target'] }
    ch_binders_for_monomer = ch_binder_msas.map { [it[0], it[1], 'binder'] }

    CREATE_BOLTZ_YAML_TARGET_MONOMER(
        ch_targets_for_monomer, 
        params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates"), 
    )
    CREATE_BOLTZ_YAML_BINDER_MONOMER(
        ch_binders_for_monomer, 
        params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates"), 
    )

    // Create 4-element tuples for BOLTZ process: (meta, yaml, target_msa, binder_msa)
    ch_target_monomer_for_boltz = CREATE_BOLTZ_YAML_TARGET_MONOMER.out
        .map { meta, yaml, msa -> [meta, yaml, msa, file("${projectDir}/assets/dummy_files/empty_binder_msa")] }
    
    ch_binder_monomer_for_boltz = CREATE_BOLTZ_YAML_BINDER_MONOMER.out
        .map { meta, yaml, msa -> [meta, yaml, file("${projectDir}/assets/dummy_files/empty_target_msa"), msa] }

    BOLTZ_TARGET_MONOMER(
        ch_target_monomer_for_boltz, 
        params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates"), 
        tuple('boltz_pulldown', 'target_monomer')
    )
    BOLTZ_BINDER_MONOMER(
        ch_binder_monomer_for_boltz, 
        params.templates ? file(params.templates) : file("${projectDir}/assets/dummy_files/empty_templates"),
        tuple('boltz_pulldown', 'binder_monomer')
    )

    // Log progress
    BOLTZ_COMPLEX.out.confidence_json.view { meta, json ->
        "Finished complex: ${meta.target} + ${meta.binder}"
    }
    BOLTZ_TARGET_MONOMER.out.confidence_json.view { meta, json ->
        "Finished target monomer: ${meta.target}"
    }
    BOLTZ_BINDER_MONOMER.out.confidence_json.view { meta, json ->
        "Finished binder monomer: ${meta.binder}"
    }

    PARSE_BOLTZ_CONFIDENCE_JSON(BOLTZ_COMPLEX.out.confidence_json)

    ch_tsv_output = PARSE_BOLTZ_CONFIDENCE_JSON.out
        .collectFile(name: "boltz_pulldown.tsv", 
                     storeDir: "${params.outdir}/boltz_pulldown", 
                     keepHeader: true, 
                     skip: 1)

    // Prepare lists of predicted structures for RMSD
    ch_binder_monomer_structs = BOLTZ_BINDER_MONOMER.out.predicted_structure
        .filter { it[0].type == 'binder' }
        .map { it[1] }.collect()
        
    ch_target_monomer_structs = BOLTZ_TARGET_MONOMER.out.predicted_structure
        .filter { it[0].type == 'target' }
        .map { it[1] }.collect()

    ch_complex_pred_structs = BOLTZ_COMPLEX.out.predicted_structure
        .map { it[1] }.collect()

    // Determine reference directories for RMSD comparisons
    def original_binder_dir = null
    def original_target_dir = null
    
    if (has_binder_structures) {
        original_binder_dir = params.binder_structures
    } else if (use_complex_for_binder) {
        original_binder_dir = params.complex_structures
    }
    
    if (has_target_structures) {
        original_target_dir = params.target_structures
    } else if (use_complex_for_target) {
        original_target_dir = params.complex_structures
    }

    // Run RMSD: binder monomers vs original binders (if originals exist)
    if (original_binder_dir) {
        ch_binder_monomer_rmsd_tsv = RMSD4ALL_BINDER_MONOMERS(
            ch_binder_monomer_structs,
            original_binder_dir,
            params.binder_chains ? params.binder_chains : '',
            "predicted-binder-monomers_vs_input-binders.tsv"
        ).rmsd_tsv
    } else {
        ch_binder_monomer_rmsd_tsv = channel.empty()
    }
    
    // Run RMSD: target monomers vs original targets (if originals exist)
    if (original_target_dir) {
        ch_target_monomer_rmsd_tsv = RMSD4ALL_TARGET_MONOMERS(
            ch_target_monomer_structs,
            original_target_dir,
            params.target_chains ? params.target_chains : '',
            "predicted-target-monomers_vs_input-targets.tsv"
        ).rmsd_tsv
    } else {
        ch_target_monomer_rmsd_tsv = channel.empty()
    }

    // Run RMSD: predicted complexes vs provided complex structures (if provided)
    if (has_complex_structures) {
        ch_complex_rmsd_tsv = RMSD4ALL_COMPLEXES(
            ch_complex_pred_structs,
            params.complex_structures,
            '',
            'predicted-complexes_vs_input-complexes.tsv'
        ).rmsd_tsv
    } else {
        ch_complex_rmsd_tsv = channel.empty()
    }

    ch_rmsd_tsv = channel.of(file("${projectDir}/assets/dummy_files/empty")).mix(
        ch_binder_monomer_rmsd_tsv, 
        ch_target_monomer_rmsd_tsv, 
        ch_complex_rmsd_tsv)

    BOLTZ_PULLDOWN_REPORTING(file("${projectDir}/assets/boltz_pulldown_reporting.qmd"), 
                             ch_tsv_output, 
                             ch_rmsd_tsv)

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

def paramsToMap(params) {
    def map = [:]
    params.each { key, value ->
        if (value instanceof Path || value instanceof File) {
            map[key] = value.toString()
        } else if (!(value instanceof Closure) && !(key in [
            'class', 'launchDir', 'projectDir', 'workDir'])) {
            map[key] = value
        }
    }
    return map
}

workflow.onComplete {
    // Write the pipeline parameters to a JSON file
    def params_json = [:]

    params_json['params'] = paramsToMap(params)

    params_json['workflow'] = [
        name: workflow.manifest.name,
        version: workflow.manifest.version,
        runName: workflow.runName,
        start: workflow.start.format('yyyy-MM-dd HH:mm:ss'),
        complete: workflow.complete.format('yyyy-MM-dd HH:mm:ss'),
        duration: workflow.duration,
        success: workflow.success
    ]

    def output_file = "${params.outdir}/params.json"
    def json_string = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params_json))
    
    new File(params.outdir).mkdirs()
    new File(output_file).text = json_string
    
    log.info "Pipeline parameters saved to: ${output_file}"
}