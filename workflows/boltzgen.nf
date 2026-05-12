#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
BoltzGen binder design workflow

Usage via main.nf:
  nextflow run main.nf --method boltzgen --config_yaml config.yaml

*/

// Default parameters
params.config_yaml = false
params.outdir = 'results'
params.design_name = false
params.protocol = 'protein-anything'
params.num_designs = 100
params.batch_size = 10
params.budget = 10
params.devices = false
params.num_workers = false
params.inverse_fold_num_sequences = false
params.alpha = false
params.filter_biased = false
params.metrics_override = []
params.additional_filters = []
params.size_buckets = []
params.refolding_rmsd_threshold = false
params.do_foldseek = false
params.foldseek_cath_names_path = false

include { BOLTZGEN_DESIGN } from '../modules/local/boltzgen/boltzgen_design'
include { BOLTZGEN_INVERSE_FOLDING } from '../modules/local/boltzgen/boltzgen_inverse_folding'
include { BOLTZGEN_FOLDING } from '../modules/local/boltzgen/boltzgen_folding'
include { BOLTZGEN_DESIGN_FOLDING } from '../modules/local/boltzgen/boltzgen_design_folding'
include { BOLTZGEN_AFFINITY } from '../modules/local/boltzgen/boltzgen_affinity'
include { BOLTZGEN_MERGE } from '../modules/local/boltzgen/boltzgen_merge'
include { BOLTZGEN_ANALYSIS } from '../modules/local/boltzgen/boltzgen_analysis'
include { BOLTZGEN_FILTERING } from '../modules/local/boltzgen/boltzgen_filtering'

include { PARSE_BOLTZGEN_CONFIG; detectParams; buildFilteringArgs } from '../modules/local/boltzgen/boltzgen_utils'
include { FOLDSEEK_SEARCH } from '../subworkflows/local/foldseek_search'
include { FOLDSEEK_PREPARE_QUERIES } from '../modules/local/foldseek/foldseek_prepare_queries'

// Function to validate design_name does not end with a number
def validateDesignName(design_name) {
    if (design_name == null) {
        return true  // No validation needed if parameter is not provided
    }

    // Check for empty or whitespace-only string
    if (design_name.trim().isEmpty()) {
        error "Invalid design_name format: '${design_name}'. Design name cannot be empty or whitespace only."
    }

    // Check if design_name ends with a digit
    def pattern = ~/\d$/
    if (design_name =~ pattern) {
        error "Invalid design_name format: '${design_name}'. Design name cannot end with a number. Use a name like 'mydesign' instead of 'mydesign_1'."
    }

    return true
}

workflow BOLTZGEN {

    main:

    // Show help message
    if (params.config_yaml == false) {
        log.info(
            """
        ==================================================================
        BOLTZGEN PIPELINE
        ==================================================================

        Required arguments:
            --config_yaml              Path to BoltzGen YAML config file

        Optional arguments:
            --outdir                      Output directory [default: ${params.outdir}]
            --design_name                 Name of the design, used for output file prefixes [default: auto-derived from config_yaml basename]
            --protocol                    Protocol type (protein-anything, peptide-anything, protein-small_molecule, nanobody-anything) [default: ${params.protocol}]
            --num_designs                 Total number of designs to generate [default: ${params.num_designs}]
            --batch_size                  Number of designs per batch [default: ${params.batch_size}]
            --budget                      Final diversity-optimized set size [default: ${params.budget}]
            --devices                     Number of GPU devices [default: unspecified]
            --num_workers                 Number of DataLoader workers [default: unspecified]
            --inverse_fold_num_sequences  Number of sequences per design in inverse folding step [default: unspecified]
            --alpha                       Trade-off for sequence diversity selection: 0.0=quality-only, 1.0=diversity-only
            --filter_biased               Remove amino-acid composition outliers (default: true, use --filter_biased false to disable)
            --metrics_override            Per-metric inverse-importance weights for ranking. Format: metric_name=weight (e.g., 'plip_hbonds_refolded=4' 'delta_sasa_refolded=2')
            --additional_filters          Extra hard filters. Format: feature>threshold or feature<threshold (e.g., 'design_ALA>0.3' 'design_GLY<0.2')
            --size_buckets                Optional constraint for maximum number of designs in size ranges. Format: min-max:count (e.g., '10-20:5' '20-30:10')
            --refolding_rmsd_threshold     Threshold used for RMSD-based filters (lower is better)

        FoldSeek (optional, enabled with --do_foldseek):
            --do_foldseek                 Enable FoldSeek structural similarity search on final designs
            --foldseek_database           Database to search [default: CATH50]
                                           Local search (default):
                                             - CATH50 (default) — CATH domain DB, combined AF2+PDB at 50% seq.id.
                                             - PDB              — Protein Data Bank
                                             - Alphafold/UniProt — Full AF2 database (~700GB)
                                             - Alphafold/UniProt50 — AF2 clustered at 50% seq.id.
                                             - Alphafold/UniProt50-minimal — AF2 50% clusters (reps only)
                                             - Alphafold/Proteome — AF2 proteomes
                                             - Alphafold/Swiss-Prot — AF2 Swiss-Prot
                                             - ESMAtlas30       — ESM metagenomic atlas at 30% seq.id.
                                             - BFMD             — Big Fantastic Multimer Database
                                             - BFVD             — Big Fantastic Virus Database
                                           Remote search (--foldseek_use_webserver true):
                                             - pdb100           — PDB (complex-aware, includes multimers)
                                             - afdb50           — AlphaFold/UniProt50 (~37M structures)
                                             - afdb-swissprot   — AlphaFold/Swiss-Prot
                                             - afdb-proteome    — AlphaFold/Proteomes
                                             - cath50           — CATH50 domain database
                                             - bfmd             — BFMD multimers
                                             - BFVD             — BFVD viruses
                                             - mgnify_esm30     — MGnify-ESM30 metagenomic atlas
            --foldseek_databases_path     Path to local databases directory (auto-downloaded if unset)
            --foldseek_use_webserver      Use FoldSeek web API instead of local search [default: false]
            --foldseek_mode               Search mode [default: 3diaa]
            --foldseek_maxaccept          Max accepted alignments per query [default: 1]
            --foldseek_gzip_output        Gzip TSV output [default: false]
            --foldseek_include_html_output  Include HTML report [default: false]
            --foldseek_cath_names_path    Path to local cath-names.txt (auto-downloaded if unset)
                                          CATH annotation is automatic when database name starts with "CATH"

        """.stripIndent()
        )
        exit(1)
    }

    def design_name = params.design_name
    
    // Set design_name from config_yaml basename if not explicitly set
    if (!params.design_name) {
        def config_file = new File(params.config_yaml)
        def config_basename = config_file.getName()
        // Remove .yaml or .yml extension
        design_name = config_basename.replaceFirst(/\.(yaml|yml)$/, '')
    }

    // Validate design_name does not end with a number
    validateDesignName(design_name)

    // Validate config_yaml exists
    def config_file = file(params.config_yaml)
    def config_dir = config_file.parent

    ch_config_yaml = Channel.fromPath(params.config_yaml).first()
    ch_config_dir = Channel.value(config_dir.toString())

    // Parse config to collect all files that need staging (entity paths + scaffold inner files)
    // Runs as a process so PyYAML is available via a container
    PARSE_BOLTZGEN_CONFIG(ch_config_yaml, ch_config_dir)

    // stdout is a single string with one absolute path per line - split, deduplicate, resolve
    ch_input_files = PARSE_BOLTZGEN_CONFIG.out
        .map { stdout_text ->
            def paths = stdout_text.trim().split('\n').findAll { it.trim() }.unique()
            if (paths.isEmpty()) return []
            return paths.collect { p -> file(p.trim()) }
        }

    // Generate batch start indices - create separate channels to avoid double consumption
    def batch_indices = (0..params.num_designs - 1).findAll { it % params.batch_size == 0 }
    
    ch_batch_n_designs = Channel.from(batch_indices)
        .map { start_idx ->
            Math.min(params.batch_size, params.num_designs - start_idx)
        }
    
    ch_batch_start_idx = Channel.from(batch_indices)

    // Create channels for constant values
    ch_design_name = Channel.value(design_name)
    ch_protocol = Channel.value(params.protocol)
    ch_devices = Channel.value(params.devices)
    ch_num_workers = Channel.value(params.num_workers)
    ch_inverse_fold_num_sequences = Channel.value(params.inverse_fold_num_sequences)

    // Phase 1: Design (Parallel)
    BOLTZGEN_DESIGN(
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        ch_batch_n_designs,
        ch_batch_start_idx,
        ch_devices,
        ch_num_workers,
    )

    // Phase 2: Inverse Folding (Parallel)
    // Extract start_index from batch_dir path (batch_dir is batch_{start_index}/)
    ch_batch_with_start = BOLTZGEN_DESIGN.out.batch_dir.map { batch_dir ->
        def start_idx = batch_dir.name.replaceAll(/batch_/, '').toInteger()
        return [batch_dir, start_idx]
    }

    BOLTZGEN_INVERSE_FOLDING(
        ch_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        ch_batch_with_start.map { _batch_dir, start_idx -> start_idx },
        ch_devices,
        ch_num_workers,
        ch_inverse_fold_num_sequences,
    )

    // Phase 3: Folding (Parallel)
    ch_folding_batch_with_start = BOLTZGEN_INVERSE_FOLDING.out.batch_dir.map { batch_dir ->
        def start_idx = batch_dir.name.replaceAll(/batch_/, '').toInteger()
        return [batch_dir, start_idx]
    }

    BOLTZGEN_FOLDING(
        ch_folding_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        ch_folding_batch_with_start.map { _batch_dir, start_idx -> start_idx },
        ch_devices,
        ch_num_workers,
    )

    // Phase 4: Design Folding (Parallel, if applicable)
    if (params.protocol in ['protein-anything', 'protein-small_molecule']) {
        ch_design_folding_batch_with_start = BOLTZGEN_FOLDING.out.batch_dir.map { batch_dir ->
            def start_idx = batch_dir.name.replaceAll(/batch_/, '').toInteger()
            return [batch_dir, start_idx]
        }

        BOLTZGEN_DESIGN_FOLDING(
            ch_design_folding_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
            ch_config_yaml,
            ch_input_files,
            ch_design_name,
            ch_protocol,
            ch_design_folding_batch_with_start.map { _batch_dir, start_idx -> start_idx },
            ch_devices,
            ch_num_workers,
        )
        
        if (params.protocol == 'protein-small_molecule') {
            ch_affinity_batch_with_start = BOLTZGEN_DESIGN_FOLDING.out.batch_dir.map { batch_dir ->
                def start_idx = batch_dir.name.replaceAll(/batch_/, '').toInteger()
                return [batch_dir, start_idx]
            }

            BOLTZGEN_AFFINITY(
                ch_affinity_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
                ch_config_yaml,
                ch_input_files,
                ch_design_name,
                ch_protocol,
                ch_affinity_batch_with_start.map { _batch_dir, start_idx -> start_idx },
                ch_devices,
                ch_num_workers,
            )
            ch_design_folded = BOLTZGEN_AFFINITY.out.batch_dir
        }
        else {
            ch_design_folded = BOLTZGEN_DESIGN_FOLDING.out.batch_dir
        }
    }
    else {
        ch_design_folded = BOLTZGEN_FOLDING.out.batch_dir
    }

    // Phase 5: Analysis (Parallel)
    ch_analysis_batch_with_start = ch_design_folded.map { batch_dir ->
        def start_idx = batch_dir.name.replaceAll(/batch_/, '').toInteger()
        return [batch_dir, start_idx]
    }

    BOLTZGEN_ANALYSIS(
        ch_analysis_batch_with_start.map { batch_dir, _start_idx -> batch_dir },
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        ch_analysis_batch_with_start.map { _batch_dir, start_idx -> start_idx },
    )

    // Phase 6: Merge
    BOLTZGEN_MERGE(
        BOLTZGEN_ANALYSIS.out.batch_dir.collect()
    )

    // Phase 7: Filtering (Single process)
    ch_merged_dir = BOLTZGEN_MERGE.out.merged_dir.first()

    // Build filtering arguments string
    def filtering_args = buildFilteringArgs(
        params.alpha,
        params.filter_biased,
        params.metrics_override,
        params.additional_filters,
        params.size_buckets,
        params.refolding_rmsd_threshold
    )

    BOLTZGEN_FILTERING(
        ch_merged_dir,
        ch_config_yaml,
        ch_input_files,
        ch_design_name,
        ch_protocol,
        Channel.value(params.budget),
        Channel.value(filtering_args),
    )

    // Phase 8: FoldSeek structural search (optional)
    if (params.do_foldseek) {
        // Collect BoltzGen output structure files, then extract design chains
        ch_ranked_dir = BOLTZGEN_FILTERING.out.final_ranked_designs_dir.first()

        ch_boltzgen_files = ch_ranked_dir.flatMap { dir ->
            def files = []
            def dirFile = new File(dir.toString())
            dirFile.eachDir { batchDir ->
                if (batchDir.name.startsWith('final_') && batchDir.name.endsWith('_designs')) {
                    batchDir.listFiles()?.each { f ->
                        if (f.isFile() && (f.name.endsWith('.cif') || f.name.endsWith('.cif.gz') ||
                            f.name.endsWith('.pdb') || f.name.endsWith('.pdb.gz'))) {
                            files << file(f.absolutePath)
                        }
                    }
                }
            }
            files
        }

        FOLDSEEK_PREPARE_QUERIES(ch_boltzgen_files)

        ch_foldseek_pdbs = FOLDSEEK_PREPARE_QUERIES.out.design_chains
        ch_foldseek_meta = Channel.of([id: 'boltzgen_foldseek'])

        FOLDSEEK_SEARCH(ch_foldseek_pdbs, ch_foldseek_meta)

        FOLDSEEK_SEARCH.out.tsv.subscribe { f -> }
    }

    emit:
    filtered_dir = BOLTZGEN_FILTERING.out.final_ranked_designs_dir
    foldseek_tsv = params.do_foldseek ? FOLDSEEK_SEARCH.out.tsv : Channel.empty()
    foldseek_tsv_annotated = params.do_foldseek ? FOLDSEEK_SEARCH.out.tsv_annotated : Channel.empty()
    foldseek_html = params.do_foldseek ? FOLDSEEK_SEARCH.out.html : Channel.empty()
}
