#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
BOLTZ_REFOLD_SCORING subworkflow

Refolds filtered AF2 initial guess designs with Boltz-2 and calculates 
RMSD metrics and BindCraft-style scores.

This subworkflow is shared between the RFD and RFD_PARTIAL workflows.
*/

include { AF2IG_SCORE_FILTER } from '../../modules/local/rfd/af2ig_score_filter'
include { COMBINE_SCORES } from '../../modules/local/common/combine_scores'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_AF2IG } from '../../modules/local/bindcraft/bindcraft_scoring'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_BOLTZ_COMPLEX } from '../../modules/local/bindcraft/bindcraft_scoring'
include { PDB_TO_FASTA } from '../../modules/local/common/pdb_to_fasta'
include { BOLTZ_COMPARE_COMPLEX } from '../../modules/local/common/boltz_compare_complex'
include { BOLTZ_COMPARE_BINDER_MONOMER } from '../../modules/local/common/boltz_compare_binder_monomer'
include { MMSEQS_COLABFOLDSEARCH } from '../../modules/local/common/mmseqs_colabfoldsearch'

workflow BOLTZ_REFOLD_SCORING {

    take:
    ch_af2ig_pdbs_with_scores  // tuple(pdbs, scores)
    ch_af2ig_pdbs              // channel of PDB files
    ch_af2ig_scores            // channel of score files
    binder_chain               // String: binder chain ID (e.g., 'A')
    target_chain               // String: target chain ID (e.g., 'B')
    refold_af2ig_filters       // String or false: filter expression
    refold_max                 // Integer or false: max designs to refold
    refold_create_target_msa   // Boolean: create target MSA
    refold_use_msa_server      // Boolean: use MSA server
    refold_target_fasta        // Path or false: target FASTA file
    refold_target_templates    // Path or false: target templates directory
    colabfold_envdb            // Path or false: ColabFold env database
    uniref30                   // Path or false: UniRef30 database
    outdir                     // String: output directory

    main:

    // Define output channel variable
    ch_combined_scores = Channel.empty()

    if (refold_af2ig_filters) {
        // Filter designs by score thresholds
        AF2IG_SCORE_FILTER(
            ch_af2ig_pdbs_with_scores.map { pdbs, scores -> scores },
            ch_af2ig_pdbs_with_scores.map { pdbs, scores -> pdbs },
            refold_af2ig_filters,
        )

        // Flatten and create metadata for each accepted PDB
        ch_filtered_with_meta = AF2IG_SCORE_FILTER.out.accepted
            .flatten()
            .map { pdb -> [[id: pdb.baseName], pdb] }

        // Limit to refold_max if specified
        ch_filtered_for_refold = refold_max
            ? ch_filtered_with_meta.take(refold_max)
            : ch_filtered_with_meta

        // Optionally create target MSAs
        if (refold_create_target_msa && !refold_use_msa_server) {
            if (refold_target_fasta) {
                // Use the refold_target_fasta file directly for MSA creation
                ch_target_fastas = ch_filtered_for_refold.map { meta, pdb ->
                    [meta, file(refold_target_fasta)]
                }
                ch_target_msas = MMSEQS_COLABFOLDSEARCH(
                    ch_target_fastas,
                    colabfold_envdb,
                    uniref30,
                )
            }
            else {
                // Create target FASTA files from PDBs for MSA creation
                ch_target_fastas = ch_filtered_for_refold.map { meta, pdb -> pdb }
                    | PDB_TO_FASTA(
                        target_chain
                    ).map { fasta ->
                        def basename = fasta.baseName.replaceAll(/_${target_chain}$/, '')
                        [[id: basename], fasta]
                    }
                ch_target_msas = MMSEQS_COLABFOLDSEARCH(
                    ch_target_fastas,
                    colabfold_envdb,
                    uniref30,
                )
            }
        }
        else if (refold_create_target_msa && refold_use_msa_server) {
            ch_target_msas = ch_filtered_for_refold.map { meta, pdb ->
                [meta, file("${projectDir}/assets/dummy_files/boltz_will_make_target_msa")]
            }
        }
        else {
            ch_target_msas = ch_filtered_for_refold.map { meta, pdb ->
                [meta, file("${projectDir}/assets/dummy_files/empty_target_msa")]
            }
        }

        // Run Boltz complex refolding with RMSD analysis
        BOLTZ_COMPARE_COMPLEX(
            ch_filtered_for_refold,
            binder_chain,
            target_chain,
            refold_create_target_msa,
            refold_use_msa_server,
            ch_target_msas.map { meta, target_msa -> target_msa },
            file("${projectDir}/assets/dummy_files/empty_binder_msa"),
            file(refold_target_templates ?: "${projectDir}/assets/dummy_files/empty_templates"),
            refold_target_fasta ? file(refold_target_fasta) : "${projectDir}/assets/dummy_files/empty",
        )

        // Run Boltz binder monomer prediction with RMSD analysis
        BOLTZ_COMPARE_BINDER_MONOMER(
            ch_filtered_for_refold.join(BOLTZ_COMPARE_COMPLEX.out.pdb).map { meta, af2ig_pdb, boltz_pdb ->
                [meta, af2ig_pdb, boltz_pdb]
            },
            binder_chain,
        )

        // Aggregate RMSD outputs
        ch_target_aligned_rmsd = BOLTZ_COMPARE_COMPLEX.out.rmsd_target_aligned
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_target_aligned_binder.tsv',
                storeDir: "${outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1,
            )

        ch_complex_rmsd = BOLTZ_COMPARE_COMPLEX.out.rmsd_complex
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_complex_vs_af2ig.tsv',
                storeDir: "${outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1,
            )

        ch_monomer_vs_af2ig_rmsd = BOLTZ_COMPARE_BINDER_MONOMER.out.rmsd_monomer_vs_af2ig
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_monomer_vs_af2ig.tsv',
                storeDir: "${outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1,
            )

        ch_monomer_vs_complex_rmsd = BOLTZ_COMPARE_BINDER_MONOMER.out.rmsd_monomer_vs_complex
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'rmsd_monomer_vs_complex.tsv',
                storeDir: "${outdir}/boltz_refold/rmsd",
                keepHeader: true,
                skip: 1,
            )

        // Aggregate confidence outputs
        ch_complex_confidence = BOLTZ_COMPARE_COMPLEX.out.confidence_tsv
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'boltz_scores_complex.tsv',
                storeDir: "${outdir}/boltz_refold",
                keepHeader: true,
                skip: 1,
            )

        ch_monomer_confidence = BOLTZ_COMPARE_BINDER_MONOMER.out.confidence_tsv
            .map { meta, tsv_file -> tsv_file }
            .collectFile(
                name: 'boltz_scores_binder_monomer.tsv',
                storeDir: "${outdir}/boltz_refold",
                keepHeader: true,
                skip: 1,
            )

        // BindCraft-score the Boltz-2 refolded complexes
        BINDCRAFT_SCORING_BOLTZ_COMPLEX(
            BOLTZ_COMPARE_COMPLEX.out.pdb.map { meta, pdb -> pdb },
            binder_chain,
            'default_4stage_multimer',
        )

        extra_scores = BINDCRAFT_SCORING_BOLTZ_COMPLEX.out.scores.collectFile(
            name: 'boltz_complex_extra_scores.tsv',
            storeDir: "${outdir}/boltz_refold",
            keepHeader: true,
            skip: 1,
        )

        // Combine all the score files into a single TSV file
        COMBINE_SCORES(
            ch_af2ig_scores.collect(),
            extra_scores.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            ch_complex_confidence.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            ch_monomer_vs_complex_rmsd.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            ch_target_aligned_rmsd.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            ch_af2ig_pdbs.collect(),
        )
        ch_combined_scores = COMBINE_SCORES.out.combined_scores
    }
    else {
        // No refolding - just score AF2IG designs directly
        BINDCRAFT_SCORING_AF2IG(
            ch_af2ig_pdbs,
            binder_chain,
            'default_4stage_multimer',
        )

        extra_scores = BINDCRAFT_SCORING_AF2IG.out.scores.collectFile(
            name: "${outdir}/af2_initial_guess/af2ig_extra_scores.tsv",
            keepHeader: true,
            skip: 1,
        )

        // Combine all the score files into a single TSV file
        COMBINE_SCORES(
            ch_af2ig_scores.collect(),
            extra_scores,
            "${projectDir}/assets/dummy_files/empty",
            "${projectDir}/assets/dummy_files/empty",
            "${projectDir}/assets/dummy_files/empty",
            ch_af2ig_pdbs.collect(),
        )
        ch_combined_scores = COMBINE_SCORES.out.combined_scores
    }

    emit:
    scores = ch_combined_scores
}
