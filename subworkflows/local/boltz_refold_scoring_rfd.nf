#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
BOLTZ_REFOLD_SCORING subworkflow

Refolds filtered AF2 initial guess designs with Boltz-2 and calculates 
RMSD metrics and BindCraft-style scores.

This subworkflow is shared between the RFD and RFD_PARTIAL workflows.

The refold + RMSD + BindCraft scoring is delegated to BOLTZ_REFOLD_CORE (shared
with the RFD3 workflow); this wrapper adds the AF2IG score filtering up front and
the COMBINE_SCORES aggregation at the end.
*/

include { AF2IG_SCORE_FILTER } from '../../modules/local/rfd/af2ig_score_filter'
include { COMBINE_SCORES } from '../../modules/local/common/combine_scores'
include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_AF2IG } from '../../modules/local/rfd/bindcraft_scoring'
include { BOLTZ_REFOLD_CORE } from './boltz_refold_core'

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
    ch_ipsae_tsv = Channel.empty()
    ch_ipsae_byres_tsv = Channel.empty()

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

        // Refold + RMSD + BindCraft scoring (shared with RFD3)
        BOLTZ_REFOLD_CORE(
            ch_filtered_with_meta,
            binder_chain,
            target_chain,
            refold_max,
            refold_create_target_msa,
            refold_use_msa_server,
            false,
            refold_target_fasta,
            refold_target_templates,
            colabfold_envdb,
            uniref30,
            'af2ig',
            'rfd/af2_initial_guess/extra_scores/',
            outdir,
        )

        ch_ipsae_tsv = BOLTZ_REFOLD_CORE.out.ipsae_tsv
        ch_ipsae_byres_tsv = BOLTZ_REFOLD_CORE.out.ipsae_byres_tsv

        // Combine all the score files into a single TSV file
        COMBINE_SCORES(
            ch_af2ig_scores.collect(),
            BOLTZ_REFOLD_CORE.out.boltz_extra_scores.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            BOLTZ_REFOLD_CORE.out.boltz_scores_complex.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            BOLTZ_REFOLD_CORE.out.rmsd_monomer_vs_complex.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
            BOLTZ_REFOLD_CORE.out.rmsd_target_aligned_binder.ifEmpty(file("${projectDir}/assets/dummy_files/empty")),
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
            'rfd/af2_initial_guess/extra_scores/',
        )

        extra_scores = BINDCRAFT_SCORING_AF2IG.out.scores.collectFile(
            name: "${outdir}/rfd/af2_initial_guess/af2ig_extra_scores.tsv",
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
    ipsae_tsv = ch_ipsae_tsv
    ipsae_byres_tsv = ch_ipsae_byres_tsv
}
