#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
BOLTZ_REFOLD_SCORING_RFD3 subworkflow

Refolds RFD3/MPNN designs with Boltz-2 and calculates
RMSD metrics and BindCraft-style scores.

Thin wrapper around BOLTZ_REFOLD_CORE; RFD3 designs are passed straight through
(no AF2IG filtering), with RMSD outputs labelled `rf3`.
*/

include { BOLTZ_REFOLD_CORE } from './boltz_refold_core'

workflow BOLTZ_REFOLD_SCORING_RFD3 {

    take:
    ch_cifs_with_meta          // tuple(meta, pdb/cif)
    binder_chain               // String: binder chain ID (e.g., 'B' in RFD3)
    target_chain               // String: target chain ID (e.g., 'A' in RFD3)
    refold_max                 // Integer or false: max designs to refold
    refold_create_target_msa   // Boolean: create target MSA
    refold_use_msa_server      // Boolean: use MSA server
    refold_alignment           // Path or false: external A3M file; bypasses MMseqs2 when set
    refold_target_fasta        // Path or false: target FASTA file
    refold_target_templates    // Path or false: target templates directory
    colabfold_envdb            // Path or false: ColabFold env database
    uniref30                   // Path or false: UniRef30 database
    outdir                     // String: output directory

    main:

    BOLTZ_REFOLD_CORE(
        ch_cifs_with_meta,
        binder_chain,
        target_chain,
        refold_max,
        refold_create_target_msa,
        refold_use_msa_server,
        refold_alignment,
        refold_target_fasta,
        refold_target_templates,
        colabfold_envdb,
        uniref30,
        'rf3',
        'boltz_refold/extra_scores/',
        outdir,
    )

    emit:
    rmsd_target_aligned_binder = BOLTZ_REFOLD_CORE.out.rmsd_target_aligned_binder
    rmsd_complex = BOLTZ_REFOLD_CORE.out.rmsd_complex
    rmsd_monomer_vs_rf3 = BOLTZ_REFOLD_CORE.out.rmsd_monomer_vs_design
    rmsd_monomer_vs_complex = BOLTZ_REFOLD_CORE.out.rmsd_monomer_vs_complex
    boltz_scores_complex = BOLTZ_REFOLD_CORE.out.boltz_scores_complex
    boltz_scores_monomer = BOLTZ_REFOLD_CORE.out.boltz_scores_monomer
    boltz_extra_scores = BOLTZ_REFOLD_CORE.out.boltz_extra_scores
    ipsae_tsv = BOLTZ_REFOLD_CORE.out.ipsae_tsv
    ipsae_byres_tsv = BOLTZ_REFOLD_CORE.out.ipsae_byres_tsv

    // Also emit all the individual outputs in case needed
    boltz_complex_pdb = BOLTZ_REFOLD_CORE.out.boltz_complex_pdb
    boltz_monomer_pdb = BOLTZ_REFOLD_CORE.out.boltz_monomer_pdb
}
