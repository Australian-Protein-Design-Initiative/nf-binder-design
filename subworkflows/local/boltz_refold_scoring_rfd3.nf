#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
BOLTZ_REFOLD_SCORING_RFD3 subworkflow

Refolds RFD3/MPNN designs with Boltz-2 and calculates 
RMSD metrics and BindCraft-style scores.
*/

include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_BOLTZ_COMPLEX } from '../../modules/local/rfd/bindcraft_scoring'
include { PDB_TO_FASTA } from '../../modules/local/common/pdb_to_fasta'
include { BOLTZ_COMPARE_COMPLEX } from '../../modules/local/common/boltz_compare_complex'
include { BOLTZ_COMPARE_BINDER_MONOMER } from '../../modules/local/common/boltz_compare_binder_monomer'
include { MMSEQS_COLABFOLDSEARCH } from '../../modules/local/common/mmseqs_colabfoldsearch'

workflow BOLTZ_REFOLD_SCORING_RFD3 {

    take:
    ch_cifs_with_meta          // tuple(meta, pdb/cif)
    binder_chain               // String: binder chain ID (e.g., 'B' in RFD3)
    target_chain               // String: target chain ID (e.g., 'A' in RFD3)
    refold_max                 // Integer or false: max designs to refold
    refold_create_target_msa   // Boolean: create target MSA
    refold_use_msa_server      // Boolean: use MSA server
    refold_target_fasta        // Path or false: target FASTA file
    refold_target_templates    // Path or false: target templates directory
    colabfold_envdb            // Path or false: ColabFold env database
    uniref30                   // Path or false: UniRef30 database
    outdir                     // String: output directory

    main:

    // Limit to refold_max if specified
    ch_filtered_for_refold = refold_max
        ? ch_cifs_with_meta.take(refold_max)
        : ch_cifs_with_meta

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
            name: 'rmsd_complex_vs_rf3.tsv',
            storeDir: "${outdir}/boltz_refold/rmsd",
            keepHeader: true,
            skip: 1,
        )

    ch_monomer_vs_af2ig_rmsd = BOLTZ_COMPARE_BINDER_MONOMER.out.rmsd_monomer_vs_af2ig
        .map { meta, tsv_file -> tsv_file }
        .collectFile(
            name: 'rmsd_monomer_vs_rf3.tsv',
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

    ch_ipsae_tsv = BOLTZ_COMPARE_COMPLEX.out.ipsae_tsv
    ch_ipsae_byres_tsv = BOLTZ_COMPARE_COMPLEX.out.ipsae_byres_tsv

    extra_scores = Channel.empty()

    emit:
    rmsd_target_aligned_binder = ch_target_aligned_rmsd
    rmsd_complex = ch_complex_rmsd
    rmsd_monomer_vs_rf3 = ch_monomer_vs_af2ig_rmsd
    rmsd_monomer_vs_complex = ch_monomer_vs_complex_rmsd
    boltz_scores_complex = ch_complex_confidence
    boltz_scores_monomer = ch_monomer_confidence
    boltz_extra_scores = extra_scores
    ipsae_tsv = ch_ipsae_tsv
    ipsae_byres_tsv = ch_ipsae_byres_tsv

    // Also emit all the individual outputs in case needed
    boltz_complex_pdb = BOLTZ_COMPARE_COMPLEX.out.pdb
    boltz_monomer_pdb = BOLTZ_COMPARE_BINDER_MONOMER.out.pdb
}
