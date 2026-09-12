#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
BOLTZ_REFOLD_CORE subworkflow

Shared core for refolding binder designs with Boltz-2 and calculating RMSD
metrics, Boltz confidence, ipSAE and BindCraft-style scores.

Both BOLTZ_REFOLD_SCORING (RFD/RFD_PARTIAL) and BOLTZ_REFOLD_SCORING_RFD3 wrap
this core. The only differences between those workflows are handled via
parameters: the RMSD output label (`af2ig` vs `rf3`), the BindCraft publish
subdirectory, and whether an external target MSA (`refold_alignment`) is used.
*/

// --boltz_refold_batch_size groups designs into a single `boltz predict` call; its default and
// the measurements behind it are in nextflow.config. Read here as null-safe, so this
// subworkflow still works if included somewhere that does not set it.

include { BINDCRAFT_SCORING as BINDCRAFT_SCORING_BOLTZ_COMPLEX } from '../../modules/local/rfd/bindcraft_scoring'
include { PDB_TO_FASTA } from '../../modules/local/common/pdb_to_fasta'
include { BOLTZ_COMPARE_COMPLEX } from '../../modules/local/common/boltz_compare_complex'
include { BOLTZ_COMPARE_BINDER_MONOMER } from '../../modules/local/common/boltz_compare_binder_monomer'
include { MMSEQS_COLABFOLDSEARCH } from '../../modules/local/common/mmseqs_colabfoldsearch'

workflow BOLTZ_REFOLD_CORE {

    take:
    ch_structures_with_meta    // tuple(meta, pdb/cif): designs to refold
    binder_chain               // String: binder chain ID in the input design
    target_chain               // String: target chain ID in the input design
    refold_max                 // Integer or false: max designs to refold
    refold_create_target_msa   // Boolean: create target MSA
    refold_use_msa_server      // Boolean: use MSA server
    refold_alignment           // Path or false: external A3M file; bypasses MMseqs2 when set
    refold_target_fasta        // Path or false: target FASTA file
    refold_target_templates    // Path or false: target templates directory
    colabfold_envdb            // Path or false: ColabFold env database
    uniref30                   // Path or false: UniRef30 database
    rmsd_label                 // String: label for RMSD-vs-design outputs (e.g. 'af2ig' or 'rf3')
    bindcraft_publish_subdir   // String: publishDir subdir for BindCraft complex scores
    outdir                     // String: output directory

    main:

    // Limit to refold_max if specified
    ch_filtered_for_refold = refold_max
        ? ch_structures_with_meta.take(refold_max)
        : ch_structures_with_meta

    // Target MSA: external file, or create via MMseqs2/server, or none
    if (refold_alignment) {
        ch_target_msas = ch_filtered_for_refold.map { meta, pdb ->
            [meta, file(refold_alignment)]
        }
    }
    else if (refold_create_target_msa && !refold_use_msa_server) {
        if (refold_target_fasta) {
            // Use the refold_target_fasta file directly for MSA creation
            ch_target_fastas = ch_filtered_for_refold.map { meta, pdb ->
                [meta, file(refold_target_fasta)]
            }
            ch_target_msas = MMSEQS_COLABFOLDSEARCH(
                ch_target_fastas,
                false,
                colabfold_envdb,
                uniref30,
                'boltz_refold/mmseqs2',
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
                false,
                colabfold_envdb,
                uniref30,
                'boltz_refold/mmseqs2',
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

    def boltz_batch_int = (params.boltz_refold_batch_size == null || params.boltz_refold_batch_size == false)
        ? 1
        : (params.boltz_refold_batch_size as int)
    if (boltz_batch_int < 1) {
        throw new Exception('--boltz_refold_batch_size must be >= 1')
    }

    // Group designs into batches, and group the target MSAs the same way so the two stay
    // paired. The MSA channel is collated rather than joined on meta: the process never
    // reads the staged MSA path (it is present only so a YAML can reference it by
    // basename), the MSA belongs to the target and so is the same for every design, and
    // collating preserves the positional pairing these two channels already relied on.
    // The extension travels with the batch because the process stages its inputs under
    // numbered names, which discards it (see BOLTZ_COMPARE_COMPLEX).
    ch_complex_batches = ch_filtered_for_refold
        .collate(boltz_batch_int, true)
        .map { batch ->
            tuple(
                batch.collect { it[0] },
                batch.collect { it[1] },
                batch.collect { it[1].extension },
            )
        }

    ch_target_msa_batches = ch_target_msas
        .map { meta, target_msa -> target_msa }
        .collate(boltz_batch_int, true)
        .map { batch -> batch.first() }

    // Run Boltz complex refolding with RMSD analysis
    use_target_msa = refold_create_target_msa || refold_alignment
    BOLTZ_COMPARE_COMPLEX(
        ch_complex_batches,
        binder_chain,
        target_chain,
        use_target_msa,
        refold_use_msa_server,
        ch_target_msa_batches,
        file("${projectDir}/assets/dummy_files/empty_binder_msa"),
        file(refold_target_templates ?: "${projectDir}/assets/dummy_files/empty_templates"),
        refold_target_fasta ? file(refold_target_fasta) : "${projectDir}/assets/dummy_files/empty",
    )

    // Split each batch back into per-design channels. Paths are rebuilt from the design id
    // rather than from output order, so a batch cannot mis-pair a design with another
    // design's results.
    def ch_complex_per_design = BOLTZ_COMPARE_COMPLEX.out.per_design_bundle.flatMap { metas, bundle ->
        def ml = metas instanceof List ? metas : [metas]
        ml.collect { m -> tuple(m, file("${bundle}/${m.id}")) }
    }
    def ch_boltz_complex_pdb = ch_complex_per_design.map { meta, d -> tuple(meta, file("${d}/${meta.id}_complex.pdb")) }

    // Run Boltz binder monomer prediction with RMSD analysis
    ch_monomer_batches = ch_filtered_for_refold
        .join(ch_boltz_complex_pdb)
        .collate(boltz_batch_int, true)
        .map { batch ->
            tuple(
                batch.collect { it[0] },
                batch.collect { it[1] },
                batch.collect { it[2] },
                batch.collect { it[1].extension },
            )
        }

    BOLTZ_COMPARE_BINDER_MONOMER(
        ch_monomer_batches,
        binder_chain,
    )

    def ch_monomer_per_design = BOLTZ_COMPARE_BINDER_MONOMER.out.per_design_bundle.flatMap { metas, bundle ->
        def ml = metas instanceof List ? metas : [metas]
        ml.collect { m -> tuple(m, file("${bundle}/${m.id}")) }
    }
    def ch_boltz_monomer_pdb = ch_monomer_per_design.map { meta, d -> tuple(meta, file("${d}/${meta.id}_monomer.pdb")) }

    // Aggregate RMSD outputs
    ch_target_aligned_rmsd = ch_complex_per_design
        .map { meta, d -> file("${d}/rmsd_target_aligned_binder.tsv") }
        .collectFile(
            name: 'rmsd_target_aligned_binder.tsv',
            storeDir: "${outdir}/boltz_refold/rmsd",
            keepHeader: true,
            skip: 1,
        )

    ch_complex_rmsd = ch_complex_per_design
        .map { meta, d -> file("${d}/rmsd_complex.tsv") }
        .collectFile(
            name: "rmsd_complex_vs_${rmsd_label}.tsv",
            storeDir: "${outdir}/boltz_refold/rmsd",
            keepHeader: true,
            skip: 1,
        )

    ch_monomer_vs_design_rmsd = ch_monomer_per_design
        .map { meta, d -> file("${d}/rmsd_monomer_vs_af2ig.tsv") }
        .collectFile(
            name: "rmsd_monomer_vs_${rmsd_label}.tsv",
            storeDir: "${outdir}/boltz_refold/rmsd",
            keepHeader: true,
            skip: 1,
        )

    ch_monomer_vs_complex_rmsd = ch_monomer_per_design
        .map { meta, d -> file("${d}/rmsd_monomer_vs_complex.tsv") }
        .collectFile(
            name: 'rmsd_monomer_vs_complex.tsv',
            storeDir: "${outdir}/boltz_refold/rmsd",
            keepHeader: true,
            skip: 1,
        )

    // Aggregate confidence outputs
    ch_complex_confidence = ch_complex_per_design
        .map { meta, d -> file("${d}/confidence.tsv") }
        .collectFile(
            name: 'boltz_scores_complex.tsv',
            storeDir: "${outdir}/boltz_refold",
            keepHeader: true,
            skip: 1,
        )

    ch_monomer_confidence = ch_monomer_per_design
        .map { meta, d -> file("${d}/confidence.tsv") }
        .collectFile(
            name: 'boltz_scores_binder_monomer.tsv',
            storeDir: "${outdir}/boltz_refold",
            keepHeader: true,
            skip: 1,
        )

    ch_ipsae_tsv = ch_complex_per_design.map { meta, d -> tuple(meta, file("${d}/ipsae.tsv")) }
    ch_ipsae_byres_tsv = ch_complex_per_design.map { meta, d -> tuple(meta, file("${d}/ipsae_byres.tsv")) }

    // BindCraft-score the Boltz-2 refolded complexes
    BINDCRAFT_SCORING_BOLTZ_COMPLEX(
        ch_boltz_complex_pdb.map { meta, pdb -> pdb },
        binder_chain,
        'default_4stage_multimer',
        bindcraft_publish_subdir,
    )

    ch_extra_scores = BINDCRAFT_SCORING_BOLTZ_COMPLEX.out.scores.collectFile(
        name: 'boltz_complex_extra_scores.tsv',
        storeDir: "${outdir}/boltz_refold",
        keepHeader: true,
        skip: 1,
    )

    emit:
    rmsd_target_aligned_binder = ch_target_aligned_rmsd
    rmsd_complex = ch_complex_rmsd
    rmsd_monomer_vs_design = ch_monomer_vs_design_rmsd
    rmsd_monomer_vs_complex = ch_monomer_vs_complex_rmsd
    boltz_scores_complex = ch_complex_confidence
    boltz_scores_monomer = ch_monomer_confidence
    boltz_extra_scores = ch_extra_scores
    ipsae_tsv = ch_ipsae_tsv
    ipsae_byres_tsv = ch_ipsae_byres_tsv
    boltz_complex_pdb = ch_boltz_complex_pdb
    boltz_monomer_pdb = ch_boltz_monomer_pdb
}
