#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
NIM-accelerated proof-of-concept variant of the RFD workflow (see workflows/rfd.nf
for the baseline). Only exists to prove Nextflow can drive real NVIDIA NIM
containers (RFdiffusion NIM, ProteinMPNN NIM) as ordinary AWS Batch tasks -
start container, do one unit of work, exit - same shape as every other step in
this pipeline, rather than as a separately-managed standing service.

Scope: RFDIFFUSION_NIM -> DL_BINDER_DESIGN_PROTEINMPNN_NIM -> THREAD_AND_RELAX.
Stops before AF2 scoring - no AlphaFold2 NIM equivalent for the baseline's
"initial guess" trick yet (separate, not-yet-made decision).

Usage:
  nextflow run main.nf --method rfd_nim --input_pdb target.pdb --rfd_n_designs=4 \
    --pmpnn_seqs_per_struct=2 -profile aws_batch_nims
*/

params.input_pdb = false
params.outdir = 'results'
params.contigs = ''
params.hotspot_res = false
params.rfd_n_designs = 2
params.pmpnn_seqs_per_struct = 1
params.pmpnn_temperature = 0.000001

include { UNIQUE_ID } from '../modules/local/common/unique_id'
include { RFDIFFUSION_NIM } from '../modules/local/rfd/rfdiffusion_nim'
include { DL_BINDER_DESIGN_PROTEINMPNN_NIM } from '../modules/local/rfd/dl_binder_design_nim'
include { THREAD_AND_RELAX } from '../modules/local/rfd/thread_and_relax'

workflow RFD_NIM {

    main:

    if (params.input_pdb == false) {
        log.info(
            """
        ==================================================================
        PROTEIN BINDER DESIGN PIPELINE - RFDiffusion NIM proof of concept
        ==================================================================
        Covers RFdiffusion + ProteinMPNN via their NVIDIA NIM containers, then
        threads/relaxes the designed sequence onto the backbone. Stops before
        AF2 scoring - see this file's header comment for why.

        Required arguments:
            --input_pdb           Input PDB file for the target

        Optional arguments:
            --outdir              Output directory [default: ${params.outdir}]
            --contigs             Contig map for RFdiffusion [default: ${params.contigs}]
            --hotspot_res         Hotspot residues, eg "A56" - chain ID required [default: ${params.hotspot_res}]
            --rfd_n_designs       Number of RFdiffusion designs [default: ${params.rfd_n_designs}]
            --pmpnn_seqs_per_struct Number of ProteinMPNN sequences per backbone [default: ${params.pmpnn_seqs_per_struct}]
            --pmpnn_temperature   Sampling temperature for ProteinMPNN [default: ${params.pmpnn_temperature}]
        """.stripIndent()
        )
        exit(1)
    }

    if (!System.getenv('NGC_API_KEY')) {
        log.error("NGC_API_KEY is not set in the environment. Both NIM containers download model weights from NVIDIA's NGC API on startup and need this to authenticate - export it locally before running (e.g. `export NGC_API_KEY=...`), not as a --param (params get written to params.json in the output bucket).")
        exit(1)
    }

    UNIQUE_ID()
    ch_unique_id = UNIQUE_ID.out.id_file.map { it.text.trim() }

    ch_input_pdb = Channel.fromPath(params.input_pdb).first()

    def hotspot_res = params.hotspot_res
    if (params.hotspot_res) {
        hotspot_res = "[${params.hotspot_res.trim().replaceAll(/^\[+/, '').replaceAll(/\]+\$/, '')}]"
    }

    ch_design_index = Channel.of(0..(params.rfd_n_designs - 1))

    RFDIFFUSION_NIM(
        ch_input_pdb,
        params.contigs,
        hotspot_res,
        ch_design_index,
        ch_unique_id,
    )

    ch_pmpnn_inputs = RFDIFFUSION_NIM.out.pdbs.flatten()
        | combine(Channel.of(0..(params.pmpnn_seqs_per_struct - 1)))

    DL_BINDER_DESIGN_PROTEINMPNN_NIM(
        ch_pmpnn_inputs.map { pdb, idx -> pdb },
        'A',
        params.pmpnn_temperature,
        ch_pmpnn_inputs.map { pdb, idx -> idx },
    )

    THREAD_AND_RELAX(
        DL_BINDER_DESIGN_PROTEINMPNN_NIM.out.backbone,
        DL_BINDER_DESIGN_PROTEINMPNN_NIM.out.fasta,
    )

    emit:
    backbones = RFDIFFUSION_NIM.out.pdbs
    sequences = DL_BINDER_DESIGN_PROTEINMPNN_NIM.out.fasta
    threaded_pdbs = THREAD_AND_RELAX.out.pdbs
}
