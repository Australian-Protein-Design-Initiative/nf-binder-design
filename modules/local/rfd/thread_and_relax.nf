// Builds a real 3D structure from ProteinMPNN NIM's sequence-only output,
// reusing the baseline dl_binder_design container's PyRosetta/FastRelax code.
process THREAD_AND_RELAX {
    container 'ghcr.io/australian-protein-design-initiative/containers/proteinmpnn_dl_binder_design:latest'
    publishDir "${params.outdir}/rfd/proteinmpnn_nim_relaxed", pattern: '*.pdb', mode: 'copy'

    input:
    path backbone_pdb
    path fasta

    output:
    path "${fasta.baseName}.pdb", emit: pdbs

    script:
    """
    set -euo pipefail
    export CUDA_VISIBLE_DEVICES=""
    export THREAD_BACKBONE_PDB="${backbone_pdb}"
    export THREAD_FASTA="${fasta}"
    export THREAD_OUTPUT_PDB="${fasta.baseName}.pdb"
    export THREAD_RELAX_XML="/app/dl_binder_design/mpnn_fr/RosettaFastRelaxUtil.xml"

    thread_and_relax.py
    """
}
