// Extract the target and per-design binder sequences from a set of design complexes,
// as the targets.fasta / binders.fasta pair FOLD_PULLDOWN_CORE takes. One task for
// the whole set, so the target is read once and every design is checked against it.
//
// With --refold_target_fasta the target comes from that file instead, which is how a
// design made against a truncated construct is refolded against the full-length
// protein. The designs' own target is still read, so the script can say whether the
// supplied one contains it.
process REFOLD_PREPARE_INPUTS {
    tag "refold_inputs"

    container 'ghcr.io/australian-protein-design-initiative/containers/nf-binder-design-utils:0.1.8'

    publishDir path: "${params.outdir}/${params.refold_publish_dir ?: 'refold'}/input", mode: 'copy'

    input:
    // One directory per design rather than one shared one: not every workflow gives its
    // designs distinct filenames. RFdiffusion3 writes every design to model.cif inside
    // its own directory, so staging them side by side fails the task outright with
    // "multiple input files for each of the following file names: designs/model.cif",
    // and anything that deduplicated the names silently would make both designs answer
    // to the id "model".
    path designs, stageAs: 'staged/d*/'
    val design_ids
    val binder_chain
    val target_chain
    val target_id
    // The target to fold against, when it is not the one in the designs. The
    // empty placeholder means "read it back out of the designs", as before.
    path target_fasta

    output:
    path 'targets.fasta', emit: targets
    path 'binders.fasta', emit: binders

    script:
    // A one-element input arrives as a single Path and a single String, not as lists.
    def design_files = (designs instanceof List) ? designs : [designs]
    def ids = (design_ids instanceof List) ? design_ids : [design_ids]
    if (design_files.size() != ids.size()) {
        throw new IllegalArgumentException(
            "REFOLD_PREPARE_INPUTS got ${design_files.size()} designs and ${ids.size()} ids."
        )
    }
    // prepare_inputs.py names each record after the file it read, so the link name is
    // what ends up in binders.fasta and in every score table downstream. The extension
    // has to survive because that is how the structure's format is recognised.
    def links = [ids, design_files].transpose().collect { id, design ->
        def name = design.name
        def ext = ['.cif.gz', '.pdb.gz', '.cif', '.pdb'].find { name.endsWith(it) } ?: ''
        "ln -s '../${design}' 'designs/${id}${ext}'"
    }.join('\n    ')
    def target_fasta_arg = target_fasta.name == 'no_refold_target.fasta' ? '' : "--target-fasta '${target_fasta}'"
    """
    set -euo pipefail

    mkdir -p designs
    ${links}

    python ${projectDir}/bin/refold/prepare_inputs.py designs/ \\
        --binder-chain '${binder_chain}' \\
        --target-chain '${target_chain}' \\
        --target-id '${target_id}' \\
        ${target_fasta_arg} \\
        --targets-out targets.fasta \\
        --binders-out binders.fasta
    """
}
