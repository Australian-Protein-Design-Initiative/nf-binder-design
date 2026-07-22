nextflow.enable.dsl = 2

// Copy a single structure under a collision-free basename so a later
// groupTuple can stage many files into one ENGENS work dir. Nextflow 24.04
// rejects duplicate input basenames (e.g. multiple *_model_0.cif from
// different batches/engines).
process ENGENS_STAGE_STRUCTURE {
    tag "${newname}"

    input:
    tuple val(id), val(newname), path(structure)

    output:
    tuple val(id), path(newname), emit: staged

    script:
    """
    cp -L -- "${structure}" "${newname}"
    """
}
