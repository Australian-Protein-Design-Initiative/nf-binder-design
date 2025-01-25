process SILENT_FROM_PDBS {
    publishDir "${params.outdir}/silent", mode: 'copy'
    
    input:
    path '*.pdb'
    
    output:
    path "designs.silent", emit: silent
    
    script:
    """
    PATH=/app/dl_binder_design/include/silent_tools:${PATH}
    silentfrompdbs *.pdb > designs.silent
    """
}