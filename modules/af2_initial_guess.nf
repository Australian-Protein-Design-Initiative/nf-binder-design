process AF2_INITIAL_GUESS {
    container 'ghcr.io/australian-protein-design-initiative/containers/af2_initial_guess:nv-cuda12'

    publishDir "${params.outdir}/af2_initial_guess", mode: 'copy'

    input:
    path 'input/*'

    output:
    path 'pdbs/*.pdb', emit: pdbs
    path 'scores/*.cs', emit: scores

    script:
    """
    mkdir -p scores/

    # Find least-used GPU (by active processes and VRAM) and set CUDA_VISIBLE_DEVICES
    if [[ -n "${params.gpu_devices}" ]]; then
        free_gpu=\$(${baseDir}/bin/find_available_gpu.py "${params.gpu_devices}" --verbose --exclude "${params.gpu_allocation_detect_process_regex}" --random-wait 2)
        export CUDA_VISIBLE_DEVICES="\$free_gpu"
        echo "Set CUDA_VISIBLE_DEVICES=\$free_gpu"
    fi

    # Get first input PDB filename without extension
    PREFIX=\$(ls input/*.pdb | head -n1 | xargs basename | sed 's/\\.pdb\$//')

    python /app/dl_binder_design/af2_initial_guess/predict.py \
        -pdbdir input/ \
        -outpdbdir pdbs/ \
        -recycle ${params.af2ig_recycle} \
        -scorefilename scores/\${PREFIX}.scores.cs
    """
    /*
usage: predict.py [-h] [-pdbdir PDBDIR] [-silent SILENT] [-outpdbdir OUTPDBDIR] [-outsilent OUTSILENT] [-runlist RUNLIST] [-checkpoint_name CHECKPOINT_NAME] [-scorefilename SCOREFILENAME]
                  [-maintain_res_numbering] [-debug] [-max_amide_dist MAX_AMIDE_DIST] [-recycle RECYCLE] [-no_initial_guess] [-force_monomer]

optional arguments:
  -h, --help            show this help message and exit
  -pdbdir PDBDIR        The name of a directory of pdbs to run through the model
  -silent SILENT        The name of a silent file to run through the model
  -outpdbdir OUTPDBDIR  The directory to which the output PDB files will be written. Only used when -pdbdir is active
  -outsilent OUTSILENT  The name of the silent file to which output structs will be written. Only used when -silent is active
  -runlist RUNLIST      The path of a list of pdb tags to run. Only used when -pdbdir is active (default: ''; Run all PDBs)
  -checkpoint_name CHECKPOINT_NAME
                        The name of a file where tags which have finished will be written (default: check.point)
  -scorefilename SCOREFILENAME
                        The name of a file where scores will be written (default: out.sc)
  -maintain_res_numbering
                        When active, the model will not renumber the residues when bad inputs are encountered (default: False)
  -debug                When active, errors will cause the script to crash and the error message to be printed out (default: False)
  -max_amide_dist MAX_AMIDE_DIST
                        The maximum distance between an amide bond's carbon and nitrogen (default: 3.0)
  -recycle RECYCLE      The number of AF2 recycles to perform (default: 3)
  -no_initial_guess     When active, the model will not use an initial guess (default: False)
  -force_monomer        When active, the model will predict the structure of a monomer (default: False)

    */
}
