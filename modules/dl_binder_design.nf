process DL_BINDER_DESIGN_PROTEINMPNN {
    publishDir "${params.outdir}/proteinmpnn", mode: 'copy'
    
    input:
    path 'input/*'
    val relax_cycles
    val seqs_per_struct
    
    output:
    path "pdbs/*", emit: pdbs
    
    script:
    """
    /app/dl_binder_design/mpnn_fr/dl_interface_design.py \
        -pdbdir input/ \
        -relax_cycles ${relax_cycles} \
        -seqs_per_struct ${seqs_per_struct} \
        -outpdbdir pdbs/
    """

/*
usage: dl_interface_design.py [-h] [-pdbdir PDBDIR] [-silent SILENT] [-outpdbdir OUTPDBDIR] [-outsilent OUTSILENT] [-runlist RUNLIST] [-checkpoint_name CHECKPOINT_NAME] [-debug]
                            [-relax_cycles RELAX_CYCLES] [-output_intermediates] [-seqs_per_struct SEQS_PER_STRUCT] [-checkpoint_path CHECKPOINT_PATH] [-temperature TEMPERATURE]
                            [-augment_eps AUGMENT_EPS] [-protein_features PROTEIN_FEATURES] [-omit_AAs OMIT_AAS] [-bias_AA_jsonl BIAS_AA_JSONL] [-num_connections NUM_CONNECTIONS]

options:
  -h, --help            show this help message and exit
  -pdbdir PDBDIR        The name of a directory of pdbs to run through the model
  -silent SILENT        The name of a silent file to run through the model
  -outpdbdir OUTPDBDIR  The directory to which the output PDB files will be written, used if the -pdbdir arg is active
  -outsilent OUTSILENT  The name of the silent file to which output structs will be written, used if the -silent arg is active
  -runlist RUNLIST      The path of a list of pdb tags to run, only active when the -pdbdir arg is active (default: ''; Run all PDBs)
  -checkpoint_name CHECKPOINT_NAME
                        The name of a file where tags which have finished will be written (default: check.point)
  -debug                When active, errors will cause the script to crash and the error message to be printed out (default: False)
  -relax_cycles RELAX_CYCLES
                        The number of relax cycles to perform on each structure (default: 1)
  -output_intermediates
                        Whether to write all intermediate sequences from the relax cycles to disk (default: False)
  -seqs_per_struct SEQS_PER_STRUCT
                        The number of sequences to generate for each structure (default: 1)
  -checkpoint_path CHECKPOINT_PATH
                        The path to the ProteinMPNN weights you wish to use, default /app/dl_binder_design/mpnn_fr/ProteinMPNN/vanilla_model_weights/v_48_020.pt
  -temperature TEMPERATURE
                        The sampling temperature to use when running ProteinMPNN (default: 0.000001)
  -augment_eps AUGMENT_EPS
                        The variance of random noise to add to the atomic coordinates (default 0)
  -protein_features PROTEIN_FEATURES
                        What type of protein features to input to ProteinMPNN (default: full)
  -omit_AAs OMIT_AAS    A string of all residue types (one letter case-insensitive) that you would not like to use for design. Letters not corresponding to residue types will be ignored
                        (default: CX)
  -bias_AA_jsonl BIAS_AA_JSONL
                        The path to a JSON file containing a dictionary mapping residue one-letter names to the bias for that residue eg. {A: -1.1, F: 0.7} (default: ; no bias)
  -num_connections NUM_CONNECTIONS
                        Number of neighbors each residue is connected to. Do not mess around with this argument unless you have a specific set of ProteinMPNN weights which expects a
                        different number of connections. (default: 48)
*/
} 