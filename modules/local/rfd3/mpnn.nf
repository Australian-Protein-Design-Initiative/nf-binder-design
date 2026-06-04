process MPNN {
    container 'ghcr.io/australian-protein-design-initiative/containers/rc-foundry:0.1.12-weights'

    publishDir path: "${params.outdir}/rfd3/mpnn", pattern: 'output/*.cif', mode: 'copy'
    publishDir path: "${params.outdir}/rfd3/mpnn", pattern: 'output/*.fa', mode: 'copy'

    input:
    path structure_file
    val mpnn_args

    output:
    path 'output/*.cif', emit: cifs
    path 'output/*.fa', emit: fastas

    script:
    """
    set -euo pipefail

    mkdir -p output

    mpnn \
        --structure_path "${structure_file}" \
        --out_directory output \
        ${mpnn_args} \
        ${task.ext.args ?: ''}
    
    ####
    # Fix unquoted [] in CIF files output by mpnn (v0.1.11) so that ChimeraX can read them
    # Issue: https://github.com/RosettaCommons/foundry/issues/128
    ####
    for f in output/*.cif; do
        [ -f "\${f}" ] || continue
        sed -i "s/\\([[:space:]]\\)\\[\\]/\\1'[]'/g" "\${f}"
    done
    """
}

/* 
Usage:

https://rosettacommons.github.io/foundry/models/mpnn/index.html
https://github.com/RosettaCommons/foundry/tree/production/models/mpnn

$ mpnn --help

Environment variable CCD_MIRROR_PATH not set. Will not be able to use function requiring this variable. To set it you may:
  (1) add the line 'export VAR_NAME=path/to/variable' to your .bashrc or .zshrc file
  (2) set it in your current shell with 'export VAR_NAME=path/to/variable'
  (3) write it to a .env file in the root of the atomworks.io repository
Environment variable PDB_MIRROR_PATH not set. Will not be able to use function requiring this variable. To set it you may:
  (1) add the line 'export VAR_NAME=path/to/variable' to your .bashrc or .zshrc file
  (2) set it in your current shell with 'export VAR_NAME=path/to/variable'
  (3) write it to a .env file in the root of the atomworks.io repository
usage: mpnn [-h] [--config_json CONFIG_JSON] [--model_type {protein_mpnn,ligand_mpnn}] [--checkpoint_path CHECKPOINT_PATH] [--is_legacy_weights {True,False}]
            [--out_directory OUT_DIRECTORY] [--write_fasta {True,False}] [--write_structures {True,False}] [--structure_path STRUCTURE_PATH] [--name NAME]
            [--seed SEED] [--batch_size BATCH_SIZE] [--number_of_batches NUMBER_OF_BATCHES] [--remove_ccds REMOVE_CCDS] [--remove_waters {True,False,None}]
            [--occupancy_threshold_sidechain OCCUPANCY_THRESHOLD_SIDECHAIN] [--occupancy_threshold_backbone OCCUPANCY_THRESHOLD_BACKBONE]
            [--undesired_res_names UNDESIRED_RES_NAMES] [--structure_noise STRUCTURE_NOISE] [--decode_type {auto_regressive,teacher_forcing}]
            [--causality_pattern {auto_regressive,unconditional,conditional,conditional_minus_self}]
            [--initialize_sequence_embedding_with_ground_truth {True,False}] [--features_to_return FEATURES_TO_RETURN] [--atomize_side_chains {True,False}]
            [--fixed_residues FIXED_RESIDUES | --designed_residues DESIGNED_RESIDUES | --fixed_chains FIXED_CHAINS | --designed_chains DESIGNED_CHAINS]
            [--bias BIAS] [--bias_per_residue BIAS_PER_RESIDUE] [--omit OMIT] [--omit_per_residue OMIT_PER_RESIDUE] [--pair_bias PAIR_BIAS]
            [--pair_bias_per_residue_pair PAIR_BIAS_PER_RESIDUE_PAIR] [--temperature TEMPERATURE] [--temperature_per_residue TEMPERATURE_PER_RESIDUE]
            [--symmetry_residues SYMMETRY_RESIDUES | --homo_oligomer_chains HOMO_OLIGOMER_CHAINS] [--symmetry_residues_weights SYMMETRY_RESIDUES_WEIGHTS]

MPNN JSON-driven inference CLI

options:
  -h, --help            show this help message and exit
  --config_json CONFIG_JSON
                        Path to existing JSON config file. When provided, all other CLI flags are parsed but ignored. (default: None)
  --model_type {protein_mpnn,ligand_mpnn}
                        Model type to use. (default: None)
  --checkpoint_path CHECKPOINT_PATH
                        Path to model checkpoint. (default: None)
  --is_legacy_weights {True,False}
                        Whether to interpret checkpoint as legacy-weight ordering. (default: None)
  --out_directory OUT_DIRECTORY
                        Output directory for CIF/FASTA. (default: None)
  --write_fasta {True,False}
                        Whether to write FASTA outputs. (default: True)
  --write_structures {True,False}
                        Whether to write designed structures (CIF). (default: True)
  --structure_path STRUCTURE_PATH
                        Path to structure file (CIF or PDB). (default: None)
  --name NAME           Optional name / label for the input. (default: None)
  --seed SEED           Random seed for sampling. (default: None)
  --batch_size BATCH_SIZE
                        Batch size for sampling. At inference, this also controls the effective repeat_sample_num passed to the pipeline. (default: 1)
  --number_of_batches NUMBER_OF_BATCHES
                        Number of batches of size batch_size to draw. (default: 1)
  --remove_ccds REMOVE_CCDS
                        Comma-separated list of CCD residue names to remove as solvents/crystallization components during parsing (overrides
                        STANDARD_PARSER_ARGS). 'None' has special behavior: use the parser default behavior. (default: [])
  --remove_waters {True,False,None}
                        If set, override the parser default for removing water-like residues (overrides STANDARD_PARSER_ARGS). 'None' has special behavior:
                        use the parser default behavior. (default: None)
  --occupancy_threshold_sidechain OCCUPANCY_THRESHOLD_SIDECHAIN
                        Sidechain occupancy threshold used in the MPNN pipeline. 'None' has special behavior: use the pipeline default behavior. (default:
                        0.0)
  --occupancy_threshold_backbone OCCUPANCY_THRESHOLD_BACKBONE
                        Backbone occupancy threshold used in the MPNN pipeline. 'None' has special behavior: use the pipeline default behavior. (default: 0.0)
  --undesired_res_names UNDESIRED_RES_NAMES
                        JSON or comma-separated list of residue names to treat as undesired in the pipeline. 'None' has special behavior: use the pipeline
                        default behavior. (default: [])
  --structure_noise STRUCTURE_NOISE
                        Structure noise (Angstroms) used in user settings. (default: 0.0)
  --decode_type {auto_regressive,teacher_forcing}
                        Decoding type for MPNN inference. - auto_regressive: use previously predicted residues for all previous positions when predicting each
                        residue. This is the default for inference. - teacher_forcing: use ground-truth residues from the structure for all previous positions
                        when predicting each residue. (default: auto_regressive)
  --causality_pattern {auto_regressive,unconditional,conditional,conditional_minus_self}
                        Causality pattern for decoding. - auto_regressive: each position attends to the sequence and decoder representation of all previously
                        decoded positions. This is the default for inference. - unconditional: each position does not attend to the sequence or decoder
                        representation of any other positions (encoder representations only). - conditional: each position attends to the sequence and decoder
                        representation of all other positions. - conditional_minus_self: each position attends to the sequence and decoder representation of
                        all other positions, except for itself (as a destination node). (default: auto_regressive)
  --initialize_sequence_embedding_with_ground_truth {True,False}
                        Whether to initialize the sequence embedding with ground truth residues from the input structure. - False: initialize the sequence
                        embedding with zeros. If doing auto-regressive decoding, initialize S_sampled with unknown residues. This is the default for
                        inference. - True: initialize the sequence embedding with the ground truth sequence from the input structure. If doing auto-regressive
                        decoding, also initialize S_sampled with the ground truth. This affects the pair bias application. (default: False)
  --features_to_return FEATURES_TO_RETURN
                        JSON dict for features_to_return; e.g. '{"input_features": ["mask_for_loss"], "decoder_features": ["log_probs"]}' (default: None)
  --atomize_side_chains {True,False}
                        Whether to atomize side chains of fixed residues. Only applicable for LigandMPNN. (default: False)
  --fixed_residues FIXED_RESIDUES
                        List of residue IDs to fix: e.g. '["A35","B40","C52"]' or "A35,B40,C52" (default: None)
  --designed_residues DESIGNED_RESIDUES
                        List of residue IDs to design: e.g. '["A35","B40","C52"]' or "A35,B40,C52" (default: None)
  --fixed_chains FIXED_CHAINS
                        List of chain IDs to fix: e.g. '["A","B"]' or "A,B" (default: None)
  --designed_chains DESIGNED_CHAINS
                        List of chain IDs to design: e.g. '["A","B"]' or "A,B" (default: None)
  --bias BIAS           Bias dict: e.g. '{"ALA": -1.0, "GLY": 0.5}' (default: None)
  --bias_per_residue BIAS_PER_RESIDUE
                        Per-residue bias dict: e.g. '{"A35": {"ALA": -2.0}}'. Overwrites --bias. (default: None)
  --omit OMIT           List of residue types to omit: e.g. '["ALA","GLY","UNK"]'. (default: ['UNK'])
  --omit_per_residue OMIT_PER_RESIDUE
                        Per-residue list of residue types to omit: e.g. '{"A35": ["ALA","GLY","UNK"]}'. Overwrites --omit. (default: None)
  --pair_bias PAIR_BIAS
                        Controls the bias applied due to residue selections at neighboring positions: '{"ALA": {"GLY": -0.5}, "GLY": {"ALA": -0.5}}' (default:
                        None)
  --pair_bias_per_residue_pair PAIR_BIAS_PER_RESIDUE_PAIR
                        Per-residue-pair dict for controlling bias due to residue selections at neighboring positions: '{"A35": {"B40": {"ALA": {"GLY":
                        -1.0}}}}' . Overwrites --pair_bias. Note that this is NOT applied symmetrically; if the outer residue ID corresponds to the first
                        token; the inner residue ID corresponds to the second token. This should be read as follows: for residue pair (i,j) (e.g.
                        ("A35","B40")), the inner dictionaries dictate that if residue i is assigned as the first token (e.g. "ALA"), then the bias for
                        assigning residue j is the innermost dict (e.g. {"GLY": -1.0} ). (default: None)
  --temperature TEMPERATURE
                        Temperature for sampling. (default: 0.1)
  --temperature_per_residue TEMPERATURE_PER_RESIDUE
                        Per-residue temperature dict: e.g. '{"A35": 0.1}'. Overwrites --temperature. (default: None)
  --symmetry_residues SYMMETRY_RESIDUES
                        Residue-based symmetry groups, each a list of residue IDs. Example: '[["A35","B35"],["A40","B40","C40"]]' (default: None)
  --homo_oligomer_chains HOMO_OLIGOMER_CHAINS
                        Homo-oligomer chain groups, each a list of chain IDs. Within each group, chains must have the same number of residues in the same
                        order; residues at matching positions across chains are treated as symmetry-equivalent. Example: '[["A","B","C"]]' (default: None)
  --symmetry_residues_weights SYMMETRY_RESIDUES_WEIGHTS
                        Optional list of symmetry weights matching the shape of symmetry_residues. Example: '[[1.0, 1.0], [1.0, 0.5, -0.5]]'. Ignored if
                        homo_oligomer_chains is used. (default: None)
*/
