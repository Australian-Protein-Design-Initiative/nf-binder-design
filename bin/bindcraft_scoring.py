# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "biopython",
#     "numpy",
#     "pyrosetta-installer",
#     "scipy",
# ]
# ///

##
# Based on code from https://github.com/martinpacesa/BindCraft
# MIT License, 2024 Martin Pacesa
##

# Run using a container like:
#  apptainer exec /scratch/projects/apdi/containers/bindcraft-nv-cuda12.sif \
#    /opt/conda/envs/BindCraft/bin/python bin/bindcraft_scoring.py \
#    /path/to/my_binder.pdb \
#    --binder-chain A \
#    --omit-aas CX \
#    --dalphaball-path /app/BindCraft/functions/DAlphaBall.gcc \
#    --dssp-path /app/BindCraft/functions/dssp

# This doesn't seem to work.
# try:
#     import pyrosetta
# except ImportError:
#     import pyrosetta_installer
#     pyrosetta_installer.install_pyrosetta(skip_if_installed=True)

import os
import math
import argparse
import json
import shutil
import sys
import contextlib
import logging
import numpy as np
import pyrosetta as pr
from collections import defaultdict, OrderedDict
from scipy.spatial import cKDTree
from Bio.PDB import (
    Selection,
    Polypeptide,
    Chain,
)
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.Superimposer import Superimposer
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Polypeptide import is_aa

from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.io import pose_from_pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects


# clean unnecessary rosetta information from PDB
def clean_pdb(pdb_file):
    # Read the pdb file and filter relevant lines
    with open(pdb_file, "r") as f_in:
        relevant_lines = [
            line
            for line in f_in
            if line.startswith(("ATOM", "HETATM", "MODEL", "TER", "END", "LINK"))
        ]

    # Write the cleaned lines back to the original pdb file
    with open(pdb_file, "w") as f_out:
        f_out.writelines(relevant_lines)


three_to_one_map = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


# identify interacting residues at the binder interface
def interface_residues(
    trajectory_pdb, binder_chain="A", target_chains_str="", atom_distance_cutoff=4.0
):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", trajectory_pdb)
    model = structure[0]

    # Check if binder chain exists
    if binder_chain not in model:
        logging.warning(
            f"Binder chain '{binder_chain}' not found in '{trajectory_pdb}'. Cannot determine interface residues."
        )
        return {}

    # Get the specified chain's atoms
    binder_atoms = Selection.unfold_entities(model[binder_chain], "A")
    binder_coords = np.array([atom.coord for atom in binder_atoms])

    # Get atoms for the target chains
    target_atoms = []
    if target_chains_str:
        target_chain_ids = target_chains_str.split(",")
        for chain_id in target_chain_ids:
            if chain_id in model:
                target_atoms.extend(Selection.unfold_entities(model[chain_id], "A"))
            else:
                logging.warning(
                    f"Target chain '{chain_id}' not found in '{trajectory_pdb}'."
                )

    if not target_atoms:
        logging.warning(
            f"No target atoms found in '{trajectory_pdb}' for chains '{target_chains_str}'."
        )
        return {}
    target_coords = np.array([atom.coord for atom in target_atoms])

    # Build KD trees for both sets of atoms
    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    # Prepare to collect interacting residues
    interacting_residues = {}

    # Query the tree for pairs of atoms within the distance cutoff
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

    # Process each binder atom's interactions
    for binder_idx, close_indices in enumerate(pairs):
        binder_residue = binder_atoms[binder_idx].get_parent()
        binder_resname = binder_residue.get_resname()

        # Convert three-letter code to single-letter code using the manual dictionary
        if binder_resname in three_to_one_map:
            aa_single_letter = three_to_one_map[binder_resname]
            for close_idx in close_indices:
                target_residue = target_atoms[close_idx].get_parent()
                interacting_residues[binder_residue.id[1]] = aa_single_letter

    return interacting_residues


# Rosetta interface scores
def score_interface(pdb_file, binder_chain="A", target_chains_str=""):
    # load pose
    pose = pr.pose_from_pdb(pdb_file)

    # analyze interface statistics
    iam = InterfaceAnalyzerMover()
    if target_chains_str:
        interface_str = f"{binder_chain}_{'_'.join(target_chains_str.split(','))}"
        iam.set_interface(interface_str)
    else:
        logging.warning("No target chains provided for interface scoring.")
        # Fallback or default behavior if needed, e.g., iam.set_interface("A_B")
        # For now, we'll just log and the mover might fail or do nothing.

    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)

    # Initialize dictionary with all amino acids
    interface_AA = {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWY"}

    # Initialize list to store PDB residue IDs at the interface
    interface_residues_set = interface_residues(
        pdb_file, binder_chain, target_chains_str=target_chains_str
    )
    interface_residues_pdb_ids = []

    # Iterate over the interface residues
    for pdb_res_num, aa_type in interface_residues_set.items():
        # Increase the count for this amino acid type
        interface_AA[aa_type] += 1

        # Append the binder_chain and the PDB residue number to the list
        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")

    # count interface residues
    interface_nres = len(interface_residues_pdb_ids)

    # Convert the list into a comma-separated string
    interface_residues_pdb_ids_str = ",".join(interface_residues_pdb_ids)

    # Calculate the percentage of hydrophobic residues at the interface of the binder
    hydrophobic_aa = set("ACFILMPVWY")
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
    if interface_nres != 0:
        interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100
    else:
        interface_hydrophobicity = 0

    # retrieve statistics
    interfacescore = iam.get_all_data()
    interface_sc = interfacescore.sc_value  # shape complementarity
    interface_interface_hbonds = (
        interfacescore.interface_hbonds
    )  # number of interface H-bonds
    interface_dG = iam.get_interface_dG()  # interface dG
    interface_dSASA = (
        iam.get_interface_delta_sasa()
    )  # interface dSASA (interface surface area)
    interface_packstat = iam.get_interface_packstat()  # interface pack stat score
    interface_dG_SASA_ratio = (
        interfacescore.dG_dSASA_ratio * 100
    )  # ratio of dG/dSASA (normalised energy for interface area size)
    buns_filter = XmlObjects.static_get_filter(
        '<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />'
    )
    interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    if interface_nres != 0:
        interface_hbond_percentage = (
            interface_interface_hbonds / interface_nres
        ) * 100  # Hbonds per interface size percentage
        interface_bunsch_percentage = (
            interface_delta_unsat_hbonds / interface_nres
        ) * 100  # Unsaturated H-bonds per percentage
    else:
        interface_hbond_percentage = None
        interface_bunsch_percentage = None

    # calculate binder energy score
    chain_design = ChainSelector(binder_chain)
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_design)
    binder_score = tem.calculate(pose)

    # calculate binder SASA fraction
    bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(chain_design)
    binder_sasa = bsasa.calculate(pose)

    if binder_sasa > 0:
        interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
    else:
        interface_binder_fraction = 0

    # calculate surface hydrophobicity
    binder_pose = {
        pose.pdb_info().chain(pose.conformation().chain_begin(i)): p
        for i, p in zip(range(1, pose.num_chains() + 1), pose.split_by_chain())
    }[binder_chain]

    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
    surface_res = layer_sel.apply(binder_pose)

    exp_apol_count = 0
    total_count = 0

    # count apolar and aromatic residues at the surface
    for i in range(1, len(surface_res) + 1):
        if surface_res[i] == True:
            res = binder_pose.residue(i)

            # count apolar and aromatic residues as hydrophobic
            if (
                res.is_apolar() == True
                or res.name() == "PHE"
                or res.name() == "TRP"
                or res.name() == "TYR"
            ):
                exp_apol_count += 1
            total_count += 1

    surface_hydrophobicity = exp_apol_count / total_count

    # output interface score array and amino acid counts at the interface
    interface_scores = OrderedDict(
        {
            "binder_score": binder_score,
            "surface_hydrophobicity": surface_hydrophobicity,
            "interface_sc": interface_sc,
            "interface_packstat": interface_packstat,
            "interface_dG": interface_dG,
            "interface_dSASA": interface_dSASA,
            "interface_dG_SASA_ratio": interface_dG_SASA_ratio,
            "interface_fraction": interface_binder_fraction,
            "interface_hydrophobicity": interface_hydrophobicity,
            "interface_nres": interface_nres,
            "interface_interface_hbonds": interface_interface_hbonds,
            "interface_hbond_percentage": interface_hbond_percentage,
            "interface_delta_unsat_hbonds": interface_delta_unsat_hbonds,
            "interface_delta_unsat_hbonds_percentage": interface_bunsch_percentage,
        }
    )

    # round to two decimal places
    interface_scores = OrderedDict(
        [
            (k, round(v, 2) if isinstance(v, float) else v)
            for k, v in interface_scores.items()
        ]
    )

    return interface_scores, interface_AA, interface_residues_pdb_ids_str


# align pdbs to have same orientation
def align_pdbs(reference_pdb, align_pdb, reference_chain_id, align_chain_id):
    # initiate poses
    reference_pose = pr.pose_from_pdb(reference_pdb)
    align_pose = pr.pose_from_pdb(align_pdb)

    align = AlignChainMover()
    align.pose(reference_pose)

    # If the chain IDs contain commas, split them and only take the first value
    reference_chain_id = reference_chain_id.split(",")[0]
    align_chain_id = align_chain_id.split(",")[0]

    # Get the chain number corresponding to the chain ID in the poses
    reference_chain = pr.rosetta.core.pose.get_chain_id_from_chain(
        reference_chain_id, reference_pose
    )
    align_chain = pr.rosetta.core.pose.get_chain_id_from_chain(
        align_chain_id, align_pose
    )

    align.source_chain(align_chain)
    align.target_chain(reference_chain)
    align.apply(align_pose)

    # Overwrite aligned pdb
    align_pose.dump_pdb(align_pdb)
    clean_pdb(align_pdb)


# calculate the rmsd without alignment
def unaligned_rmsd(reference_pdb, align_pdb, reference_chain_id, align_chain_id):
    reference_pose = pr.pose_from_pdb(reference_pdb)
    align_pose = pr.pose_from_pdb(align_pdb)

    # Define chain selectors for the reference and align chains
    reference_chain_selector = ChainSelector(reference_chain_id)
    align_chain_selector = ChainSelector(align_chain_id)

    # Apply selectors to get residue subsets
    reference_chain_subset = reference_chain_selector.apply(reference_pose)
    align_chain_subset = align_chain_selector.apply(align_pose)

    # Convert subsets to residue index vectors
    reference_residue_indices = get_residues_from_subset(reference_chain_subset)
    align_residue_indices = get_residues_from_subset(align_chain_subset)

    # Create empty subposes
    reference_chain_pose = pr.Pose()
    align_chain_pose = pr.Pose()

    # Fill subposes
    pose_from_pose(reference_chain_pose, reference_pose, reference_residue_indices)
    pose_from_pose(align_chain_pose, align_pose, align_residue_indices)

    # Calculate RMSD using the RMSDMetric
    rmsd_metric = RMSDMetric()
    rmsd_metric.set_comparison_pose(reference_chain_pose)
    rmsd = rmsd_metric.calculate(align_chain_pose)

    return round(rmsd, 2)


# Relax designed structure
def pr_relax(pdb_file, relaxed_pdb_path):
    if not os.path.exists(relaxed_pdb_path):
        # Generate pose
        pose = pr.pose_from_pdb(pdb_file)
        start_pose = pose.clone()

        ### Generate movemaps
        mmf = MoveMap()
        mmf.set_chi(True)  # enable sidechain movement
        mmf.set_bb(
            True
        )  # enable backbone movement, can be disabled to increase speed by 30% but makes metrics look worse on average
        mmf.set_jump(False)  # disable whole chain movement

        # Run FastRelax
        fastrelax = FastRelax()
        scorefxn = pr.get_fa_scorefxn()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf)  # set MoveMap
        fastrelax.max_iter(200)  # default iterations is 2500
        fastrelax.min_type("lbfgs_armijo_nonmonotone")
        fastrelax.constrain_relax_to_start_coords(True)
        fastrelax.apply(pose)

        # Align relaxed structure to original trajectory
        align = AlignChainMover()
        align.source_chain(0)
        align.target_chain(0)
        align.pose(start_pose)
        align.apply(pose)

        # Copy B factors from start_pose to pose
        for resid in range(1, pose.total_residue() + 1):
            if pose.residue(resid).is_protein():
                # Get the B factor of the first heavy atom in the residue
                bfactor = start_pose.pdb_info().bfactor(resid, 1)
                for atom_id in range(1, pose.residue(resid).natoms() + 1):
                    pose.pdb_info().bfactor(resid, atom_id, bfactor)

        # output relaxed and aligned PDB
        pose.dump_pdb(relaxed_pdb_path)
        clean_pdb(relaxed_pdb_path)


def get_sequence_from_pdb(pdb_file, chain_id):
    """
    Extracts the amino acid sequence from a specific chain in a PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_file)
    model = structure[0]
    chain = model[chain_id]

    seq = ""
    for residue in chain:
        if is_aa(residue):
            resname = residue.get_resname()
            if resname in three_to_one_map:
                seq += three_to_one_map[resname]

    return seq


def extinction_coefficient(sequence):
    """
    Calculate the reduced extinction coefficient per 1% solution for a protein sequence.
    Returns the extinction coefficient as a float.
    """
    analysis = ProteinAnalysis(sequence)
    extinction_coefficient_reduced = analysis.molar_extinction_coefficient()[0]
    molecular_weight = round(analysis.molecular_weight() / 1000, 2)
    if molecular_weight == 0:
        return 0.0
    extinction_coefficient_reduced_1 = round(
        extinction_coefficient_reduced / molecular_weight * 0.01, 2
    )
    return extinction_coefficient_reduced_1


# analyze sequence composition of design
def validate_design_sequence(sequence, num_clashes, omit_AAs):
    note_array = []

    # Check if protein contains clashes after relaxation
    if num_clashes > 0:
        note_array.append("Relaxed structure contains clashes.")

    # Check if the sequence contains disallowed amino acids
    if omit_AAs:
        restricted_AAs = list(omit_AAs)
        for restricted_AA in restricted_AAs:
            if restricted_AA in sequence:
                note_array.append("Contains: " + restricted_AA + "!")

    # Use the extinction_coefficient function
    extinction_coefficient_reduced_1 = extinction_coefficient(sequence)

    # Check if the absorption is high enough
    if extinction_coefficient_reduced_1 <= 2:
        note_array.append(
            f"Absorption value is {extinction_coefficient_reduced_1}, consider adding tryptophane to design."
        )

    # Join the notes into a single string
    notes = " ".join(note_array)

    return notes


# temporary function, calculate RMSD of input PDB and trajectory target
def target_pdb_rmsd(
    trajectory_pdb,
    trajectory_target_chains_string,
    reference_pdb,
    reference_target_chains_string,
):
    # Parse the PDB files
    parser = PDBParser(QUIET=True)
    structure_trajectory = parser.get_structure("trajectory", trajectory_pdb)
    structure_reference = parser.get_structure("reference", reference_pdb)

    # Extract residues from specified chains in the reference PDB
    reference_chain_ids = reference_target_chains_string.split(",")
    residues_reference = []
    model_reference = structure_reference[0]
    for chain_id in map(str.strip, reference_chain_ids):
        if chain_id in model_reference:
            for residue in model_reference[chain_id]:
                if is_aa(residue, standard=True):
                    residues_reference.append(residue)
        else:
            logging.warning(
                f"Chain '{chain_id}' not found in reference PDB '{reference_pdb}'"
            )

    # Extract residues from specified chains in the trajectory PDB
    trajectory_chain_ids = trajectory_target_chains_string.split(",")
    residues_trajectory = []
    model_trajectory = structure_trajectory[0]
    for chain_id in map(str.strip, trajectory_chain_ids):
        if chain_id in model_trajectory:
            for residue in model_trajectory[chain_id]:
                if is_aa(residue, standard=True):
                    residues_trajectory.append(residue)
        else:
            logging.warning(
                f"Chain '{chain_id}' not found in trajectory PDB '{trajectory_pdb}'"
            )

    if not residues_reference or not residues_trajectory:
        logging.warning(
            "Could not calculate target RMSD due to missing chains or residues."
        )
        return 0.0

    # Ensure that both structures have the same number of residues for comparison
    min_length = min(len(residues_reference), len(residues_trajectory))
    residues_reference = residues_reference[:min_length]
    residues_trajectory = residues_trajectory[:min_length]

    # Collect CA atoms from the two sets of residues
    atoms_reference = [
        residue["CA"] for residue in residues_reference if "CA" in residue
    ]
    atoms_trajectory = [
        residue["CA"] for residue in residues_trajectory if "CA" in residue
    ]

    if (
        not atoms_reference
        or not atoms_trajectory
        or len(atoms_reference) != len(atoms_trajectory)
    ):
        logging.warning(
            "Cannot calculate RMSD due to missing C-alpha atoms or unequal lengths."
        )
        return 0.0

    # Calculate RMSD using structural alignment
    sup = Superimposer()
    sup.set_atoms(atoms_reference, atoms_trajectory)
    rmsd = sup.rms

    return round(rmsd, 2)


# detect C alpha clashes for deformed trajectories
def calculate_clash_score(pdb_file, threshold=2.4, only_ca=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    atoms = []
    atom_info = []  # Detailed atom info for debugging and processing

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == "H":  # Skip hydrogen atoms
                        continue
                    if only_ca and atom.get_name() != "CA":
                        continue
                    atoms.append(atom.coord)
                    atom_info.append(
                        (chain.id, residue.id[1], atom.get_name(), atom.coord)
                    )

    tree = cKDTree(atoms)
    pairs = tree.query_pairs(threshold)

    valid_pairs = set()
    for i, j in pairs:
        chain_i, res_i, name_i, coord_i = atom_info[i]
        chain_j, res_j, name_j, coord_j = atom_info[j]

        # Exclude clashes within the same residue
        if chain_i == chain_j and res_i == res_j:
            continue

        # Exclude directly sequential residues in the same chain for all atoms
        if chain_i == chain_j and abs(res_i - res_j) == 1:
            continue

        # If calculating sidechain clashes, only consider clashes between different chains
        if not only_ca and chain_i == chain_j:
            continue

        valid_pairs.add((i, j))

    return len(valid_pairs)


# calculate secondary structure percentage of design
def calc_ss_percentage(
    pdb_file, dssp_path, chain_id="B", target_chains_str="", atom_distance_cutoff=4.0
):
    # Parse the structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]  # Consider only the first model in the structure

    # TODO: Can we replace this with MDAnalysis.analysis.dssp.pydssp_numpy
    #       so we don't need a dssp binary?
    # Calculate DSSP for the model
    dssp = DSSP(model, pdb_file, dssp=dssp_path)

    # Prepare to count residues
    ss_counts = defaultdict(int)
    ss_interface_counts = defaultdict(int)
    plddts_interface = []
    plddts_ss = []

    # Get chain and interacting residues once
    chain = model[chain_id]
    interacting_residues = set(
        interface_residues(
            pdb_file, chain_id, target_chains_str, atom_distance_cutoff
        ).keys()
    )

    for residue in chain:
        residue_id = residue.id[1]
        if (chain_id, residue_id) in dssp:
            ss = dssp[(chain_id, residue_id)][2]  # Get the secondary structure
            ss_type = "loop"
            if ss in ["H", "G", "I"]:
                ss_type = "helix"
            elif ss == "E":
                ss_type = "sheet"

            ss_counts[ss_type] += 1

            if ss_type != "loop":
                # calculate secondary structure normalised pLDDT
                avg_plddt_ss = sum(atom.bfactor for atom in residue) / len(residue)
                plddts_ss.append(avg_plddt_ss)

            if residue_id in interacting_residues:
                ss_interface_counts[ss_type] += 1

                # calculate interface pLDDT
                avg_plddt_residue = sum(atom.bfactor for atom in residue) / len(residue)
                plddts_interface.append(avg_plddt_residue)

    # Calculate percentages
    total_residues = sum(ss_counts.values())
    total_interface_residues = sum(ss_interface_counts.values())

    percentages = calculate_percentages(
        total_residues, ss_counts["helix"], ss_counts["sheet"]
    )
    interface_percentages = calculate_percentages(
        total_interface_residues,
        ss_interface_counts["helix"],
        ss_interface_counts["sheet"],
    )

    i_plddt = (
        round(sum(plddts_interface) / len(plddts_interface) / 100, 2)
        if plddts_interface
        else 0
    )
    ss_plddt = round(sum(plddts_ss) / len(plddts_ss) / 100, 2) if plddts_ss else 0

    return (*percentages, *interface_percentages, i_plddt, ss_plddt)


def calculate_percentages(total, helix, sheet):
    helix_percentage = round((helix / total) * 100, 2) if total > 0 else 0
    sheet_percentage = round((sheet / total) * 100, 2) if total > 0 else 0
    loop_percentage = (
        round(((total - helix - sheet) / total) * 100, 2) if total > 0 else 0
    )

    return helix_percentage, sheet_percentage, loop_percentage


# TODO: Can we also run something like:
#    colabdesign.af.model.mk_af_model(protocol="binder",
#                                     use_multimer=True,
#                                     use_initial_guess=True,
#                                     use_initial_atom_pos=True)
# to get the ipTM ? (possibly instead of af2_initial_guess ?)
def score_pdb(
    pdb_file,
    binder_chain,
    dssp_path,
    omit_aas,
    reference_pdb=None,
    reference_target_chains=None,
):
    """
    Calculates a comprehensive set of scores for a given binder-target PDB file.
    """
    scores = OrderedDict()

    # Create a temporary directory for relaxed files
    temp_dir = "temp_relaxed"
    os.makedirs(temp_dir, exist_ok=True)

    base_name = os.path.basename(pdb_file)
    relaxed_pdb_path = os.path.join(temp_dir, base_name.replace(".pdb", "_relaxed.pdb"))

    # 1. Relax the structure
    logging.info(f"Relaxing {pdb_file}...")
    pr_relax(pdb_file, relaxed_pdb_path)

    # Determine target chains in the relaxed model, to be used by other functions
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", relaxed_pdb_path)
    model_chains = [c.id for c in structure[0]]
    model_target_chains_list = [c for c in model_chains if c != binder_chain]
    model_target_chains_str = ",".join(model_target_chains_list)

    # 2. Calculate clash scores
    logging.info("Calculating clash scores...")
    scores["unrelaxed_clashes"] = calculate_clash_score(pdb_file)
    scores["relaxed_clashes"] = calculate_clash_score(relaxed_pdb_path)

    # 3. Score the interface
    logging.info("Scoring interface...")
    interface_scores, interface_AA, interface_residues_str = score_interface(
        relaxed_pdb_path, binder_chain, target_chains_str=model_target_chains_str
    )
    scores.update(interface_scores)
    scores["interface_AAs"] = interface_AA
    scores["interface_residues"] = interface_residues_str

    # 4. Calculate secondary structure
    logging.info("Calculating secondary structure...")
    ss_metrics = calc_ss_percentage(
        relaxed_pdb_path,
        dssp_path,
        binder_chain,
        target_chains_str=model_target_chains_str,
    )
    (
        scores["binder_helix%"],
        scores["binder_betasheet%"],
        scores["binder_loop%"],
        scores["interface_helix%"],
        scores["interface_betasheet%"],
        scores["interface_loop%"],
        scores["i_plddt"],
        scores["ss_plddt"],
    ) = ss_metrics

    # 5. Calculate RMSDs if a reference is provided
    if reference_pdb:
        logging.info("Calculating RMSDs...")
        scores["hotspot_rmsd"] = unaligned_rmsd(
            reference_pdb, relaxed_pdb_path, binder_chain, binder_chain
        )
        if reference_target_chains:
            if model_target_chains_str:
                scores["target_rmsd"] = target_pdb_rmsd(
                    trajectory_pdb=relaxed_pdb_path,
                    trajectory_target_chains_string=model_target_chains_str,
                    reference_pdb=reference_pdb,
                    reference_target_chains_string=reference_target_chains,
                )
            else:
                logging.warning(
                    f"No target chains found in {relaxed_pdb_path} to calculate RMSD against."
                )

    # 6. Analyze the sequence
    logging.info("Analyzing sequence...")
    sequence = get_sequence_from_pdb(relaxed_pdb_path, binder_chain)
    scores["sequence"] = sequence
    scores["extinction_coefficient"] = extinction_coefficient(sequence)
    scores["notes"] = validate_design_sequence(
        sequence, scores["relaxed_clashes"], omit_aas
    )

    # Clean up the temporary directory
    # shutil.rmtree(temp_dir) # Maybe don't delete immediately for debugging.

    return scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Score a binder-target PDB file.")
    parser.add_argument(
        "pdb_files",
        metavar="PDB_FILE",
        type=str,
        nargs="*",
        help="Path(s) to the PDB file(s) to score.",
    )
    parser.add_argument(
        "--pdbs-path",
        type=str,
        help="If specified, read all *.pdb files from this directory and ignore positional arguments.",
    )
    parser.add_argument(
        "--binder-chain",
        default="A",
        type=str,
        required=True,
        help="Chain ID of the binder.",
    )
    parser.add_argument(
        "--omit-aas",
        default="CX",
        type=str,
        help="Amino acids to omit from the sequence. (default: CX)",
    )
    parser.add_argument(
        "--reference-pdb",
        type=str,
        help="(Optional) Path to a reference PDB for RMSD calculations.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        help='(Optional) Comma-separated list of target chains for RMSD, e.g., "A,C"',
    )
    parser.add_argument(
        "--dalphaball-path",
        type=str,
        required=True,
        help="Explicitly specify the path to the dalphaball executable.",
    )
    parser.add_argument(
        "--dssp-path",
        type=str,
        required=True,
        help="Explicitly specify the path to the dssp executable.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="-",
        help="Path to save the output. Defaults to stdout (-).",
    )
    parser.add_argument(
        "--format",
        choices=["json", "tsv"],
        default="json",
        help="Output format (default: json).",
    )
    parser.add_argument(
        "--unpack-dicts",
        action="store_true",
        help="For TSV format, unpack dictionaries into separate columns.",
    )

    args = parser.parse_args()

    # Validate the dalphaball path
    dalphaball_path = args.dalphaball_path
    if not os.path.isfile(dalphaball_path) or not os.access(dalphaball_path, os.X_OK):
        raise FileNotFoundError(
            f"The dalphaball executable was not found or is not executable at '{dalphaball_path}'. "
            "Please check the path provided via --dalphaball-path."
        )

    # Validate the dssp path
    dssp_path = args.dssp_path
    if (
        not dssp_path
        or not os.path.isfile(dssp_path)
        or not os.access(dssp_path, os.X_OK)
    ):
        raise FileNotFoundError(
            f"The DSSP executable was not found or is not executable at '{dssp_path}'. "
            "Please provide a valid path using --dssp-path."
        )

    # Initialise PyRosetta
    with contextlib.redirect_stdout(sys.stderr):
        pr.init(
            f"-ignore_unrecognized_res -ignore_zero_occupancy -mute all -holes:dalphaball {dalphaball_path} -corrections::beta_nov16 true -relax:default_repeats 1"
        )

    # Set up logging to stderr
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    # Determine which PDB files to process
    if args.pdbs_path:
        pdb_files = [
            os.path.join(args.pdbs_path, f)
            for f in os.listdir(args.pdbs_path)
            if f.endswith(".pdb") and os.path.isfile(os.path.join(args.pdbs_path, f))
        ]
    else:
        pdb_files = args.pdb_files

    if not pdb_files:
        logging.error("No PDB files specified.")
        sys.exit(1)

    # Run scoring for all PDBs
    results = []
    for pdb_file in pdb_files:
        scores = score_pdb(
            pdb_file=pdb_file,
            binder_chain=args.binder_chain,
            dssp_path=dssp_path,
            omit_aas=args.omit_aas,
            reference_pdb=args.reference_pdb,
            reference_target_chains=args.target_chains,
        )
        results.append(OrderedDict([("filename", pdb_file), ("scores", scores)]))

    # Format output
    if args.format == "json":
        output_data = json.dumps({"results": results}, indent=4)
    elif args.format == "tsv":
        rows = []
        ordered_headers = []
        for result in results:
            row = OrderedDict([("filename", result["filename"])])
            scores = result["scores"]
            if args.unpack_dicts:
                for key, value in scores.items():
                    if isinstance(value, dict):
                        for sub_key, sub_value in value.items():
                            row[f"{key}_{sub_key}"] = sub_value
                    else:
                        row[key] = value
            else:
                for key, value in scores.items():
                    if isinstance(value, dict):
                        row[key] = json.dumps(value)
                    else:
                        row[key] = value

            for key in row.keys():
                if key not in ordered_headers:
                    ordered_headers.append(key)
            rows.append(row)

        lines = ["\t".join(ordered_headers)]
        for row in rows:
            values = [str(row.get(h, "")) for h in ordered_headers]
            lines.append("\t".join(values))
        output_data = "\n".join(lines) + "\n"

    # Write output
    if args.output == "-":
        print(output_data)
    else:
        with open(args.output, "w") as f:
            f.write(output_data)
        logging.info(f"Scores saved to {args.output}")
