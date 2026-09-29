#!/usr/bin/env python3
"""Thread a ProteinMPNN-designed sequence onto a backbone and relax it.
Runs inside the same container as the baseline dl_binder_design tool, reusing
its PyRosetta FastRelax config, so it doesn't need its own MPNN call.
"""
import os

from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta import core, protocols

BACKBONE_PDB = os.environ["THREAD_BACKBONE_PDB"]
FASTA_PATH = os.environ["THREAD_FASTA"]
OUTPUT_PDB = os.environ["THREAD_OUTPUT_PDB"]
RELAX_XML = os.environ["THREAD_RELAX_XML"]

ALPHA_1 = list("ARNDCQEGHILKMFPSTWYV")
ALPHA_3 = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]
AA_1_3 = dict(zip(ALPHA_1, ALPHA_3))


def read_designed_sequence(fasta_path):
    records = []
    header, seq = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line, []
            elif line:
                seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    for header, seq in records:
        if not header.startswith(">input"):
            return seq
    raise ValueError(f"No designed sequence found in {fasta_path}")


init("-beta_nov16 -in:file:silent_struct_type binary -mute all"
     " -use_terminal_residues true -mute basic.io.database core.scoring")

sequence = read_designed_sequence(FASTA_PATH)
pose = pose_from_pdb(BACKBONE_PDB)

rsd_set = pose.residue_type_set_for_pose(core.chemical.FULL_ATOM_t)
for resi, aa in enumerate(sequence, start=1):
    new_res = core.conformation.ResidueFactory.create_residue(
        rsd_set.name_map(AA_1_3[aa])
    )
    pose.replace_residue(resi, new_res, True)

objs = protocols.rosetta_scripts.XmlObjects.create_from_file(RELAX_XML)
fast_relax = objs.get_mover("FastRelax")
fast_relax.apply(pose)

pose.dump_pdb(OUTPUT_PDB)
