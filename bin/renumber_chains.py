#!/usr/bin/env python3
import argparse
from Bio.PDB import PDBParser, PDBIO
import string

# See: https://github.com/RosettaCommons/RFdiffusion/issues/214#issuecomment-2292225175
#      https://github.com/RosettaCommons/RFdiffusion/issues/214#issuecomment-2402728749
#
# Residues and chains for partial diffusion need to be renumbered to be contiguous,
# and the diffusable chains (binder) need to be the first chains (chain A), before the
# non-diffusable chains.

# This especially applies to structures output by BindCraft, which have the binder
# as chain B.


def renumber_residues(structure, ignore_chains=False, start_resnum=1):
    if ignore_chains:
        # Get all residues across all chains and sort by original residue number
        all_residues = []
        for model in structure:
            for chain in model:
                all_residues.extend(chain.get_residues())
        all_residues = sorted(all_residues, key=lambda r: r.id[1])
        
        # Renumber sequentially starting from start_resnum
        for new_number, residue in enumerate(all_residues, start=start_resnum):
            new_id = (' ', new_number, ' ')
            residue.id = new_id
    else:
        # Original chain-specific renumbering
        for model in structure:
            for chain in model:
                residues = sorted(chain.get_residues(), key=lambda r: r.id[1])
                for new_number, residue in enumerate(residues, start=start_resnum):
                    new_id = (' ', new_number, ' ')
                    residue.id = new_id

def reorder_chains(structure, binder_chains):
    if not binder_chains:
        return
    
    # Convert single chain to list
    if isinstance(binder_chains, str):
        binder_chains = [binder_chains]
    
    # Get all chains in the structure
    all_chains = []
    for model in structure:
        all_chains.extend(model.get_chains())
    
    # Separate binder chains and other chains
    binder_chains_list = [chain for chain in all_chains if chain.id in binder_chains]
    other_chains = [chain for chain in all_chains if chain.id not in binder_chains]
    
    # Generate new chain IDs using all uppercase letters
    new_chain_ids = list(string.ascii_uppercase)
    
    # First assign temporary unique IDs to avoid conflicts
    temp_ids = [f"TEMP{i}" for i in range(len(all_chains))]
    for chain, temp_id in zip(all_chains, temp_ids):
        chain.id = temp_id
    
    # Now assign the final chain IDs
    current_index = 0
    for chain in binder_chains_list:
        chain.id = new_chain_ids[current_index]
        current_index += 1
    for chain in other_chains:
        chain.id = new_chain_ids[current_index]
        current_index += 1

def main():
    parser = argparse.ArgumentParser(description='Renumber PDB residues sequentially within each chain')
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('--ignore-chains', action='store_true',
                       help='Number residues sequentially across all chains')
    parser.add_argument('--binder-chains', nargs='+', default=None,
                       help='Chain(s) to move to the front (e.g., --binder-chains B will rename chain B to A)')
    parser.add_argument('--start-resnum', type=int, default=1,
                       help='Residue number to start counting from (default: 1)')
    args = parser.parse_args()

    # Parse the input PDB file
    parser = PDBParser()
    structure = parser.get_structure('input', args.pdb_file)

    # Reorder chains if binder_chains is specified
    if args.binder_chains:
        reorder_chains(structure, args.binder_chains)

    # Renumber residues
    renumber_residues(structure, args.ignore_chains, args.start_resnum)

    # Output the modified structure with chains in order
    io = PDBIO()
    io.set_structure(structure)
    
    # Sort chains alphabetically before saving
    for model in structure:
        model.child_list.sort(key=lambda x: x.id)
    
    io.save(sys.stdout)

if __name__ == '__main__':
    import sys
    main()
