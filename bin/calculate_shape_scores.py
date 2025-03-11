#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "mdanalysis",
#     "numpy",
#     "biopython",
# ]
# ///

import sys
import argparse
import logging
from typing import List, Optional, Dict, Tuple
import MDAnalysis as mda
from MDAnalysis.lib.distances import self_distance_array
import numpy as np
from Bio import PDB
from Bio.PDB.Polypeptide import protein_letters_3to1, is_aa


def calculate_rg(all_atoms: mda.AtomGroup) -> float:
    return all_atoms.radius_of_gyration()

def calculate_dmax(all_atoms: mda.AtomGroup) -> float:
    positions = all_atoms.positions
    return np.max(self_distance_array(positions))

def calculate_asphericity(all_atoms: mda.AtomGroup) -> float:
    positions = all_atoms.positions
    mass = all_atoms.masses[:, np.newaxis]
    total_mass = np.sum(mass)
    center_of_mass = np.sum(mass * positions, axis=0) / total_mass
    relative_positions = positions - center_of_mass
    rg_tensor = np.dot((mass * relative_positions).T, relative_positions) / total_mass
    
    eigenvalues = np.sort(np.linalg.eigvalsh(rg_tensor))
    return ((eigenvalues[2] - eigenvalues[1]) ** 2 + 
            (eigenvalues[1] - eigenvalues[0]) ** 2 + 
            (eigenvalues[0] - eigenvalues[2]) ** 2) / (2 * (sum(eigenvalues)) ** 2)

def calculate_stokes_radius(all_atoms: mda.AtomGroup) -> float:
    return (5 / 3) * all_atoms.radius_of_gyration()

def find_chains(universe: mda.Universe) -> List[str]:
    """Find all chain IDs in the structure."""
    # Try to get chainIDs first
    chains = set(universe.atoms.chainIDs)
    # Filter out empty chain IDs
    chains = [c for c in chains if c and c.strip()]
    
    # If no chainIDs found, try segids
    if not chains:
        chains = set(universe.atoms.segids)
        chains = [c for c in chains if c and c.strip()]
    
    return sorted(chains)

def count_residues(atoms: mda.AtomGroup) -> int:
    """Count the number of residues in an atom group."""
    return len(atoms.residues)

def get_chain_sequences(pdb_file: str, chain_ids: Optional[List[str]] = None) -> Dict[str, str]:
    """Extract one-letter amino acid sequences for specified chains in a PDB file."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    sequences = {}
    
    # Process each chain in the structure
    for chain in structure.get_chains():
        chain_id = chain.id.strip()
        if not chain_id:
            continue
            
        # Skip if specific chains requested and this one isn't in the list
        if chain_ids and chain_id not in chain_ids:
            continue
            
        # Extract amino acid sequence
        sequence = ""
        for residue in chain:
            if is_aa(residue):
                try:
                    # Use the protein_letters_3to1 dictionary instead of three_to_one function
                    res_name = residue.get_resname().upper()
                    one_letter = protein_letters_3to1.get(res_name, "X")
                    sequence += one_letter
                except (KeyError, AttributeError):
                    # Skip non-standard amino acids
                    sequence += "X"
        
        if sequence:
            sequences[chain_id] = sequence
            
    return sequences

def process_structures(pdb_files: List[str], chain: Optional[str] = None) -> List[tuple]:
    results = []
    for pdb_file in pdb_files:
        logging.info(f"Processing {pdb_file}")
        
        u = mda.Universe(pdb_file)
        
        # Get all chains in the structure
        all_chains = find_chains(u)
        chain_str = ",".join(all_chains)
        logging.info(f"Found chains: {chain_str}")
        
        # Extract chain IDs for sequence calculation
        chain_ids = None
        if chain:
            chain_ids = [c.strip() for c in chain.split(",")]
        
        # Get sequences for the chains
        sequences = get_chain_sequences(pdb_file, chain_ids)
        if chain_ids:
            sequence_str = ",".join([sequences.get(c, "") for c in chain_ids])
        else:
            sequence_str = ",".join([sequences.get(c, "") for c in all_chains])
        
        # Select atoms based on chain if specified, otherwise select all atoms
        if chain and chain_ids:  # Ensure chain_ids is not None before iterating
            # Handle comma-separated list of chains
            selection_parts = []
            for chain_id in chain_ids:
                selection_parts.append(f"(segid {chain_id} or chainID {chain_id})")
            
            selection = " or ".join(selection_parts)
            logging.info(f"Selection query: {selection}")
            
            atoms = u.select_atoms(selection)
            if len(atoms) == 0:
                logging.warning(f"No atoms found for chain(s) '{chain}' in {pdb_file}, skipping")
                continue
            logging.info(f"Selected {len(atoms)} atoms from chain(s) '{chain}'")
            chain_str = chain  # Use the specified chain(s) for output
            res_count = count_residues(atoms)
        else:
            atoms = u.select_atoms("all")
            # Count residues in all chains
            res_count = count_residues(atoms)
        
        logging.info(f"Residue count: {res_count}")
        
        # Calculate properties
        Rg = calculate_rg(atoms)
        Dmax = calculate_dmax(atoms)
        asphericity = calculate_asphericity(atoms)
        Rh = calculate_stokes_radius(atoms)
        
        # Log results
        chain_info = f" (chain(s) {chain})" if chain else ""
        logging.info(f"Radius of Gyration (Rg, Å){chain_info}: {Rg:.2f} Å")
        logging.info(f"Maximum Distance (Dmax, Å){chain_info}: {Dmax:.2f} Å")
        logging.info(f"Asphericity{chain_info}: {asphericity:.4f}")
        logging.info(f"Approximate Stokes Radius (Rh, Å){chain_info}: {Rh:.2f} Å\n")
        
        results.append((pdb_file, Rg, Dmax, asphericity, Rh, chain_str, res_count, sequence_str))
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Calculate structural properties of PDB files"
    )
    parser.add_argument("files", nargs="+", help="PDB files to process")
    parser.add_argument(
        "-o", "--output", help="Output file (default: stdout)", default="-"
    )
    parser.add_argument(
        "--chain", help="Chain ID(s) to analyze (comma-separated for multiple chains)", default=None
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", 
        help="Enable verbose output with detailed information about calculations"
    )
    args = parser.parse_args()

    # Set logging level based on verbose flag
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=log_level, format="%(message)s", stream=sys.stderr)

    results = process_structures(args.files, args.chain)

    # Prepare output
    header = "filename\trg\tdmax\tasphericity\tapprox_rh\tchain\tlength\tsequence\n"
    output_lines = [header]
    output_lines.extend(
        f"{fname}\t{rg:.2f}\t{dmax:.2f}\t{asph:.4f}\t{rh:.2f}\t{chain}\t{length}\t{seq}\n"
        for fname, rg, dmax, asph, rh, chain, length, seq in results
    )

    if args.output == "-":
        sys.stdout.writelines(output_lines)
    else:
        with open(args.output, "w") as f:
            f.writelines(output_lines)


if __name__ == "__main__":
    main()
