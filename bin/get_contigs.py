#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "biopython",
# ]
# ///

from typing import List, Tuple, Optional
import argparse
import sys
import json
from pathlib import Path
from Bio import PDB


def get_chain_ranges(structure: PDB.Structure.Structure) -> List[Tuple[str, int, int]]:
    """Get the residue ranges for each chain.

    Args:
        structure: BioPython Structure object

    Returns:
        List of (chain_id, start_res, end_res) tuples
    """
    ranges = []
    for model in structure:
        for chain in model:
            residues = [r for r in chain.get_residues()]
            if not residues:
                continue
            
            # Initialize variables for tracking ranges
            current_ranges = []
            start = None
            prev_res = None
            
            for residue in residues:
                res_id = residue.get_id()[1]
                
                # If this is the first residue, start a new range
                if start is None:
                    start = res_id
                # If there's a gap in residue numbering, end current range and start new one
                elif res_id != prev_res + 1:
                    current_ranges.append((chain.id, start, prev_res))
                    start = res_id
                
                prev_res = res_id
            
            # Add the last range
            if start is not None:
                current_ranges.append((chain.id, start, prev_res))
            
            ranges.extend(current_ranges)

    # Sort by chain ID and start residue
    ranges.sort(key=lambda x: (x[0], x[1]))
    return ranges


def main():
    parser = argparse.ArgumentParser(
        description="""Get the contigs string for RFdiffusion from a PDB file.
If no binder chain is specified, we return the contigs for each of chains in the target.
If a binder chain is specified, we return the length of the binder chain and the contigs for the target chains.
"""
    )
    parser.add_argument("pdb_file", type=Path, help="Input PDB file")
    parser.add_argument(
        "binder_chain", 
        nargs="?", 
        default=None, 
        help="Chain ID of the binder (optional)"
    )
    parser.add_argument(
        "--target_contigs",
        default=None,
        help="Target contigs string - if not provided, will auto-detect from PDB",
    )
    args = parser.parse_args()

    if not args.pdb_file.exists():
        print(f"Error: PDB file {args.pdb_file} does not exist", file=sys.stderr)
        sys.exit(1)

    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("pdb", args.pdb_file)
    except Exception as e:
        print(f"Error parsing PDB file: {e}", file=sys.stderr)
        sys.exit(1)

    ranges = get_chain_ranges(structure)

    # Use provided target_contigs if specified, otherwise detect from PDB
    if args.target_contigs:
        target_contigs = args.target_contigs
    else:
        if args.binder_chain:
            # Find target chain ranges (all non-binder chains)
            target_ranges = [r for r in ranges if r[0] != args.binder_chain]
        else:
            # No binder chain specified, use all chains as targets
            target_ranges = ranges
            
        if not target_ranges:
            print(f"Error: No target chains found in {args.pdb_file}", file=sys.stderr)
            sys.exit(1)
        target_contigs = "/".join([f"{r[0]}{r[1]}-{r[2]}" for r in target_ranges])

    # Construct final contigs string
    if args.binder_chain:
        # Find binder chain length
        try:
            binder_chain = next(r for r in ranges if r[0] == args.binder_chain)
        except StopIteration:
            print(
                f"Error: Binder chain {args.binder_chain} not found in {args.pdb_file}",
                file=sys.stderr,
            )
            sys.exit(1)
        binder_length = binder_chain[2] - binder_chain[1] + 1
        print(f"[{binder_length}-{binder_length}/0 {target_contigs}]")
    else:
        # No binder chain specified, just output the target contigs
        print(f"[{target_contigs}/0]")


if __name__ == "__main__":
    main()
