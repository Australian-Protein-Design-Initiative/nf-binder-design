#!/usr/bin/env python

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
            # Get first and last residue numbers
            start = residues[0].get_id()[1]  # Residue ID is (hetflag, resseq, icode)
            end = residues[-1].get_id()[1]
            ranges.append((chain.id, start, end))

    # Sort by chain ID
    ranges.sort(key=lambda x: x[0])
    return ranges


def main():
    parser = argparse.ArgumentParser(
        description="Get the contigs string for RFdiffusion from a PDB file"
    )
    parser.add_argument("pdb_file", type=Path, help="Input PDB file")
    parser.add_argument("binder_chain", help="Chain ID of the binder")
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

    # Use provided target_contigs if specified, otherwise detect from PDB
    if args.target_contigs:
        target_contigs = args.target_contigs
    else:
        # Find target chain ranges (all non-binder chains)
        target_ranges = [r for r in ranges if r[0] != args.binder_chain]
        if not target_ranges:
            print(f"Error: No target chains found in {args.pdb_file}", file=sys.stderr)
            sys.exit(1)
        target_contigs = "/".join([f"{r[0]}{r[1]}-{r[2]}" for r in target_ranges])

    # Construct final contigs string
    print(f"[{binder_length}-{binder_length}/0 {target_contigs}]")


if __name__ == "__main__":
    main()
