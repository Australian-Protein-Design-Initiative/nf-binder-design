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
from pathlib import Path
from collections import defaultdict
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
            residues = [r for r in chain.get_residues() if "CA" in r]
            if not residues:
                continue

            # Initialize variables for tracking ranges
            current_ranges = []
            start = None
            prev_res: Optional[int] = None

            for residue in residues:
                res_id = residue.get_id()[1]

                # If this is the first residue, start a new range
                if start is None:
                    start = res_id
                # If there's a gap in residue numbering, end current range and start new one
                elif prev_res is not None and res_id != prev_res + 1:
                    current_ranges.append((chain.id, start, prev_res))
                    start = res_id

                prev_res = res_id

            # Add the last range
            if start is not None and prev_res is not None:
                current_ranges.append((chain.id, start, prev_res))

            ranges.extend(current_ranges)

    # Sort by chain ID and start residue
    ranges.sort(key=lambda x: (x[0], x[1]))
    return ranges


def format_chimerax_selection(ranges: List[Tuple[str, int, int]]) -> str:
    """Format ranges as a ChimeraX selection string.

    Args:
        ranges: List of (chain_id, start_res, end_res) tuples

    Returns:
        ChimeraX selection string like "@A:1-100 or @B:50-75"
    """
    if not ranges:
        return ""

    # Group ranges by chain ID
    chain_groups = defaultdict(list)
    for chain_id, start, end in ranges:
        chain_groups[chain_id].append((start, end))

    # Format each chain's ranges
    chain_selections = []
    for chain_id in sorted(chain_groups.keys()):
        ranges_list = sorted(chain_groups[chain_id])
        # Combine consecutive or overlapping ranges on the same chain
        range_parts = []
        for start, end in ranges_list:
            range_parts.append(f"{start}-{end}")
        chain_selection = f"@{chain_id}:" + ",".join(range_parts)
        chain_selections.append(chain_selection)

    return " or ".join(chain_selections)


def main():
    parser = argparse.ArgumentParser(
        description="""Get the contigs string for RFdiffusion from a PDB file.
If no binder chain is specified, we return the contigs for each of chains in the target.
If a binder chain is specified, we return the length of the binder chain and the contigs for the target chains.
Use --chimerax to output ChimeraX selection strings instead of RFdiffusion contigs.
"""
    )
    parser.add_argument("pdb_file", type=Path, help="Input PDB file")
    parser.add_argument(
        "binder_chain",
        nargs="?",
        default=None,
        help="Chain ID of the binder (optional)",
    )
    parser.add_argument(
        "--target_contigs",
        default=None,
        help="Target contigs string - if not provided, will auto-detect from PDB",
    )
    parser.add_argument(
        "--chimerax",
        action="store_true",
        help="Output ChimeraX selection strings instead of RFdiffusion contigs",
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

    # ChimeraX output format
    if args.chimerax:
        if args.binder_chain:
            # Find binder and target ranges
            binder_ranges = [r for r in ranges if r[0] == args.binder_chain]
            target_ranges = [r for r in ranges if r[0] != args.binder_chain]
            
            if not binder_ranges:
                print(
                    f"Error: Binder chain {args.binder_chain} not found in {args.pdb_file}",
                    file=sys.stderr,
                )
                sys.exit(1)
            
            if target_ranges:
                binder_sel = format_chimerax_selection(binder_ranges)
                target_sel = format_chimerax_selection(target_ranges)
                print("# ChimeraX contig selection")
                print(f"select {binder_sel}")
                print(f"select {target_sel}")
            else:
                binder_sel = format_chimerax_selection(binder_ranges)
                print("# ChimeraX contig selection")
                print(f"select {binder_sel}")
        else:
            # No binder specified, output all chains
            all_sel = format_chimerax_selection(ranges)
            print("# ChimeraX contig selection")
            print(f"select {all_sel}")
        return

    # RFdiffusion contigs output format (original behaviour)
    # Use provided target_contigs if specified, otherwise detect from PDB
    if args.target_contigs:
        target_contigs_str = args.target_contigs
        # If binder is also specified, format with binder length
        if args.binder_chain:
            try:
                binder_chain = next(r for r in ranges if r[0] == args.binder_chain)
                binder_length = binder_chain[2] - binder_chain[1] + 1
                print(f"[{binder_length}-{binder_length}/0 {target_contigs_str}]")
            except StopIteration:
                print(
                    f"Error: Binder chain {args.binder_chain} not found in {args.pdb_file}",
                    file=sys.stderr,
                )
                sys.exit(1)
        else:
            print(f"[{target_contigs_str}/0]")
        return

    if args.binder_chain:
        # Find target chain ranges (all non-binder chains)
        target_ranges = [r for r in ranges if r[0] != args.binder_chain]
    else:
        # No binder chain specified, use all chains as targets
        target_ranges = ranges

    if not target_ranges:
        print(f"Error: No target chains found in {args.pdb_file}", file=sys.stderr)
        sys.exit(1)

    # Group ranges by chain ID
    chain_groups = defaultdict(list)
    for chain_id, start, end in target_ranges:
        chain_groups[chain_id].append(f"{chain_id}{start}-{end}")

    # Format each chain group into a single slash-separated string, sorted by chain ID for consistency
    formatted_chain_contigs = [
        "/".join(chain_groups[chain_id]) for chain_id in sorted(chain_groups.keys())
    ]

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
        target_contigs_str = " ".join(formatted_chain_contigs)
        print(f"[{binder_length}-{binder_length}/0 {target_contigs_str}]")
    else:
        # No binder chain specified, output each chain group with a /0 suffix
        target_contigs_parts = [f"{contig}/0" for contig in formatted_chain_contigs]
        print(f"[{' '.join(target_contigs_parts)}]")


if __name__ == "__main__":
    main()
