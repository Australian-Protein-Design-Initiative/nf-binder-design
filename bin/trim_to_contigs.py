#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "biopython",
# ]
# ///

# Given a PDB and an RFDiffusion style contigs string, trim each chain to the specified contigs
# Return the PDB to stdout (--output -), or --output {file} if specified

import argparse
import sys
import re
from pathlib import Path
from collections import defaultdict
from Bio import PDB
from io import StringIO


class ContigSelect(PDB.Select):
    """
    Selects residues based on a list of contigs.
    """

    def __init__(self, ranges_to_keep: defaultdict[str, list[tuple[int, int]]]):
        self.ranges_to_keep = ranges_to_keep

    def accept_chain(self, chain: PDB.Chain.Chain) -> bool:
        return chain.get_id() in self.ranges_to_keep

    def accept_residue(self, residue: PDB.Residue.Residue) -> bool:
        chain_id = residue.get_parent().get_id()
        res_id = residue.get_id()[1]

        if chain_id in self.ranges_to_keep:
            for start, end in self.ranges_to_keep[chain_id]:
                if start <= res_id <= end:
                    return True
        return False


def parse_contigs(contig_string: str) -> list[tuple[str, int, int]]:
    """
    Parses an RFdiffusion-style contig string into a list of residue ranges.

    Supports:
    - v1 style: '[A18-132/0 65-120]' (space-separated, brackets, /0 attached)
    - v3 style: 'A18-132,/0,65-120' (comma-separated, /0 as separate element)

    Binder-length-only segments (e.g. 65-120 with no chain ID) are ignored.
    Returns e.g. [('A', 18, 132)] for trimming the target chain.
    """
    contig_string = contig_string.strip()
    if contig_string.startswith("[") and contig_string.endswith("]"):
        contig_string = contig_string[1:-1]

    # Normalise to tokens: split on comma and space so v3 "A18-132,/0,65-120" and v1 "A18-132/0 65-120" both work
    parts = re.split(r"[,\s]+", contig_string)

    residue_ranges = []
    regex = re.compile(r"^([A-Za-z]+)(-?\d+)-(-?\d+)$")
    # Also match segment that has / in it (v1): "A18-132/0" -> take "A18-132"
    for part in parts:
        part = part.strip()
        if not part or part in ("/0", "0"):
            continue
        # v1 may have "A18-132/0" as one token after space split; split on / and take first
        if "/" in part:
            part = part.split("/")[0]
        if not part:
            continue
        match = regex.match(part)
        if match:
            chain_id = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
            residue_ranges.append((chain_id, start, end))
        # num-num only (e.g. 65-120) is binder length, no chain in input structure -> skip

    return residue_ranges


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Trim a PDB file based on an RFdiffusion contig string."
    )
    parser.add_argument("pdb_file", type=Path, help="Input PDB file.")
    parser.add_argument(
        "contigs", type=str, help="RFdiffusion-style contig string."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        help="Output PDB file path. Default is stdout ('-').",
    )

    args = parser.parse_args()

    if not args.pdb_file.exists():
        print(f"Error: PDB file {args.pdb_file} does not exist", file=sys.stderr)
        sys.exit(1)

    pdb_parser = PDB.PDBParser(QUIET=True)
    try:
        structure = pdb_parser.get_structure("pdb_structure", args.pdb_file)
    except Exception as e:
        print(f"Error parsing PDB file: {e}", file=sys.stderr)
        sys.exit(1)

    parsed_ranges = parse_contigs(args.contigs)

    ranges_to_keep = defaultdict(list)
    for chain_id, start, end in parsed_ranges:
        ranges_to_keep[chain_id].append((start, end))

    io = PDB.PDBIO()
    io.set_structure(structure)
    selector = ContigSelect(ranges_to_keep)

    if args.output == "-":
        with StringIO() as output_handle:
            io.save(output_handle, selector)
            print(output_handle.getvalue(), end="")
    else:
        with open(args.output, "w") as output_handle:
            io.save(output_handle, selector)


if __name__ == "__main__":
    main()