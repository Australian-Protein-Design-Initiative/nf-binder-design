#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Convert a structure to PDB, keeping chain ids.

Some scoring tools read PDB only - BindCraft's interface metrics go through
PyRosetta's pose_from_pdb - while most folding engines now write mmCIF. This is
the adapter between the two, and nothing more: no renumbering, no chain renaming,
no filtering.

Chain ids longer than one character cannot be written to PDB at all, so they are
rejected here with a clear message rather than silently truncated into a clash.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx


def drop_unmodelled_atoms(structure, source: str):
    """Drop atoms whose coordinates are NaN, keeping the rest of the structure.

    An engine that emits a row for an atom it did not place writes NaN
    coordinates for it - RFdiffusion3/ProteinMPNN do this for the terminal OXT.
    The PDB writer rejects NaN outright, which would fail the whole structure
    over an atom that was never modelled; the residue is still described by its
    remaining atoms, so only those atoms are dropped.
    """
    finite = ~np.isnan(structure.coord).any(axis=-1)
    if finite.all():
        return structure
    names = sorted(set(structure.atom_name[~finite]))
    print(
        f"WARNING: {source}: dropping {int((~finite).sum())} unmodelled atom(s) "
        f"with NaN coordinates ({', '.join(names)})",
        file=sys.stderr,
    )
    return structure[finite]


def read_structure(path: Path):
    suffixes = path.suffixes
    opener = gzip.open if suffixes and suffixes[-1] == ".gz" else open
    suffix = suffixes[-2] if (suffixes and suffixes[-1] == ".gz" and len(suffixes) >= 2) else path.suffix

    if suffix in (".pdb", ".pdb1", ".ent"):
        with opener(path, "rt") as handle:
            return pdb.PDBFile.read(handle).get_structure(model=1)
    if suffix in (".cif", ".mmcif"):
        with opener(path, "rt") as handle:
            return pdbx.get_structure(pdbx.CIFFile.read(handle), model=1)
    raise SystemExit(f"ERROR: unsupported structure format '{suffix}' ({path})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structure", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    structure = read_structure(args.structure)
    structure = drop_unmodelled_atoms(structure, str(args.structure))

    too_long = sorted({c for c in structure.chain_id if len(c) > 1})
    if too_long:
        print(
            f"ERROR: PDB holds one character per chain id, but {args.structure} uses "
            f"{', '.join(too_long)}. Score this structure from mmCIF instead.",
            file=sys.stderr,
        )
        return 1

    out = pdb.PDBFile()
    out.set_structure(structure)
    out.write(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
