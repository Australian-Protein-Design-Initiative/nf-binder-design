#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///
"""
Split a monomer-mode AF2 complex prediction back into chains at the chain breaks.

The monomer models have no chain concept, so a complex folded via
af2_monomer_features_from_msas.py comes back as a single chain whose residue numbering
carries the +offset jumps we put into `residue_index`. Everything downstream
(FOLD_PARSE_CONFIDENCE, ipsae, rmsd4all, pose_rmsd) selects by chain, so that numbering
has to become real chains before any of it will work.

Chains are cut on the layout emitted by the features builder, never by looking for
jumps in the numbering: a fold can legitimately contain a large gap, and guessing would
silently mis-split it. Each chain is renumbered 1..L so the output matches what every
other engine emits.

Uses AF2's own alphafold.common.protein for IO, so the atom/chain conventions and the
mmCIF writer are identical to the ones AF2 used to write the input (gemmi is not in
the container). Must run inside the AF2 container.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, "/app/alphafold")

import numpy as np  # noqa: E402
from alphafold.common import protein  # noqa: E402


def split_protein(prot: protein.Protein, layout: list) -> protein.Protein:
    total = sum(c["length"] for c in layout)
    n = prot.aatype.shape[0]
    if n != total:
        shape = ", ".join("{}={}".format(c["chain_id"], c["length"]) for c in layout)
        raise SystemExit(f"{n} residues but layout expects {total} ({shape})")

    chain_index = np.concatenate(
        [np.full(c["length"], i, dtype=np.int32) for i, c in enumerate(layout)]
    )
    residue_index = np.concatenate(
        [np.arange(1, c["length"] + 1, dtype=np.int32) for c in layout]
    )
    return protein.Protein(
        atom_positions=prot.atom_positions,
        aatype=prot.aatype,
        atom_mask=prot.atom_mask,
        residue_index=residue_index,
        chain_index=chain_index,
        b_factors=prot.b_factors,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--layout", required=True,
                    help="chain layout JSON from af2_monomer_features_from_msas.py")
    ap.add_argument("--outdir", help="output directory (default: overwrite in place)")
    ap.add_argument("--also-mmcif", action="store_true",
                    help="write a sibling .cif for every .pdb processed")
    ap.add_argument("structures", nargs="+", help="PDB files written by AF2")
    args = ap.parse_args()

    spec = json.loads(Path(args.layout).read_text())
    layout = spec["chains"]
    ids = [c["chain_id"] for c in layout]
    # protein.to_pdb assigns chain letters from chain_index in PDB_CHAIN_IDS order,
    # so the layout must already be in that order for the letters to come out right.
    if ids != list(protein.PDB_CHAIN_IDS)[: len(ids)]:
        raise SystemExit(f"layout chain ids {ids} are not in PDB_CHAIN_IDS order")

    n = 0
    for s in args.structures:
        src = Path(s)
        if not src.is_file() or src.suffix.lower() != ".pdb":
            continue
        prot = split_protein(protein.from_pdb_string(src.read_text()), layout)
        dst_dir = Path(args.outdir) if args.outdir else src.parent
        dst_dir.mkdir(parents=True, exist_ok=True)
        (dst_dir / src.name).write_text(protein.to_pdb(prot))
        if args.also_mmcif:
            (dst_dir / (src.stem + ".cif")).write_text(
                protein.to_mmcif(prot, src.stem, "Monomer")
            )
        print(f"{src.name}: {prot.aatype.shape[0]} residues -> chains {'/'.join(ids)}",
              file=sys.stderr)
        n += 1
    if n == 0:
        raise SystemExit("no .pdb inputs matched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
