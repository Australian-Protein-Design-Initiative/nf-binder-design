#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "biopython",
# ]
# ///

"""
Turn a directory of binder design complexes into the targets.fasta / binders.fasta
pair that the fold_pulldown machinery takes.

Every design holds the same target, so the target FASTA gets one record and the
binder FASTA one record per design, named after the design file. The pulldown then
co-folds each design's binder sequence against the target with whichever oracles
were selected, and the pairs come back keyed by the design id.

Chain extraction is delegated to bin/pdb_to_fasta.py, so the chain handling (CIF and
PDB, gzip, non-standard residues) stays in one place.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List

_PDB_TO_FASTA = Path(__file__).resolve().parents[1] / "pdb_to_fasta.py"


_pdb_to_fasta = None


def _load_pdb_to_fasta():
    global _pdb_to_fasta
    if _pdb_to_fasta is None:
        spec = importlib.util.spec_from_file_location("pdb_to_fasta", _PDB_TO_FASTA)
        _pdb_to_fasta = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_pdb_to_fasta)
    return _pdb_to_fasta


def _design_id(path: Path) -> str:
    name = path.name
    for suffix in (".cif.gz", ".pdb.gz", ".cif", ".pdb"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def collect_sequences(
    structures: List[Path], binder_chain: str, target_chain: str
) -> Dict[str, Dict[str, str]]:
    """Return {design_id: {"binder": seq, "target": seq}} for every readable design."""
    mod = _load_pdb_to_fasta()
    out: Dict[str, Dict[str, str]] = {}
    for path in structures:
        chains = mod.get_chain_sequences(str(path), [binder_chain, target_chain])
        design_id = _design_id(path)
        binder = chains.get(binder_chain)
        target = chains.get(target_chain)
        if not binder:
            logging.warning("%s: no chain %s, skipping", path.name, binder_chain)
            continue
        if not target:
            logging.warning("%s: no chain %s, skipping", path.name, target_chain)
            continue
        out[design_id] = {"binder": binder, "target": target}
    return out


def consensus_target(sequences: Dict[str, Dict[str, str]]) -> str:
    """The target sequence shared by the designs.

    Preserved contig segments are copied from the input, so every design should carry
    an identical target. A disagreement means the contig did not pin the target the
    way the run assumed, which would make the pulldown fold a target that only some
    designs were built against - so warn loudly and go with the majority.
    """
    counts = Counter(v["target"] for v in sequences.values())
    if len(counts) > 1:
        logging.warning(
            "designs disagree on the target sequence (%d distinct, lengths %s); "
            "using the most common one",
            len(counts),
            sorted({len(s) for s in counts}),
        )
    return counts.most_common(1)[0][0]


def write_fasta(path: Path, records: List[tuple]) -> None:
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "structures",
        nargs="+",
        help="Design complexes (.pdb/.cif, optionally gzipped), or directories of them",
    )
    parser.add_argument("--binder-chain", required=True, help="Binder chain id in the design")
    parser.add_argument("--target-chain", required=True, help="Target chain id in the design")
    parser.add_argument(
        "--target-id", default="target", help="Name for the target record [default: target]"
    )
    parser.add_argument("--targets-out", default="targets.fasta")
    parser.add_argument("--binders-out", default="binders.fasta")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    mod = _load_pdb_to_fasta()
    structures = [Path(p) for p in mod.expand_structure_inputs(args.structures)]
    if not structures:
        logging.error("no structures found in %s", args.structures)
        return 1

    sequences = collect_sequences(structures, args.binder_chain, args.target_chain)
    if not sequences:
        logging.error(
            "no design had both chain %s and chain %s", args.binder_chain, args.target_chain
        )
        return 1

    write_fasta(Path(args.targets_out), [(args.target_id, consensus_target(sequences))])
    write_fasta(
        Path(args.binders_out),
        [(design_id, sequences[design_id]["binder"]) for design_id in sorted(sequences)],
    )
    logging.info("wrote %d binder records against 1 target", len(sequences))
    return 0


if __name__ == "__main__":
    sys.exit(main())
