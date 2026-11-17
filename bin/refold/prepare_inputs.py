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

--target-fasta replaces that read-back target with a sequence given on the command
line, which is how a design made against a truncated construct is refolded against
the full-length protein. It must hold exactly one record: the pulldown folds every
binder against every target, and one row per design downstream depends on there
being a single target. The design's own target is still read, and a warning says so
when it is not a subsequence of the one supplied - the usual sign that the chain or
the construct is not the one the run assumed.

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


def read_single_target(path: Path) -> str:
    """The one sequence in an external target FASTA.

    More than one record is rejected rather than silently taking the first: the
    pulldown folds every binder against every target, so a second record would double
    the fold jobs and give two rows per design where the score tables expect one.
    """
    records: List[tuple] = []
    name = None
    chunks: List[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))

    records = [(n, s) for n, s in records if s]
    if not records:
        raise SystemExit(f"--target-fasta {path}: no sequence records")
    if len(records) > 1:
        raise SystemExit(
            f"--target-fasta {path}: {len(records)} records "
            f"({', '.join(n for n, _ in records[:5])}). Refolding takes exactly one "
            "target; run it once per target to compare two of them."
        )
    return records[0][1].upper()


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
    parser.add_argument(
        "--target-fasta",
        help="FASTA holding the target to fold against, instead of the one read back out "
        "of the designs. Exactly one record; its sequence is used and --target-id still "
        "names it.",
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

    design_target = consensus_target(sequences)
    if args.target_fasta:
        target_seq = read_single_target(Path(args.target_fasta))
        if design_target not in target_seq:
            logging.warning(
                "--target-fasta sequence (%d residues) does not contain the target the "
                "designs were built against (%d residues) as a subsequence. Refolding "
                "will still run, but the binder was not designed against this target; "
                "check --refold_target_chain and that the FASTA is the same protein.",
                len(target_seq),
                len(design_target),
            )
        else:
            logging.info(
                "refolding against the supplied target (%d residues), which contains the "
                "designs' own target (%d residues)",
                len(target_seq),
                len(design_target),
            )
    else:
        target_seq = design_target

    write_fasta(Path(args.targets_out), [(args.target_id, target_seq)])
    write_fasta(
        Path(args.binders_out),
        [(design_id, sequences[design_id]["binder"]) for design_id in sorted(sequences)],
    )
    logging.info("wrote %d binder records against 1 target", len(sequences))
    return 0


if __name__ == "__main__":
    sys.exit(main())
