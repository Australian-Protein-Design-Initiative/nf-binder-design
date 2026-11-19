#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "biopython",
# ]
# ///

"""
Turn a set of target-binder complexes into the binders.fasta (and optionally
targets.fasta) that the fold_pulldown machinery takes, plus a chains.tsv recording
which chain of each complex is the binder and which is the target.

The complexes serve two purposes downstream: they are where the binder sequences come
from, and they are the reference each predicted complex is compared against once the
pulldown has folded it. chains.tsv carries the per-complex chain ids over to that
comparison, which is why the chain resolution happens here, once, rather than being
re-derived later.

Chain extraction is delegated to bin/pdb_to_fasta.py, so the chain handling (CIF and
PDB, gzip, non-standard residues) stays in one place. Its --chains output is not used
directly because it names every record after the file alone, which collides as soon as
one complex contributes more than one binder chain.

Target chain resolution (--target-chains):
  <id>   that chain, in every complex
  auto   per complex, the chain whose sequence best matches a --targets record
  unset  per complex, the single chain that is not a binder chain

--targets-out writes the target sequences read back out of the complexes, for a run
that gave no --targets of its own. They are deduplicated by sequence: binders that
share a target then fold against one target record rather than one each.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

_PDB_TO_FASTA = Path(__file__).resolve().parents[1] / "pdb_to_fasta.py"

_pdb_to_fasta = None

STRUCTURE_SUFFIXES = (".cif.gz", ".pdb.gz", ".cif", ".pdb")


def _load_pdb_to_fasta():
    global _pdb_to_fasta
    if _pdb_to_fasta is None:
        spec = importlib.util.spec_from_file_location("pdb_to_fasta", _PDB_TO_FASTA)
        _pdb_to_fasta = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_pdb_to_fasta)
    return _pdb_to_fasta


def complex_id(path: Path) -> str:
    """The stem a complex's records are named after, minus any structure suffix."""
    name = path.name
    for suffix in STRUCTURE_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    """(id, sequence) per record, without pulling in Biopython for a two-line parse."""
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:].split()[0] if len(line) > 1 else ""
            chunks = []
        elif line and header is not None:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def identity(a: str, b: str) -> float:
    """Ungapped identity over the shorter sequence.

    Deliberately not an alignment: this only has to pick the target chain out of the
    handful of chains in one complex, where the right answer is near-identical to a
    --targets record and the wrong answer is the binder, which is a different molecule
    of a very different length. An alignment would add a dependency to separate cases
    that are already far apart.
    """
    if not a or not b:
        return 0.0
    short, long_ = (a, b) if len(a) <= len(b) else (b, a)
    best = 0
    # Slide the shorter over the longer: a truncated construct is a substring of the
    # full-length target, not aligned at residue 1.
    for offset in range(len(long_) - len(short) + 1):
        matches = sum(1 for i, c in enumerate(short) if long_[offset + i] == c)
        best = max(best, matches)
    return best / len(short)


def resolve_chains(
    chains: Dict[str, str],
    binder_chains: Sequence[str],
    target_chains: Optional[Sequence[str]],
    target_seqs: Sequence[str],
    cid: str,
) -> Tuple[List[str], str]:
    """Return (binder chain ids present, target chain id) for one complex."""
    present_binders = [c for c in binder_chains if c in chains]
    missing = [c for c in binder_chains if c not in chains]
    if missing:
        logging.warning(
            "%s: no chain %s (has %s)", cid, ",".join(missing), ",".join(sorted(chains))
        )
    if not present_binders:
        raise ValueError(
            f"{cid}: none of the binder chains {','.join(binder_chains)} are present "
            f"(has {','.join(sorted(chains))})"
        )

    if target_chains and target_chains != ["auto"]:
        target = target_chains[0]
        if target not in chains:
            raise ValueError(
                f"{cid}: no target chain {target} (has {','.join(sorted(chains))})"
            )
        return present_binders, target

    candidates = [c for c in sorted(chains) if c not in present_binders]
    if not candidates:
        raise ValueError(f"{cid}: no chain left to be the target")

    if target_chains == ["auto"]:
        if not target_seqs:
            raise ValueError(
                f"{cid}: --target-chains auto needs --targets to match against"
            )
        scored = [
            (max(identity(chains[c], t) for t in target_seqs), c) for c in candidates
        ]
        scored.sort(reverse=True)
        best, chain = scored[0]
        logging.info("%s: --target-chains auto chose %s (identity %.2f)", cid, chain, best)
        if best < 0.5:
            logging.warning(
                "%s: best target match is chain %s at identity %.2f, which is low "
                "enough that the superposition reference may be the wrong chain",
                cid,
                chain,
                best,
            )
        return present_binders, chain

    if len(candidates) > 1:
        raise ValueError(
            f"{cid}: {len(candidates)} chains ({','.join(candidates)}) are not binder "
            "chains, so the target is ambiguous. Give --target_chains, or "
            "--target_chains auto to pick by sequence."
        )
    return present_binders, candidates[0]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
    )
    parser.add_argument("structures", nargs="+", help="Complex PDB/CIF files or directories")
    parser.add_argument(
        "--binder-chains",
        required=True,
        help="Comma-separated binder chain id(s) in the complexes",
    )
    parser.add_argument(
        "--target-chains",
        default=None,
        help="Target chain id, or 'auto' to pick by sequence. Default: the one "
        "non-binder chain.",
    )
    parser.add_argument(
        "--targets",
        default=None,
        help="FASTA of targets, used by --target-chains auto to identify the chain",
    )
    parser.add_argument("--binders-out", default="binders.fasta")
    parser.add_argument("--targets-out", default=None, help="Write targets read from the complexes")
    parser.add_argument("--chains-out", default="chains.tsv")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s - %(message)s",
        stream=sys.stderr,
    )

    mod = _load_pdb_to_fasta()
    structures = [Path(p) for p in mod.expand_structure_inputs(args.structures)]
    if not structures:
        logging.error("no structure files found in %s", " ".join(args.structures))
        return 1

    binder_chains = [c.strip() for c in args.binder_chains.split(",") if c.strip()]
    if not binder_chains:
        logging.error("--binder-chains is empty")
        return 1
    target_chains = (
        [c.strip() for c in args.target_chains.split(",") if c.strip()]
        if args.target_chains
        else None
    )
    if target_chains and len(target_chains) > 1:
        logging.error(
            "--target_chains takes one chain id or 'auto', got '%s'. A pulldown folds "
            "two-chain pairs (target as chain A, binder as chain B), so a multi-chain "
            "target cannot be predicted and there is nothing to superimpose it onto.",
            args.target_chains,
        )
        return 1

    target_seqs = [s for _, s in read_fasta(Path(args.targets))] if args.targets else []

    binder_records: List[Tuple[str, str]] = []
    chain_rows: List[Tuple[str, str, str, str]] = []
    complex_targets: List[Tuple[str, str]] = []

    for path in structures:
        cid = complex_id(path)
        chains = mod.get_chain_sequences(str(path))
        if not chains:
            logging.warning("%s: no sequences, skipping", path.name)
            continue
        try:
            present_binders, target_chain = resolve_chains(
                chains, binder_chains, target_chains, target_seqs, cid
            )
        except ValueError as exc:
            logging.error("%s", exc)
            return 1

        for chain in present_binders:
            # The id is the join key for every downstream table, so it has to be unique
            # and stable. Only suffix when more than one binder chain was asked for, so
            # the common single-chain case keeps the complex's own name.
            rec_id = cid if len(binder_chains) == 1 else f"{cid}_{chain}"
            binder_records.append((rec_id, chains[chain]))
            chain_rows.append((rec_id, path.name, chain, target_chain))

        complex_targets.append((cid, chains[target_chain]))

    if not binder_records:
        logging.error("no binder sequences extracted")
        return 1

    seen = set()
    duplicates = {r for r, _ in binder_records if r in seen or seen.add(r)}
    if duplicates:
        logging.error(
            "duplicate binder ids %s: complexes must have distinct filenames",
            ",".join(sorted(duplicates)),
        )
        return 1

    if args.targets:
        target_ids = {rid for rid, _ in read_fasta(Path(args.targets))}
        clashes = sorted(target_ids & {rid for rid, _ in binder_records})
        if clashes:
            logging.error(
                "binder and target records share the id(s) %s. Every pair stages both "
                "chains' files side by side, so the names must differ; rename the "
                "--targets record(s) or the complex file(s).",
                ",".join(clashes),
            )
            return 1

    Path(args.binders_out).write_text(
        "".join(f">{rid}\n{seq}\n" for rid, seq in binder_records)
    )
    Path(args.chains_out).write_text(
        "binder_id\tcomplex\tbinder_chain\ttarget_chain\n"
        + "".join("\t".join(row) + "\n" for row in chain_rows)
    )
    logging.info("wrote %d binder record(s)", len(binder_records))

    if args.targets_out:
        # Deduplicate by sequence: the usual case is many binders against one target,
        # which must become one target record and not one per complex, or the pulldown
        # folds every binder against every copy of the same molecule.
        by_seq: Dict[str, str] = {}
        for cid, seq in complex_targets:
            by_seq.setdefault(seq, cid)
        if len(by_seq) > 1:
            logging.error(
                "the complexes hold %d distinct target sequences (%s). A pulldown "
                "folds every binder against every target, so deriving them here would "
                "fold %d pairs, most of them binder-against-someone-else's-target. "
                "Pass --targets to say which target(s) to fold against.",
                len(by_seq),
                ", ".join(sorted(by_seq.values())),
                len(by_seq) * len(binder_records),
            )
            return 1
        # Named 'target', not after the complex it came from: that name is already the
        # binder record's id, and the pulldown stages both chains' files side by side.
        Path(args.targets_out).write_text(
            "".join(f">target\n{seq}\n" for seq in by_seq)
        )
        logging.info("wrote %d target record(s)", len(by_seq))

    return 0


if __name__ == "__main__":
    sys.exit(main())
