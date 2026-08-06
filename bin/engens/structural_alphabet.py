#!/usr/bin/env python
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "biotite<=1.6",
#     "numpy",
# ]
# ///
#
# Encode protein structures as structural-alphabet sequences (3Di, Protein Blocks).
# Usage:
#   uv run bin/engens/structural_alphabet.py structures/ --alphabet 3di -o 3di.fasta
#   uv run bin/engens/structural_alphabet.py a.pdb b.pdb --alphabet pb --npz codes.npz -o -
#
# TODO: --alphabet 3dn (Zerefa et al., Bioinformatics 2025 / spetti/structure_comparison)
# is the planned extension point. Blockers before implementing:
# - Licence: the upstream repo has no LICENSE file, so vendoring is not safe.
# - Reference data: transition_mtx.npy (~8 MB), MI_centers.npy, graph_clusters_blosum.npy
#   should be baked into the EnGens container (runtime download fails on air-gapped HPC).
# - Performance: upstream getOneHot uses list.index() over ~4001 bins per neighbour;
#   needs a dict lookup / vectorising before running over large ensembles.
# Once sequences exist, the same one-hot -> substitution-matrix embedding path applies.

import argparse
import io
import logging
import sys
from pathlib import Path
from typing import Iterable, Optional

import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger("structural_alphabet")

STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif", ".pdb1"}
ALPHABETS = ("3di", "pb")


def collect_structure_paths(inputs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for item in inputs:
        p = Path(item)
        if p.is_dir():
            found = sorted(
                q
                for q in p.iterdir()
                if q.suffix.lower() in STRUCTURE_SUFFIXES
                or (len(q.suffixes) >= 2 and q.suffixes[-2].lower() in STRUCTURE_SUFFIXES)
            )
            if not found:
                log.warning("No structure files found in directory %s", p)
            paths.extend(found)
        elif p.is_file():
            paths.append(p)
        else:
            log.warning("Skipping missing path %s", p)
    return paths


def encode_structure(
    path: Path,
    alphabet: str,
) -> tuple[str, np.ndarray, np.ndarray]:
    """Return (sequence_string, code_array, chain_start_residue_indices)."""
    atoms = strucio.load_structure(str(path), model=1)
    if isinstance(atoms, struc.AtomArrayStack):
        atoms = atoms[0]
    atoms = atoms[struc.filter_amino_acids(atoms)]
    if atoms.array_length() == 0:
        raise ValueError(f"No amino-acid atoms in {path.name}")

    if alphabet == "3di":
        from biotite.structure.alphabet import to_3di

        sequences, chain_starts = to_3di(atoms)
    elif alphabet == "pb":
        from biotite.structure.alphabet import to_protein_blocks

        sequences, chain_starts = to_protein_blocks(atoms)
    else:
        raise ValueError(f"Unsupported alphabet: {alphabet}")

    symbols: list[str] = []
    codes: list[np.ndarray] = []
    for seq in sequences:
        symbols.extend(seq.symbols)
        codes.append(np.asarray(seq.code, dtype=np.int32))

    if not codes:
        raise ValueError(f"No chains encoded for {path.name}")

    code_arr = np.concatenate(codes)
    seq_str = "".join(symbols)
    # Residue-indexed chain starts (0-based within the concatenated sequence).
    res_starts = np.cumsum([0] + [len(c) for c in codes[:-1]], dtype=np.int32)
    return seq_str, code_arr, res_starts


def count_undefined(alphabet: str, seq: str) -> int:
    """Count residues in the alphabet's undefined/fallback state.

    biotite's 3Di Encoder fills missing-backbone residues with the masked fill
    value that decodes to 'd', which is also a valid 3Di letter, so this is only
    a soft signal (cannot distinguish true 'd' from missing atoms).
    Protein Blocks uses 'x' as the undefined symbol.
    """
    if alphabet == "3di":
        return seq.count("d") + seq.count("D")
    if alphabet == "pb":
        return seq.count("x") + seq.count("X")
    return 0


def write_fasta(handle: io.TextIOBase, labels: Iterable[str], seqs: Iterable[str]) -> None:
    for label, seq in zip(labels, seqs):
        handle.write(f">{label}\n")
        for i in range(0, len(seq), 80):
            handle.write(seq[i : i + 80] + "\n")


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Encode protein structures as structural-alphabet sequences.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Structure files (.pdb/.cif) and/or directories containing them",
    )
    parser.add_argument(
        "--alphabet",
        choices=ALPHABETS,
        default="3di",
        help="Structural alphabet to encode (default: 3di). "
        "'3dn' is planned but not yet implemented (see module TODO).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="FASTA output path, or '-' for stdout (default: -)",
    )
    parser.add_argument(
        "--npz",
        default=None,
        metavar="PATH",
        help="Optional .npz with arrays: labels, codes (object array of int arrays), "
        "chain_starts (object array of int arrays), sequences",
    )
    args = parser.parse_args(argv)

    if args.alphabet == "3dn":
        log.error(
            "Alphabet '3dn' is not implemented yet (licence / reference-data / "
            "performance blockers; see the TODO at the top of this script)."
        )
        return 2

    paths = collect_structure_paths(args.inputs)
    if not paths:
        log.error("No structure files to encode")
        return 1

    labels: list[str] = []
    seqs: list[str] = []
    codes: list[np.ndarray] = []
    chain_starts: list[np.ndarray] = []
    n_undefined_total = 0

    for path in paths:
        try:
            seq, code, starts = encode_structure(path, args.alphabet)
        except Exception as exc:  # noqa: BLE001
            log.warning("Failed to encode %s: %s", path.name, exc)
            continue
        n_undef = count_undefined(args.alphabet, seq)
        n_undefined_total += n_undef
        if n_undef:
            log.warning(
                "%s: %d/%d residues are the undefined/fallback state for %s",
                path.name,
                n_undef,
                len(seq),
                args.alphabet,
            )
        labels.append(path.stem)
        seqs.append(seq)
        codes.append(code)
        chain_starts.append(starts)
        log.info("%s: length %d (%d chain(s))", path.name, len(seq), len(starts))

    if not seqs:
        log.error("No structures encoded successfully")
        return 1

    if n_undefined_total:
        log.warning(
            "Total undefined/fallback residues across ensemble: %d",
            n_undefined_total,
        )

    if args.output == "-":
        write_fasta(sys.stdout, labels, seqs)
    else:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8") as fh:
            write_fasta(fh, labels, seqs)
        log.info("Wrote FASTA to %s", out)

    if args.npz:
        np.savez_compressed(
            args.npz,
            labels=np.asarray(labels, dtype=object),
            sequences=np.asarray(seqs, dtype=object),
            codes=np.asarray(codes, dtype=object),
            chain_starts=np.asarray(chain_starts, dtype=object),
            alphabet=np.asarray(args.alphabet),
        )
        log.info("Wrote arrays to %s", args.npz)

    return 0


if __name__ == "__main__":
    sys.exit(main())
