#!/usr/bin/env python3
# /// script
# requires-python = ">=3.8"
# ///
"""
Assemble an AF2 multimer precomputed-MSA tree for a target+binder pulldown pair.

Layout written (matches AlphaFold2 multimer --use_precomputed_msas):

  <pair_id>/
    msas/
      A/   # target: from jackhmmer dir, ColabFold/mmseqs2 a3m, or query-only
      B/   # binder: always query-only (pulldown binders have no useful homologs)

Required per-chain files (created as query-only when absent):
  uniref90_hits.sto, mgnify_hits.sto, bfd_uniref_hits.a3m (or bfd_uniclust_hits.a3m),
  uniprot_hits.sto, pdb_hits.hhr

features.pkl is deliberately omitted so AF2 rebuilds features for the pair.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path
from typing import List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

STO_NAMES = ("uniref90_hits.sto", "mgnify_hits.sto", "uniprot_hits.sto")
A3M_CANDIDATES = ("bfd_uniref_hits.a3m", "bfd_uniclust_hits.a3m")
HHR_NAME = "pdb_hits.hhr"


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    parts: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header = line[1:].strip()
            parts = []
        else:
            parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def parse_a3m(path: Path) -> List[Tuple[str, str]]:
    """Return (header, sequence) pairs; sequences keep a3m match columns."""
    return parse_fasta(path)  # same >header / seq format


def write_query_sto(path: Path, seq_id: str, sequence: str) -> None:
    path.write_text(
        "# STOCKHOLM 1.0\n"
        f"{seq_id} {sequence}\n"
        "//\n"
    )


def write_query_a3m(path: Path, seq_id: str, sequence: str) -> None:
    path.write_text(f">{seq_id}\n{sequence}\n")


def write_empty_hhr(path: Path, seq_id: str, sequence: str) -> None:
    path.write_text(
        f"Query         {seq_id}\n"
        f"Match_columns {len(sequence)}\n"
        "No_of_seqs    1 out of 1\n"
        "Neff          1.0\n"
        "Searched_HMMs 0\n"
        "Date          \n"
        "Command       \n\n"
        " No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM\n"
    )


def a3m_to_stockholm(records: List[Tuple[str, str]]) -> str:
    """Convert a3m records to Stockholm (match columns only: upper-case + '-')."""
    lines = ["# STOCKHOLM 1.0"]
    for i, (hdr, seq) in enumerate(records):
        seq_id = hdr.split()[0] if hdr else f"seq{i}"
        # Drop a3m insertions (lower-case) for Stockholm match columns
        match = "".join(c for c in seq if c.isupper() or c == "-")
        lines.append(f"{seq_id} {match}")
    lines.append("//")
    return "\n".join(lines) + "\n"


def is_dummy_path(path: Optional[Path]) -> bool:
    if path is None or not path.exists():
        return True
    name = path.name
    return name in ("empty", "empty_target_msa", "dummy") or (
        path.is_file() and path.stat().st_size == 0
    )


def find_chain_dir(msa_root: Path) -> Optional[Path]:
    """Locate a usable target MSA directory (msas/A or flat msas/)."""
    if not msa_root.exists():
        return None
    for candidate in (msa_root / "msas" / "A", msa_root / "A"):
        if candidate.is_dir() and any(candidate.iterdir()):
            return candidate
    flat = msa_root / "msas"
    if flat.is_dir() and any(flat.iterdir()):
        files = [p for p in flat.iterdir() if p.is_file()]
        if files:
            return flat
    return None


def write_chain_from_a3m(
    a3m_path: Path,
    dest_dir: Path,
    seq_id: str,
    sequence: str,
) -> None:
    """Materialise AF2 per-chain MSA files from a ColabFold/mmseqs2 a3m."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    records = parse_a3m(a3m_path)
    if not records:
        log.warning("empty a3m %s; falling back to query-only for chain A", a3m_path)
        copy_or_query(None, dest_dir, seq_id, sequence, force_query_only=True)
        return

    # Full MSA as the hhblits-style a3m source AF2 always merges.
    shutil.copy2(a3m_path, dest_dir / A3M_CANDIDATES[0])
    # Same alignment as Stockholm for uniref90 (depth preserved once in merge).
    (dest_dir / "uniref90_hits.sto").write_text(a3m_to_stockholm(records))
    # Query-only for the other sources so we do not triple-count the same hits.
    write_query_sto(dest_dir / "mgnify_hits.sto", seq_id, sequence)
    write_query_sto(dest_dir / "uniprot_hits.sto", seq_id, sequence)
    write_empty_hhr(dest_dir / HHR_NAME, seq_id, sequence)
    log.info(
        "chain A from a3m %s (%d sequences) -> %s",
        a3m_path.name,
        len(records),
        dest_dir,
    )


def copy_or_query(
    src_dir: Optional[Path],
    dest_dir: Path,
    seq_id: str,
    sequence: str,
    *,
    force_query_only: bool,
) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)

    for name in STO_NAMES:
        dest = dest_dir / name
        src = src_dir / name if src_dir and not force_query_only else None
        if src is not None and src.is_file():
            shutil.copy2(src, dest)
        else:
            write_query_sto(dest, seq_id, sequence)

    a3m_copied = False
    if src_dir and not force_query_only:
        for name in A3M_CANDIDATES:
            src = src_dir / name
            if src.is_file():
                shutil.copy2(src, dest_dir / name)
                a3m_copied = True
                break
    if not a3m_copied:
        write_query_a3m(dest_dir / A3M_CANDIDATES[0], seq_id, sequence)

    hhr_dest = dest_dir / HHR_NAME
    hhr_src = src_dir / HHR_NAME if src_dir and not force_query_only else None
    if hhr_src is not None and hhr_src.is_file():
        shutil.copy2(hhr_src, hhr_dest)
    else:
        write_empty_hhr(hhr_dest, seq_id, sequence)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pair-fasta", required=True, type=Path, help="Two-record FASTA (target then binder)")
    p.add_argument("--pair-id", required=True, help="Output directory name (= meta.id)")
    p.add_argument(
        "--target-msa-dir",
        type=Path,
        default=None,
        help="Optional prior AF2 MSA dir for the target (jackhmmer monomer/multimer layout)",
    )
    p.add_argument(
        "--target-a3m",
        type=Path,
        default=None,
        help="Optional ColabFold/mmseqs2 a3m for the target (used when no jackhmmer MSA dir)",
    )
    p.add_argument("-o", "--output-dir", type=Path, required=True, help="Parent dir; writes <pair-id>/ under it")
    args = p.parse_args()

    records = parse_fasta(args.pair_fasta)
    if len(records) < 2:
        log.error("pair FASTA needs >= 2 records, got %d", len(records))
        return 1
    (hdr_a, seq_a), (hdr_b, seq_b) = records[0], records[1]
    id_a = hdr_a.split()[0] if hdr_a else "A"
    id_b = hdr_b.split()[0] if hdr_b else "B"

    out_root = args.output_dir / args.pair_id
    if out_root.exists():
        shutil.rmtree(out_root)
    msas = out_root / "msas"

    src_a = None
    if args.target_msa_dir and not is_dummy_path(args.target_msa_dir):
        src_a = find_chain_dir(args.target_msa_dir)
        if src_a is None:
            log.warning("no usable MSA files under %s", args.target_msa_dir)

    a3m_a = args.target_a3m if args.target_a3m and not is_dummy_path(args.target_a3m) else None
    if a3m_a is not None and not a3m_a.is_file():
        log.warning("target a3m not a file: %s", a3m_a)
        a3m_a = None

    if src_a is not None:
        copy_or_query(src_a, msas / "A", id_a, seq_a, force_query_only=False)
        a_src = str(src_a)
    elif a3m_a is not None:
        write_chain_from_a3m(a3m_a, msas / "A", id_a, seq_a)
        a_src = f"a3m:{a3m_a.name}"
    else:
        copy_or_query(None, msas / "A", id_a, seq_a, force_query_only=True)
        a_src = "query-only"

    copy_or_query(None, msas / "B", id_b, seq_b, force_query_only=True)

    uni_a = msas / "A" / "uniprot_hits.sto"
    if not uni_a.is_file():
        write_query_sto(uni_a, id_a, seq_a)

    log.info("assembled %s (A from %s, B query-only)", out_root, a_src)
    return 0


if __name__ == "__main__":
    sys.exit(main())
