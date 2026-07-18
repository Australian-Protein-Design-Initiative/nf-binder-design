#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Derive a single merged a3m alignment from an AlphaFold2 jackhmmer/hhblits
"msas/" directory (as written by modules/local/af2/alphafold2_jackhmmer_msa.nf)
for Boltz-2 and RF3, which each read one a3m per chain directly.

Ground truth (probed against alphafold/data/pipeline.py in the AF2 CUDA-12
custom container, 2026-07-17): a monomer msas/ dir contains
bfd_uniref_hits.a3m (already a3m), uniref90_hits.sto and mgnify_hits.sto
(Stockholm), plus pdb_hits.hhr (template hits, not an alignment - not used
here).

Each source's alignment is independently built against the same query, so a3m
match columns (upper-case letters + '-') all have length == len(query) - the
sources can therefore be concatenated as separate a3m "blocks" sharing one
query, the same way ColabFold's own combined a3m works internally.
AlphaFold's a3m/Stockholm parsers discard raw insertion characters (they keep
only a per-residue deletion *count*, not the deleted residues themselves), so
this merge is not byte-for-byte reversible - re-emitted sequences carry match
columns only (no lower-case inserts). This is a documented quality caveat,
not a bug: the merged a3m still carries the same match-column
co-evolutionary signal Boltz/RF3 use for their own MSA featurization.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Set

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

# See colabfold_a3m_to_af2_msas.py for why this is needed (PYTHONPATH is
# empty in the AF2 container; only run_alphafold.py itself gets /app/alphafold
# on sys.path, via Python's script-relative sys.path[0] insertion).
sys.path.insert(0, "/app/alphafold")

from alphafold.data import parsers  # noqa: E402


def msa_to_a3m_lines(msa, seen: Set[str]) -> List[str]:
    lines: List[str] = []
    for seq, desc in zip(msa.sequences, msa.descriptions):
        if seq in seen:
            continue
        seen.add(seq)
        lines.append(f">{desc}")
        lines.append(seq)
    return lines


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--msas-dir", required=True, help="AF2 msas/ directory (contains *.sto/*.a3m/*.hhr)")
    parser.add_argument("--output", required=True, help="Merged a3m output path")
    args = parser.parse_args()

    msas_dir = Path(args.msas_dir)
    bfd_a3m_path = msas_dir / "bfd_uniref_hits.a3m"
    uniref90_sto_path = msas_dir / "uniref90_hits.sto"
    mgnify_sto_path = msas_dir / "mgnify_hits.sto"

    bfd_msa = parsers.parse_a3m(bfd_a3m_path.read_text()) if bfd_a3m_path.exists() else None
    uniref90_msa = parsers.parse_stockholm(uniref90_sto_path.read_text()) if uniref90_sto_path.exists() else None
    mgnify_msa = parsers.parse_stockholm(mgnify_sto_path.read_text()) if mgnify_sto_path.exists() else None

    for label, msa, path in (
        ("bfd_uniref_hits.a3m", bfd_msa, bfd_a3m_path),
        ("uniref90_hits.sto", uniref90_msa, uniref90_sto_path),
        ("mgnify_hits.sto", mgnify_msa, mgnify_sto_path),
    ):
        if msa is None:
            log.warning("Expected AF2 MSA source not found, skipping: %s", path)

    sources = [m for m in (bfd_msa, uniref90_msa, mgnify_msa) if m is not None]
    if not sources:
        log.error("No AF2 MSA sources found under %s", msas_dir)
        return 1

    # Query-first invariant: the first source's first record is the query
    # itself (jackhmmer/hhblits always emit the query as the first hit).
    seen: Set[str] = set()
    query_seq, query_desc = sources[0].sequences[0], sources[0].descriptions[0]
    seen.add(query_seq)
    lines = [f">{query_desc}", query_seq]

    for msa in sources:
        lines.extend(msa_to_a3m_lines(msa, seen))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    log.info("Wrote merged a3m with %d sequences (from %d source(s)) to %s", len(seen), len(sources), out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
