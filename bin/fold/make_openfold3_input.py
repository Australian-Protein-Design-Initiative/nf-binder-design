#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
OpenFold3 query JSON for fold.nf's OPENFOLD3_FOLD (precomputed-MSA mode).

OpenFold3 is always run with --use-msa-server false, feeding MSAs from the
pipeline's shared MSA stage. Its precomputed-MSA conventions (v0.5.0) drive the
layout written next to the JSON:

  - Each chain's MSAs live in a directory, and the directory name is the chain's
    representative ID. Directories are named by sequence hash, so identical
    chains share one and OpenFold3 treats the complex as a homomer.
  - Files are selected by stem against MSASettings.max_seq_counts. The unpaired
    MSA is colabfold_main.a3m (in the default aln_order). For multimers the
    species-tagged render from `msa_taxonomy.py --tool openfold3` is written as
    uniprot_hits.a3m, which is in the default msas_to_pair but NOT aln_order, so
    OpenFold3 uses it only for online cross-chain pairing.
  - The first a3m row must be the query: ColabFold '#' lines and NULs are dropped
    and row 0 is checked against the chain sequence.

Seeds are not part of the query JSON (OpenFold3 takes them from the runner YAML),
so one JSON serves every batch.
"""

import argparse
import hashlib
import json
import logging
import string
import sys
from pathlib import Path
from typing import Dict, List, Optional

from make_af3_input import CHAIN_IDS, _match_per_chain, clean_a3m, parse_fasta_records

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

MAIN_MSA_NAME = "colabfold_main.a3m"
PAIRING_MSA_NAME = "uniprot_hits.a3m"


def sanitised_name(name: str) -> str:
    """Query key, used by OpenFold3 as the output directory and file prefix."""
    allowed = set(string.ascii_letters + string.digits + "_-.")
    return "".join(c for c in name.replace(" ", "_") if c in allowed)


def msa_dir_name(seq: str) -> str:
    return "msa_" + hashlib.sha1(seq.upper().encode()).hexdigest()[:12]


def make_openfold3_input(
    fasta_path: Path,
    name: str,
    out_dir: Path,
    main_a3m_paths: Optional[List[Path]] = None,
    pairing_a3m_paths: Optional[List[Path]] = None,
) -> dict:
    sequences = parse_fasta_records(fasta_path)
    if not sequences:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    if len(sequences) > len(CHAIN_IDS):
        raise ValueError(f"{len(sequences)} chains exceeds the {len(CHAIN_IDS)} supported")

    n = len(sequences)
    main = _match_per_chain(main_a3m_paths, n, "--a3m")
    pairing = _match_per_chain(pairing_a3m_paths, n, "--pairing-a3m")

    written: Dict[str, int] = {}
    chains = []
    for i, seq in enumerate(sequences):
        cid = CHAIN_IDS[i]
        dname = msa_dir_name(seq)
        if dname not in written:
            d = out_dir / dname
            d.mkdir(parents=True, exist_ok=True)
            (d / MAIN_MSA_NAME).write_text(clean_a3m(main[i], seq, cid))
            if pairing[i] is not None:
                (d / PAIRING_MSA_NAME).write_text(clean_a3m(pairing[i], seq, cid))
            written[dname] = i
        chains.append({
            "molecule_type": "protein",
            "chain_ids": [cid],
            "sequence": seq,
            "main_msa_file_paths": [dname],
        })

    return {"queries": {sanitised_name(name): {"chains": chains}}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="FASTA file (one record per chain)")
    parser.add_argument("--name", required=True, help="Query name (sanitised; output dir name)")
    parser.add_argument(
        "--a3m",
        nargs="+",
        default=None,
        help="Optional unpaired a3m(s): one per chain in record order (query-only if omitted)",
    )
    parser.add_argument(
        "--pairing-a3m",
        dest="pairing_a3m",
        nargs="+",
        default=None,
        help="Optional species-tagged a3m(s) from msa_taxonomy.py --tool openfold3: one per chain",
    )
    parser.add_argument("-o", "--output", required=True, help="Output JSON path (MSA dirs are written alongside)")
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    spec = make_openfold3_input(
        fasta_path=Path(args.fasta),
        name=args.name,
        out_dir=out_path.parent,
        main_a3m_paths=[Path(p) for p in args.a3m] if args.a3m else None,
        pairing_a3m_paths=[Path(p) for p in args.pairing_a3m] if args.pairing_a3m else None,
    )
    with open(out_path, "w") as f:
        json.dump(spec, f, indent=2)
    (query,) = spec["queries"].values()
    log.info("Wrote OpenFold3 query JSON (%d chain(s)) to %s", len(query["chains"]), out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
