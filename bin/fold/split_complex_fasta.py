#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Split a multi-record complex FASTA into one single-record FASTA per chain, for
per-chain MSA search in fold.nf's multimer path (see
plans/fold-nf-multimer-paired-msa.md §4.1 step 1).

Each record becomes chain A, B, C, ... in file order. Output files are named
`<name>.chain_<CHAIN>.fasta` so each chain's stem is unique and stable - the AF2
jackhmmer MSA stage names its output directory after the FASTA stem, and the
chain letter keeps per-chain searches of the same complex from colliding.
"""

import argparse
import string
import sys
from pathlib import Path
from typing import List, Tuple


def parse_fasta_records(fasta_path: Path) -> List[Tuple[str, str]]:
    """Return (header, sequence) per record, in file order."""
    records: List[Tuple[str, str]] = []
    header = None
    seq: List[str] = []
    for line in fasta_path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = s[1:]
            seq = []
        else:
            seq.append(s)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="Multi-record complex FASTA")
    parser.add_argument("--name", required=True, help="Complex name (output file stem prefix)")
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    records = parse_fasta_records(Path(args.fasta))
    if not records:
        print(f"No FASTA records found in {args.fasta}", file=sys.stderr)
        return 1

    chain_ids = list(string.ascii_uppercase)
    if len(records) > len(chain_ids):
        print(f"{args.fasta} has more than {len(chain_ids)} chains; not supported", file=sys.stderr)
        return 1

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for chain_id, (header, seq) in zip(chain_ids, records):
        out = outdir / f"{args.name}.chain_{chain_id}.fasta"
        # Keep the original header but guarantee a chain-tagged, unique query id
        # as the first token (AF2 keys its output dir on the FASTA stem, not the
        # header, but a clean id helps downstream logs).
        out.write_text(f">{args.name}.chain_{chain_id}\n{seq}\n")
    print(f"Wrote {len(records)} per-chain FASTA(s) to {outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
