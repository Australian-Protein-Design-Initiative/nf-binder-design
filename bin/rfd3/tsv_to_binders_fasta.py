#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pandas",
# ]
# ///

"""
Write binders.fasta from combined_scores.tsv using id and sequence columns.
Optionally add score columns to FASTA headers (like pdb_to_fasta.py).
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def main() -> int:
    parser = argparse.ArgumentParser(description="Write binders.fasta from combined_scores TSV")
    parser.add_argument("tsv", type=Path, help="combined_scores.tsv")
    parser.add_argument("-o", "--output", type=Path, default=Path("binders.fasta"), help="Output FASTA path")
    parser.add_argument(
        "--scores",
        default="pair_pae_min,plddt",
        help="Comma-separated score columns to add to headers (default: pair_pae_min,plddt)",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    if "sequence" not in df.columns:
        log.error("TSV has no 'sequence' column")
        return 1
    id_col = "id" if "id" in df.columns else "filename"
    if id_col == "filename":
        df["_header"] = df["filename"].str.replace(r"\.cif$", "", regex=True)
    else:
        df["_header"] = df[id_col].astype(str)

    score_cols = [c.strip() for c in args.scores.split(",") if c.strip()]
    for c in score_cols:
        if c in df.columns:
            df["_header"] = df["_header"] + "|" + c + "=" + df[c].astype(str)
        else:
            log.warning("Score column %s not in TSV", c)

    with open(args.output, "w") as f:
        for _, row in df.iterrows():
            seq = row.get("sequence")
            if pd.isna(seq) or not str(seq).strip():
                continue
            header = row["_header"].strip()
            if not header:
                continue
            f.write(f">{header}\n")
            f.write(f"{seq}\n")

    log.info("Wrote %s", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
