#!/usr/bin/env python3
# /// script
# requires-python = ">=3.8"
# ///

"""
Merge per-tool fold score TSVs into the master fold_scores.tsv.

AF2/RF3/Protenix per-tool tables are already in the canonical schema (produced
by bin/parse_fold_confidence.py) and are concatenated as-is. Boltz keeps its own
native boltz_fold_scores.tsv (from bin/parse_boltz_confidence.py, shared with
boltz_pulldown.nf), so its columns are mapped onto the canonical schema here.

Canonical schema: see COLUMNS below. plddt is 0-1; asymmetric per-chain-pair
scores are omitted; blanks where an engine does not report a metric.
"""

import argparse
import csv
import sys

COLUMNS = [
    "tool", "id", "model", "original_file", "predictions_file",
    "ranking_score", "ptm", "iptm", "plddt", "pae", "pde", "has_clash",
    "ipsae", "ipsae_d0chn", "ipsae_d0dom", "pdockq", "pdockq2", "lis",
]

# Boltz native (bin/parse_boltz_confidence.py) column -> canonical column.
BOLTZ_MAP = {
    "id": "id",
    "model": "model",
    "original_file": "original_file",
    "predictions_file": "predictions_file",
    "confidence_score": "ranking_score",
    "ptm": "ptm",
    "iptm": "iptm",
    "complex_plddt": "plddt",
    "complex_pde": "pde",
    "ipsae_min": "ipsae",
}


def _read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _canonical_rows(rows_in):
    return [{c: (r.get(c, "") or "") for c in COLUMNS} for r in rows_in]


def _boltz_rows(rows_in):
    rows = []
    for r in rows_in:
        out = {c: "" for c in COLUMNS}
        out["tool"] = "boltz"
        for native, canon in BOLTZ_MAP.items():
            if native in r and r[native] not in (None, ""):
                out[canon] = r[native]
        rows.append(out)
    return rows


def _rows_for(path):
    """Auto-detect a per-tool TSV's schema and normalize it to canonical rows."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    rows_in = _read_tsv(path)
    if header and header[0] == "tool":
        return _canonical_rows(rows_in)          # af2 / rf3 / protenix (already canonical)
    if "confidence_score" in header:
        return _boltz_rows(rows_in)              # native boltz_fold_scores.tsv
    sys.stderr.write(f"merge_fold_scores: unrecognized schema, skipping {path}\n")
    return []


def main():
    p = argparse.ArgumentParser(description="Merge per-tool fold score TSVs into the master table")
    p.add_argument("--input", action="append", default=[],
                   help="Path to a per-tool scores TSV (schema auto-detected); repeatable")
    p.add_argument("-o", "--output", help="Output TSV (default stdout)")
    args = p.parse_args()

    rows = []
    for path in args.input:
        rows.extend(_rows_for(path))

    # Stable, readable order: by tool then id then model.
    rows.sort(key=lambda r: (r["tool"], r["id"], str(r["model"])))

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in COLUMNS])
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()
