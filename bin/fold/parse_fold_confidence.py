#!/usr/bin/env python3
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "numpy",
# ]
# ///

"""
Parse a single fold.nf prediction's confidence output into ONE normalized TSV
row (written to stdout, with header) for the master fold_scores table.

Handles AF2 / RF3 / Protenix (Boltz keeps its own bin/parse_boltz_confidence.py,
which is shared with boltz_pulldown.nf). Emits the canonical, cross-engine
column schema: equivalent scores share a column name, plddt is rescaled to 0-1,
and asymmetric per-chain-pair scores are intentionally dropped (only overall /
"main" values are reported). Columns an engine does not report are left blank.

Sources per tool (verified against example fold-multimer results):
  af2      --pkl result_model_N.pkl (ptm, iptm, ranking_confidence, plddt[0-100])
           optional --ipsae-tsv (ipsae.py output; Type==min row)
  rf3      --json *_summary_confidences.json
  protenix --json *_summary_confidence_sample_N.json (plddt is 0-100)
"""

import argparse
import csv
import json
import sys

# Canonical column order for every fold engine's row (see plans/fold-nf-scores-tsv.md).
COLUMNS = [
    "tool", "id", "model", "original_file", "predictions_file",
    "ranking_score", "ptm", "iptm", "plddt", "pae", "pde", "has_clash",
    "ipsae", "ipsae_d0chn", "ipsae_d0dom", "pdockq", "pdockq2", "lis",
]


def _num(x):
    """Coerce to float, or None if missing/non-numeric."""
    if x is None:
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def parse_af2(args):
    import pickle

    import numpy as np

    with open(args.pkl, "rb") as f:
        d = pickle.load(f)
    plddt = d.get("plddt")
    mean_plddt = float(np.mean(plddt)) / 100.0 if plddt is not None else None
    row = {
        "ranking_score": _num(d.get("ranking_confidence")),
        "ptm": _num(d.get("ptm")),
        "iptm": _num(d.get("iptm")),
        "plddt": mean_plddt,
    }
    if args.ipsae_tsv:
        row.update(_read_ipsae_min(args.ipsae_tsv))
    return row


def _read_ipsae_min(path):
    """Return the ipsae.py Type==min row's interface metrics, normalized.

    ipsae.py emits a leading blank line before the header and is whitespace
    (not strictly tab) delimited, so parse on any-whitespace after dropping
    blank lines.
    """
    with open(path) as f:
        lines = [ln for ln in f if ln.strip()]
    if not lines:
        return {}
    header = lines[0].split()
    rows = [dict(zip(header, ln.split())) for ln in lines[1:]]
    min_row = next((r for r in rows if r.get("Type") == "min"), None)
    if min_row is None:
        return {}
    return {
        "ipsae": _num(min_row.get("ipSAE")),
        "ipsae_d0chn": _num(min_row.get("ipSAE_d0chn")),
        "ipsae_d0dom": _num(min_row.get("ipSAE_d0dom")),
        "pdockq": _num(min_row.get("pDockQ")),
        "pdockq2": _num(min_row.get("pDockQ2")),
        "lis": _num(min_row.get("LIS")),
    }


def parse_rf3(args):
    with open(args.json) as f:
        d = json.load(f)
    return {
        "ranking_score": _num(d.get("ranking_score")),
        "ptm": _num(d.get("ptm")),
        "iptm": _num(d.get("iptm")),
        "plddt": _num(d.get("overall_plddt")),  # already 0-1
        "pae": _num(d.get("overall_pae")),
        "pde": _num(d.get("overall_pde")),
        "has_clash": d.get("has_clash"),
    }


def parse_protenix(args):
    with open(args.json) as f:
        d = json.load(f)
    plddt = _num(d.get("plddt"))
    return {
        "ranking_score": _num(d.get("ranking_score")),
        "ptm": _num(d.get("ptm")),
        "iptm": _num(d.get("iptm")),
        "plddt": plddt / 100.0 if plddt is not None else None,  # protenix is 0-100
        "pde": _num(d.get("gpde")),
        "has_clash": d.get("has_clash"),
    }


PARSERS = {"af2": parse_af2, "rf3": parse_rf3, "protenix": parse_protenix}


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def main():
    p = argparse.ArgumentParser(description="Parse a fold prediction's confidence to a normalized TSV row")
    p.add_argument("--tool", required=True, choices=list(PARSERS))
    p.add_argument("--id", required=True)
    p.add_argument("--model", required=True, help="per-structure index label")
    p.add_argument("--original-file", default="", help="engine-native structure filename")
    p.add_argument("--predictions-file", default="", help="renamed name in fold/predictions/")
    p.add_argument("--json", help="confidence/summary JSON (rf3, protenix)")
    p.add_argument("--pkl", help="AF2 result_model_N.pkl")
    p.add_argument("--ipsae-tsv", help="AF2 ipsae.py output TSV (optional)")
    p.add_argument("--no-header", action="store_true", help="omit the header line (for concatenation)")
    args = p.parse_args()

    scores = PARSERS[args.tool](args)
    row = {c: "" for c in COLUMNS}
    row.update({
        "tool": args.tool,
        "id": args.id,
        "model": args.model,
        "original_file": args.original_file,
        "predictions_file": args.predictions_file,
    })
    for k, v in scores.items():
        row[k] = _fmt(v)

    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    if not args.no_header:
        w.writerow(COLUMNS)
    w.writerow([row[c] for c in COLUMNS])


if __name__ == "__main__":
    main()
