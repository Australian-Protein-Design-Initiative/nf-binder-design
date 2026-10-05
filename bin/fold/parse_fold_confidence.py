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

Handles AF2 / RF3 / Protenix / AF3 / OpenFold3 / ESMFold2 (Boltz keeps its own bin/parse_boltz_confidence.py,
which is shared with boltz_pulldown.nf). Emits the canonical, cross-engine
column schema: equivalent scores share a column name, plddt is rescaled to 0-1,
and asymmetric per-chain-pair scores are intentionally dropped (only overall /
"main" values are reported). Columns an engine does not report are left blank.

Sources per tool (verified against example fold-multimer results):
  af2 / af2_mono
           --pkl result_model_N.pkl (ptm, iptm, ranking_confidence, plddt[0-100])
           optional --ipsae-tsv (ipsae.py output; Type==min row)
           af2_mono reuses the AF2 pickle parser; the TSV tool column is af2_mono
  rf3      --json *_summary_confidences.json
           optional --ipsae-tsv from *_confidences.json + *_model.cif
  protenix --json *_summary_confidence_sample_N.json (plddt is 0-100)
           optional --ipsae-tsv from *_full_data_sample_N.json + *_sample_N.cif
  af3      --json *_seed-S_sample-N_summary_confidences.json (no pLDDT in it)
           optional --full-json *_confidences.json (atom_plddts [0-100], pae)
           optional --ipsae-tsv from *_confidences.json + *_model.cif
  openfold3
           --json *_seed_S_sample_N_confidences_aggregated.json
           optional --full-json *_confidences.json (per-atom plddt, pae)
           optional --ipsae-tsv from *_confidences.json + *_model.cif
  esmfold2 / esmfold2_fast
           --json *_seed_S_sample_N_summary_confidences.json (plddt already 0-1;
           no ranking score - ESMFold2 reports none)
           optional --full-json *_confidences.json (pae)
           optional --ipsae-tsv from *_confidences.json + *_model.cif
"""

import argparse
import csv
import json
import sys

# Canonical column order for every fold engine's row (see plans/fold-nf-scores-tsv.md).
# batch/msa_depth distinguish repeated (tool, id, model) rows across
# --n_predictions batches and --msa_subsample depth jobs; blank when fold.nf
# ran a single, unbatched, full-depth job.
COLUMNS = [
    "tool", "id", "model", "batch", "msa_depth", "original_file", "predictions_file",
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
    row = {
        "ranking_score": _num(d.get("ranking_score")),
        "ptm": _num(d.get("ptm")),
        "iptm": _num(d.get("iptm")),
        "plddt": _num(d.get("overall_plddt")),  # already 0-1
        "pae": _num(d.get("overall_pae")),
        "pde": _num(d.get("overall_pde")),
        "has_clash": d.get("has_clash"),
    }
    if args.ipsae_tsv:
        row.update(_read_ipsae_min(args.ipsae_tsv))
    return row


def parse_protenix(args):
    with open(args.json) as f:
        d = json.load(f)
    plddt = _num(d.get("plddt"))
    row = {
        "ranking_score": _num(d.get("ranking_score")),
        "ptm": _num(d.get("ptm")),
        "iptm": _num(d.get("iptm")),
        "plddt": plddt / 100.0 if plddt is not None else None,  # protenix is 0-100
        "pde": _num(d.get("gpde")),
        "has_clash": d.get("has_clash"),
    }
    if args.ipsae_tsv:
        row.update(_read_ipsae_min(args.ipsae_tsv))
    return row


def parse_af3(args):
    import numpy as np

    with open(args.json) as f:
        d = json.load(f)
    has_clash = d.get("has_clash")
    row = {
        "ranking_score": _num(d.get("ranking_score")),
        "ptm": _num(d.get("ptm")),
        "iptm": _num(d.get("iptm")),
        "has_clash": bool(has_clash) if has_clash is not None else None,
    }
    if args.full_json:
        with open(args.full_json) as f:
            full = json.load(f)
        plddts = full.get("atom_plddts")
        if plddts:
            row["plddt"] = float(np.mean(plddts)) / 100.0  # af3 is 0-100
        pae = full.get("pae")
        if pae:
            row["pae"] = float(np.nanmean(np.asarray(pae, dtype=float)))
    if args.ipsae_tsv:
        row.update(_read_ipsae_min(args.ipsae_tsv))
    return row


def _plddt_fraction(x):
    """Scale a mean pLDDT to 0-1, whichever scale the engine reports."""
    x = _num(x)
    if x is None:
        return None
    return x / 100.0 if x > 1.0 else x


def parse_openfold3(args):
    import numpy as np

    with open(args.json) as f:
        d = json.load(f)
    has_clash = _num(d.get("has_clash"))
    row = {
        "ranking_score": _num(d.get("sample_ranking_score")),
        "ptm": _num(d.get("ptm")),
        "iptm": _num(d.get("iptm")),
        "plddt": _plddt_fraction(d.get("avg_plddt")),
        "pde": _num(d.get("gpde")),
        "has_clash": bool(has_clash) if has_clash is not None else None,
    }
    if args.full_json:
        with open(args.full_json) as f:
            full = json.load(f)
        pae = full.get("pae")
        if pae:
            row["pae"] = float(np.nanmean(np.asarray(pae, dtype=float)))
    if args.ipsae_tsv:
        row.update(_read_ipsae_min(args.ipsae_tsv))
    return row


def parse_esmfold2(args):
    import numpy as np

    with open(args.json) as f:
        d = json.load(f)
    row = {
        # ESMFold2 reports no ranking score; rank on iptm / plddt / ipsae instead.
        "ptm": _num(d.get("ptm")),
        "iptm": _num(d.get("iptm")),
        "plddt": _plddt_fraction(d.get("plddt")),  # already 0-1 from run_esmfold2.py
    }
    if args.full_json:
        with open(args.full_json) as f:
            full = json.load(f)
        pae = full.get("pae")
        if pae:
            row["pae"] = float(np.nanmean(np.asarray(pae, dtype=float)))
    if args.ipsae_tsv:
        row.update(_read_ipsae_min(args.ipsae_tsv))
    return row


PARSERS = {
    "af2": parse_af2,
    "af2_mono": parse_af2,  # same pickle schema; distinct tool tag in the TSV
    "rf3": parse_rf3,
    "protenix": parse_protenix,
    "af3": parse_af3,
    "openfold3": parse_openfold3,
    "esmfold2": parse_esmfold2,
    "esmfold2_fast": parse_esmfold2,  # same run_esmfold2.py output; distinct tool tag
}


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
    p.add_argument("--batch", default="", help="fold.nf batch index (meta.fold_batch / af2_run); blank if unbatched")
    p.add_argument("--msa-depth", default="", help="fold.nf MSA subsample depth tag (meta.msa_depth_tag); blank if full-depth")
    p.add_argument("--original-file", default="", help="engine-native structure filename")
    p.add_argument("--predictions-file", default="", help="renamed name in fold/predictions/")
    p.add_argument("--json", help="confidence/summary JSON (rf3, protenix, af3, openfold3, esmfold2[_fast])")
    p.add_argument("--full-json", help="AF3 / OpenFold3 / ESMFold2 *_confidences.json (per-atom pLDDT, PAE)")
    p.add_argument("--pkl", help="AF2 result_model_N.pkl")
    p.add_argument("--ipsae-tsv", help="ipsae.py output TSV (optional; Type==min row)")
    p.add_argument("--no-header", action="store_true", help="omit the header line (for concatenation)")
    args = p.parse_args()

    scores = PARSERS[args.tool](args)
    row = {c: "" for c in COLUMNS}
    row.update({
        "tool": args.tool,
        "id": args.id,
        "model": args.model,
        "batch": args.batch,
        "msa_depth": args.msa_depth,
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
