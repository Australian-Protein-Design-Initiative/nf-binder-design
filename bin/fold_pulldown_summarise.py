#!/usr/bin/env python3
# /// script
# requires-python = ">=3.8"
# ///
"""
Join fold_scores.tsv with pairs.tsv and emit:

  fold_pulldown_scores.tsv  - one row per predicted structure (+ target, binder)
  fold_pulldown_summary.tsv - one row per (target, binder, tool) with aggregate
                              iptm/ipsae stats, within-tool z-scores, and
                              cross-tool consensus_z.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


def _f(v: str) -> Optional[float]:
    if v is None or v == "":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _mean(xs: List[float]) -> Optional[float]:
    return statistics.mean(xs) if xs else None


def _median(xs: List[float]) -> Optional[float]:
    return statistics.median(xs) if xs else None


def _max(xs: List[float]) -> Optional[float]:
    return max(xs) if xs else None


def _sd(xs: List[float]) -> Optional[float]:
    if len(xs) < 2:
        return 0.0 if xs else None
    return statistics.stdev(xs)


def _fmt(v: Optional[float]) -> str:
    if v is None:
        return ""
    return f"{v:.6g}"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scores", required=True, help="fold_scores.tsv from FOLD_MERGE_SCORES")
    p.add_argument("--pairs", required=True, help="pairs.tsv with id,target,binder")
    p.add_argument("--scores-out", required=True)
    p.add_argument("--summary-out", required=True)
    args = p.parse_args()

    pairs: Dict[str, Tuple[str, str]] = {}
    with open(args.pairs) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pairs[row["id"]] = (row["target"], row["binder"])

    score_rows: List[dict] = []
    with open(args.scores) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            tid, bid = pairs.get(row.get("id", ""), ("", ""))
            # Fallback: split id on _and_ if pairs miss a row
            if not tid and "_and_" in row.get("id", ""):
                parts = row["id"].split("_and_", 1)
                tid, bid = parts[0], parts[1]
            row = dict(row)
            row["target"] = tid
            row["binder"] = bid
            score_rows.append(row)

    out_cols = ["target", "binder"] + [c for c in fieldnames if c not in ("target", "binder")]
    # Keep canonical order: target, binder after id if present
    if "id" in out_cols:
        out_cols = ["id", "target", "binder"] + [c for c in out_cols if c not in ("id", "target", "binder")]

    with open(args.scores_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in score_rows:
            w.writerow(row)

    # Aggregate per (target, binder, tool)
    groups: Dict[Tuple[str, str, str], List[dict]] = defaultdict(list)
    for row in score_rows:
        key = (row.get("target", ""), row.get("binder", ""), row.get("tool", ""))
        groups[key].append(row)

    summary_rows = []
    # First pass: raw aggregates
    for (target, binder, tool), rows in sorted(groups.items()):
        iptms = [v for v in (_f(r.get("iptm", "")) for r in rows) if v is not None]
        ipsaes = [v for v in (_f(r.get("ipsae", "")) for r in rows) if v is not None]
        summary_rows.append({
            "target": target,
            "binder": binder,
            "tool": tool,
            "n": len(rows),
            "iptm_mean": _mean(iptms),
            "iptm_median": _median(iptms),
            "iptm_max": _max(iptms),
            "iptm_sd": _sd(iptms),
            "ipsae_mean": _mean(ipsaes),
            "ipsae_median": _median(ipsaes),
            "ipsae_max": _max(ipsaes),
            "ipsae_sd": _sd(ipsaes),
        })

    # Within-tool z-scores on iptm_mean and ipsae_mean
    by_tool: Dict[str, List[dict]] = defaultdict(list)
    for r in summary_rows:
        by_tool[r["tool"]].append(r)

    def add_z(rows: List[dict], src: str, dst: str) -> None:
        vals = [r[src] for r in rows if r[src] is not None]
        if len(vals) < 2:
            mu, sd = (vals[0] if vals else 0.0), 0.0
        else:
            mu = statistics.mean(vals)
            sd = statistics.stdev(vals)
        for r in rows:
            v = r[src]
            if v is None or sd == 0.0:
                r[dst] = 0.0 if v is not None else None
            else:
                r[dst] = (v - mu) / sd

    for tool, rows in by_tool.items():
        add_z(rows, "iptm_mean", "iptm_z")
        add_z(rows, "ipsae_mean", "ipsae_z")

    # consensus_z = mean of available per-tool z (prefer iptm_z, fall back ipsae_z)
    # Computed per (target, binder) across tools
    pair_tools: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in summary_rows:
        pair_tools[(r["target"], r["binder"])].append(r)

    consensus: Dict[Tuple[str, str], float] = {}
    for pair, rows in pair_tools.items():
        zs = []
        for r in rows:
            z = r.get("iptm_z")
            if z is None:
                z = r.get("ipsae_z")
            if z is not None:
                zs.append(z)
        consensus[pair] = statistics.mean(zs) if zs else None

    for r in summary_rows:
        r["consensus_z"] = consensus.get((r["target"], r["binder"]))

    sum_cols = [
        "target", "binder", "tool", "n",
        "iptm_mean", "iptm_median", "iptm_max", "iptm_sd", "iptm_z",
        "ipsae_mean", "ipsae_median", "ipsae_max", "ipsae_sd", "ipsae_z",
        "consensus_z",
    ]
    with open(args.summary_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sum_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in summary_rows:
            out = {}
            for c in sum_cols:
                if c == "n":
                    out[c] = r["n"]
                elif c in ("target", "binder", "tool"):
                    out[c] = r.get(c, "")
                else:
                    out[c] = _fmt(r.get(c))
            w.writerow(out)

    return 0


if __name__ == "__main__":
    sys.exit(main())
