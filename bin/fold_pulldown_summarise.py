#!/usr/bin/env python3
# /// script
# requires-python = ">=3.8"
# ///
"""
Join fold_scores.tsv with pairs.tsv and emit:

  fold_pulldown_scores.tsv  - one row per predicted structure (+ target, binder)
  fold_pulldown_summary.tsv - one row per (target, binder, tool) with aggregate
                              iptm/ipsae stats, per-target within-tool z-scores,
                              and cross-tool consensus_z.

Three choices govern the ranking, all overridable:

  --z-scope target   Standardise within (target, tool). Raw co-folding scores are
                     not comparable across targets, and target difficulty varies
                     more than design quality does within a target, so pooling
                     targets makes a complex rank partly on which target it was
                     paired with. Use `global` for the old single-pool behaviour.
  --z-stat max       Standardise the max over samples. One diffusion sample is
                     noisy; taking the max is why --n_predictions > 1 is worth
                     paying for. Use `mean` for the old behaviour.
  --consensus-metric ipsae
                     Build consensus_z from ipSAE rather than ipTM. ipSAE
                     (Dunbrack 2025) isolates the interface from whole-complex
                     confidence, which ipTM mixes together. Use `iptm` for the
                     old behaviour.

Each row records the basis used (`z_basis`), the size of the pool its z-score was
computed over (`n_pool`), and whether that pool was smaller than --min-pool
(`z_pool_small`), because a z-score over a handful of complexes is a rank label
rather than a distance: with k complexes the largest possible absolute z is
(k-1)/sqrt(k).
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
    # A single replicate has no spread to report; 0.0 previously implied
    # "measured and tight", indistinguishable from a real zero-variance pool.
    if len(xs) < 2:
        return None
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
    p.add_argument(
        "--consensus-metric", choices=("ipsae", "iptm"), default="ipsae",
        help="Metric whose per-tool z-scores are averaged into consensus_z [default: ipsae]",
    )
    p.add_argument(
        "--z-stat", choices=("max", "mean"), default="max",
        help="Per-complex statistic over samples that is standardised [default: max]",
    )
    p.add_argument(
        "--z-scope", choices=("target", "global"), default="target",
        help="Pool for standardisation: within (target, tool), or all targets together "
             "[default: target]",
    )
    p.add_argument(
        "--min-pool", type=int, default=10,
        help="Pools holding fewer than this many complexes are flagged z_pool_small "
             "in the summary and warned about on stderr [default: 10]",
    )
    args = p.parse_args()

    pairs: Dict[str, Tuple[str, str]] = {}
    with open(args.pairs) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pairs[row["id"]] = (row["target"], row["binder"])

    score_rows: List[dict] = []
    missing_ids: List[str] = []
    with open(args.scores) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            row_id = row.get("id", "")
            if row_id not in pairs:
                # Guessing target/binder by splitting on '_and_' is ambiguous
                # (see FoldIds.groovy pair-id collision check) and silently
                # mis-attributes a complex to the wrong target/binder - fail
                # fast instead so a missing pairs.tsv row is never invented.
                missing_ids.append(row_id)
                continue
            tid, bid = pairs[row_id]
            row = dict(row)
            row["target"] = tid
            row["binder"] = bid
            score_rows.append(row)

    if missing_ids:
        print(
            f"ERROR: {len(missing_ids)} score row id(s) missing from {args.pairs}: "
            f"{sorted(set(missing_ids))}",
            file=sys.stderr,
        )
        return 1

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

    # Standardise within each pool. The pool is (target, tool) by default: raw
    # co-folding scores are not comparable across targets, so pooling them makes a
    # complex rank partly on its target's difficulty rather than on the design.
    pools: Dict[Tuple[str, ...], List[dict]] = defaultdict(list)
    for r in summary_rows:
        key = (r["tool"],) if args.z_scope == "global" else (r["target"], r["tool"])
        pools[key].append(r)

    def add_z(rows: List[dict], src: str, dst: str) -> int:
        """Standardise rows[*][src] into rows[*][dst]; return the count of
        complexes that actually had src (mu/sd are computed over those only -
        the pool can hold more complexes than contributed a value)."""
        vals = [r[src] for r in rows if r[src] is not None]
        n_used = len(vals)
        if n_used < 2:
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
        return n_used

    iptm_src = f"iptm_{args.z_stat}"
    ipsae_src = f"ipsae_{args.z_stat}"
    z_basis = f"{args.consensus_metric}_{args.z_stat}/{args.z_scope}"

    primary_metric = args.consensus_metric
    secondary_metric = "iptm" if primary_metric == "ipsae" else "ipsae"
    primary_src = f"{primary_metric}_{args.z_stat}"

    # Whether the whole (target, tool) [or (tool,) for z_scope=global] pool
    # lacks the primary metric entirely - only then is the secondary metric
    # substituted for consensus_z, and for every complex in that pool, not
    # per-complex - mixing pools per-row would let a complex whose primary
    # metric merely failed (rather than being unsupported by the tool) beat
    # complexes that were scored on the primary metric throughout the pool.
    pool_lacks_primary: Dict[Tuple[str, ...], bool] = {}

    for key, rows in sorted(pools.items()):
        iptm_n = add_z(rows, iptm_src, "iptm_z")
        ipsae_n = add_z(rows, ipsae_src, "ipsae_z")
        primary_n = ipsae_n if primary_metric == "ipsae" else iptm_n
        pool_lacks_primary[key] = primary_n == 0
        small = primary_n < args.min_pool
        for r in rows:
            r["n_pool"] = primary_n
            r["iptm_n_pool"] = iptm_n
            r["ipsae_n_pool"] = ipsae_n
            r["z_basis"] = z_basis
            # Carried in the table as well as on stderr: a Nextflow task's stderr
            # ends up in the work directory, where nothing that reads the summary
            # will see it.
            r["z_pool_small"] = small
        if small:
            # With k complexes the largest possible |z| is (k-1)/sqrt(k), so a small
            # pool yields z-scores that carry only the ordering.
            print(
                f"Warning: z-score pool {'/'.join(key)} holds {primary_n} complex(es) "
                f"with {primary_metric}; |z| cannot exceed "
                f"{max(primary_n - 1, 0) / max(primary_n, 1) ** 0.5:.3f}. "
                "Treat these z-scores as rank labels, not distances.",
                file=sys.stderr,
            )

    # consensus_z = mean over tools of the chosen metric's z for that complex.
    # The secondary metric is substituted only where its whole pool lacked the
    # primary metric (see pool_lacks_primary above); otherwise a missing
    # primary value for one complex just leaves that complex's contribution NA.
    for r in summary_rows:
        pool_key = (r["tool"],) if args.z_scope == "global" else (r["target"], r["tool"])
        if pool_lacks_primary.get(pool_key):
            z = r.get(f"{secondary_metric}_z")
            used = secondary_metric if z is not None else None
        else:
            z = r.get(f"{primary_metric}_z")
            used = primary_metric if z is not None else None
        r["_consensus_z_contrib"] = z
        r["_consensus_z_metric"] = used

    pair_tools: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in summary_rows:
        pair_tools[(r["target"], r["binder"])].append(r)

    consensus: Dict[Tuple[str, str], Optional[float]] = {}
    n_tools: Dict[Tuple[str, str], int] = {}
    metrics_used: Dict[Tuple[str, str], str] = {}
    for pair, rows in pair_tools.items():
        contribs = [(r["tool"], r["_consensus_z_contrib"], r["_consensus_z_metric"]) for r in rows]
        zs = [z for _tool, z, _m in contribs if z is not None]
        consensus[pair] = statistics.mean(zs) if zs else None
        n_tools[pair] = len(zs)
        metrics_used[pair] = ",".join(sorted({m for _t, z, m in contribs if z is not None}))

    for r in summary_rows:
        pair = (r["target"], r["binder"])
        r["consensus_z"] = consensus.get(pair)
        r["n_tools"] = n_tools.get(pair, 0)
        r["consensus_z_metric"] = metrics_used.get(pair, "")

    sum_cols = [
        "target", "binder", "tool", "n",
        "iptm_mean", "iptm_median", "iptm_max", "iptm_sd", "iptm_z", "iptm_n_pool",
        "ipsae_mean", "ipsae_median", "ipsae_max", "ipsae_sd", "ipsae_z", "ipsae_n_pool",
        "consensus_z", "consensus_z_metric", "n_tools", "z_basis", "n_pool", "z_pool_small",
    ]
    with open(args.summary_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sum_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in summary_rows:
            out = {}
            for c in sum_cols:
                if c in ("n", "n_pool", "iptm_n_pool", "ipsae_n_pool", "n_tools"):
                    out[c] = r[c]
                elif c == "z_pool_small":
                    out[c] = "True" if r.get(c) else "False"
                elif c in ("target", "binder", "tool", "z_basis", "consensus_z_metric"):
                    out[c] = r.get(c, "")
                else:
                    out[c] = _fmt(r.get(c))
            w.writerow(out)

    return 0


if __name__ == "__main__":
    sys.exit(main())
