#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Unit tests for bin/fold_pulldown_summarise.py - the target x binder summary and the
ranking it implies.

Run from the repo root (host python has no pytest):

    uv run --with pytest pytest tests/bin/test_fold_pulldown_summarise.py -q

The assertions lock the three decisions that decide an ordering: which metric
consensus_z is built from, whether the standardised statistic is the max or the mean
over samples, and whether the standardisation pool is per target or global. A
regression in any of them silently reorders a shortlist, which is the kind of defect
that only shows up in the wet lab.
"""

import csv
import io
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "fold_pulldown_summarise.py"

SCORE_COLS = ["id", "tool", "model", "iptm", "ipsae"]


def _write_tsv(path, cols, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def _run(tmp_path, score_rows, pair_rows, *extra):
    scores = tmp_path / "fold_scores.tsv"
    pairs = tmp_path / "pairs.tsv"
    _write_tsv(scores, SCORE_COLS, score_rows)
    _write_tsv(pairs, ["id", "target", "binder"], pair_rows)
    out_scores = tmp_path / "out_scores.tsv"
    out_summary = tmp_path / "out_summary.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--scores", str(scores), "--pairs", str(pairs),
         "--scores-out", str(out_scores), "--summary-out", str(out_summary), *extra],
        capture_output=True, text=True, check=True,
    )
    summary = list(csv.DictReader(io.StringIO(out_summary.read_text()), delimiter="\t"))
    return summary, proc.stderr


def _by_pair(summary, tool="boltz"):
    return {(r["target"], r["binder"]): r for r in summary if r["tool"] == tool}


# --- fixtures ------------------------------------------------------------------

def _two_sample_rows():
    """One target, two binders, two samples each.

    A's samples are (iptm 0.5, 0.9) and B's are (0.7, 0.7). The means tie at 0.7,
    so mean-standardisation cannot separate them, while the max ranks A above B.
    ipsae is deliberately ordered the opposite way to iptm, so a test that passes
    under either metric is not testing the metric choice.
    """
    return (
        [
            {"id": "T_and_A", "tool": "boltz", "model": "0", "iptm": "0.5", "ipsae": "0.20"},
            {"id": "T_and_A", "tool": "boltz", "model": "1", "iptm": "0.9", "ipsae": "0.30"},
            {"id": "T_and_B", "tool": "boltz", "model": "0", "iptm": "0.7", "ipsae": "0.60"},
            {"id": "T_and_B", "tool": "boltz", "model": "1", "iptm": "0.7", "ipsae": "0.80"},
        ],
        [
            {"id": "T_and_A", "target": "T", "binder": "A"},
            {"id": "T_and_B", "target": "T", "binder": "B"},
        ],
    )


def _two_target_rows():
    """Two targets of very different difficulty, two binders each.

    Target HARD scores 0.10/0.20; target EASY scores 0.80/0.90. Within each target
    binder Y beats binder X. Pooled globally, both EASY binders outrank both HARD
    binders and the within-target ordering is swamped by target difficulty.
    """
    rows, pairs = [], []
    for target, (x, y) in {"HARD": ("0.10", "0.20"), "EASY": ("0.80", "0.90")}.items():
        for binder, val in (("X", x), ("Y", y)):
            cid = f"{target}_and_{binder}"
            rows.append({"id": cid, "tool": "boltz", "model": "0", "iptm": val, "ipsae": val})
            pairs.append({"id": cid, "target": target, "binder": binder})
    return rows, pairs


# --- consensus metric ----------------------------------------------------------

def test_consensus_defaults_to_ipsae(tmp_path):
    """ipSAE is the metric with measured discrimination, so it is the default."""
    summary, _ = _run(tmp_path, *_two_sample_rows())
    pairs = _by_pair(summary)
    assert pairs[("T", "A")]["z_basis"] == "ipsae_max/target"
    # ipsae favours B (0.80 vs 0.30), iptm favours A (0.9 vs 0.7).
    assert float(pairs[("T", "B")]["consensus_z"]) > float(pairs[("T", "A")]["consensus_z"])


def test_consensus_metric_iptm_restores_old_ordering(tmp_path):
    summary, _ = _run(tmp_path, *_two_sample_rows(), "--consensus-metric", "iptm")
    pairs = _by_pair(summary)
    assert pairs[("T", "A")]["z_basis"] == "iptm_max/target"
    assert float(pairs[("T", "A")]["consensus_z"]) > float(pairs[("T", "B")]["consensus_z"])


def test_both_z_columns_are_always_emitted(tmp_path):
    """Choosing a metric for consensus_z must not stop the other being reported."""
    summary, _ = _run(tmp_path, *_two_sample_rows())
    for row in summary:
        assert row["iptm_z"] != ""
        assert row["ipsae_z"] != ""


# --- statistic over samples ----------------------------------------------------

def test_z_stat_max_separates_where_mean_ties(tmp_path):
    rows, pairs_in = _two_sample_rows()
    summary, _ = _run(tmp_path, rows, pairs_in, "--consensus-metric", "iptm")
    pairs = _by_pair(summary)
    # means tie at 0.7; maxima are 0.9 vs 0.7
    assert float(pairs[("T", "A")]["iptm_mean"]) == pytest.approx(0.7)
    assert float(pairs[("T", "B")]["iptm_mean"]) == pytest.approx(0.7)
    assert float(pairs[("T", "A")]["iptm_z"]) > float(pairs[("T", "B")]["iptm_z"])


def test_z_stat_mean_is_degenerate_on_a_tie(tmp_path):
    rows, pairs_in = _two_sample_rows()
    summary, _ = _run(tmp_path, rows, pairs_in, "--consensus-metric", "iptm", "--z-stat", "mean")
    pairs = _by_pair(summary)
    assert pairs[("T", "A")]["z_basis"] == "iptm_mean/target"
    # sd of a tied pool is 0, so add_z assigns 0.0 to every row
    assert float(pairs[("T", "A")]["iptm_z"]) == pytest.approx(0.0)
    assert float(pairs[("T", "B")]["iptm_z"]) == pytest.approx(0.0)


# --- standardisation pool ------------------------------------------------------

def test_z_scope_target_preserves_within_target_ordering(tmp_path):
    summary, _ = _run(tmp_path, *_two_target_rows())
    pairs = _by_pair(summary)
    for target in ("HARD", "EASY"):
        assert float(pairs[(target, "Y")]["consensus_z"]) > float(pairs[(target, "X")]["consensus_z"])
    # A hard target's best binder is not penalised for its target's difficulty.
    assert float(pairs[("HARD", "Y")]["consensus_z"]) == pytest.approx(
        float(pairs[("EASY", "Y")]["consensus_z"])
    )
    assert all(r["n_pool"] == "2" for r in summary)


def test_z_scope_global_lets_target_difficulty_dominate(tmp_path):
    summary, _ = _run(tmp_path, *_two_target_rows(), "--z-scope", "global")
    pairs = _by_pair(summary)
    assert all(r["n_pool"] == "4" for r in summary)
    # The hard target's best binder now ranks below the easy target's worst.
    assert float(pairs[("HARD", "Y")]["consensus_z"]) < float(pairs[("EASY", "X")]["consensus_z"])


# --- small-pool guard ----------------------------------------------------------

def _pool_of(n, target="T", tool="boltz"):
    """n complexes against one target, all distinct, as (score rows, pair rows)."""
    rows, pairs_in = [], []
    for i in range(n):
        cid = f"{target}_and_B{i}"
        rows.append({"id": cid, "tool": tool, "model": "0",
                     "iptm": f"{0.5 + i / 100:.2f}", "ipsae": f"{0.4 + i / 100:.2f}"})
        pairs_in.append({"id": cid, "target": target, "binder": f"B{i}"})
    return rows, pairs_in


def test_small_pool_warns_with_the_reachable_bound(tmp_path):
    summary, stderr = _run(tmp_path, *_two_sample_rows())
    assert "holds 2 complex(es)" in stderr
    # (k-1)/sqrt(k) for k=2
    assert "0.707" in stderr
    # The warning has to survive into the table: a Nextflow task's stderr goes to
    # the work directory, so anything reading the summary would never see it.
    assert all(r["z_pool_small"] == "True" for r in summary)


def test_no_warning_once_the_pool_is_large_enough(tmp_path):
    rows, pairs_in = _pool_of(12)
    summary, stderr = _run(tmp_path, rows, pairs_in)
    assert "Warning" not in stderr
    assert all(r["n_pool"] == "12" for r in summary)
    assert all(r["z_pool_small"] == "False" for r in summary)


def test_z_pool_small_tracks_the_min_pool_threshold(tmp_path):
    """The column reflects --min-pool, not a hardcoded size."""
    rows, pairs_in = _pool_of(12)
    summary, stderr = _run(tmp_path, rows, pairs_in, "--min-pool", "20")
    assert "holds 12 complex(es)" in stderr
    assert all(r["z_pool_small"] == "True" for r in summary)

    summary, stderr = _run(tmp_path, *_two_sample_rows(), "--min-pool", "2")
    assert "Warning" not in stderr
    assert all(r["z_pool_small"] == "False" for r in summary)


def test_z_pool_small_is_per_pool_not_per_run(tmp_path):
    """A small target and a large one in the same run must be flagged differently."""
    big_rows, big_pairs = _pool_of(12, target="BIG")
    small_rows, small_pairs = _pool_of(2, target="SMALL")
    summary, stderr = _run(tmp_path, big_rows + small_rows, big_pairs + small_pairs)
    flags = {r["target"]: r["z_pool_small"] for r in summary}
    assert flags == {"BIG": "False", "SMALL": "True"}
    assert "SMALL/boltz" in stderr
    assert "BIG/boltz" not in stderr


# --- robustness ----------------------------------------------------------------

def test_missing_metric_falls_back_to_the_other(tmp_path):
    """A tool that emits no ipsae must still reach consensus_z, via iptm."""
    rows = [
        {"id": "T_and_A", "tool": "af2", "model": "0", "iptm": "0.9", "ipsae": ""},
        {"id": "T_and_B", "tool": "af2", "model": "0", "iptm": "0.4", "ipsae": ""},
    ]
    pairs_in = [
        {"id": "T_and_A", "target": "T", "binder": "A"},
        {"id": "T_and_B", "target": "T", "binder": "B"},
    ]
    summary, _ = _run(tmp_path, rows, pairs_in)
    pairs = _by_pair(summary, tool="af2")
    assert pairs[("T", "A")]["ipsae_z"] == ""
    assert float(pairs[("T", "A")]["consensus_z"]) > float(pairs[("T", "B")]["consensus_z"])
