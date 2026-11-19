#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Tests for bin/pulldown/complex_structures.py's chain resolution.

    uv run --with pytest pytest tests/bin/test_complex_structures.py -q
"""

import importlib.util
from pathlib import Path

import pytest

_BIN = Path(__file__).resolve().parents[2] / "bin"


def _load():
    spec = importlib.util.spec_from_file_location(
        "complex_structures", _BIN / "pulldown" / "complex_structures.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cs = _load()

TARGET = "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR"
BINDER = "EAALLAAIAA"


def test_default_takes_the_one_non_binder_chain():
    chains = {"A": BINDER, "B": TARGET}
    assert cs.resolve_chains(chains, ["A"], None, [], "c1") == (["A"], "B")


def test_explicit_target_chain_wins():
    chains = {"A": BINDER, "B": TARGET, "C": TARGET}
    assert cs.resolve_chains(chains, ["A"], ["C"], [], "c1") == (["A"], "C")


def test_ambiguous_target_is_an_error_rather_than_a_guess():
    # Two non-binder chains and no --target_chains: picking one silently would
    # superimpose on an arbitrary chain.
    chains = {"A": BINDER, "B": TARGET, "C": TARGET}
    with pytest.raises(ValueError, match="ambiguous"):
        cs.resolve_chains(chains, ["A"], None, [], "c1")


def test_auto_picks_the_chain_matching_the_target_fasta():
    # Chain naming reversed relative to the usual convention: auto has to go by
    # sequence, not by position or by "the other one".
    chains = {"A": TARGET, "B": BINDER}
    assert cs.resolve_chains(chains, ["B"], ["auto"], [TARGET], "c1") == (["B"], "A")


def test_auto_matches_a_truncated_construct_inside_a_full_length_target():
    full = "MKWVTFISLL" + TARGET + "GGSGGS"
    chains = {"A": TARGET, "B": BINDER}
    assert cs.resolve_chains(chains, ["B"], ["auto"], [full], "c1") == (["B"], "A")


def test_auto_without_targets_is_an_error():
    chains = {"A": TARGET, "B": BINDER}
    with pytest.raises(ValueError, match="needs --targets"):
        cs.resolve_chains(chains, ["B"], ["auto"], [], "c1")


def test_missing_binder_chain_is_an_error():
    chains = {"A": BINDER, "B": TARGET}
    with pytest.raises(ValueError, match="binder chains"):
        cs.resolve_chains(chains, ["Z"], None, [], "c1")


def test_several_binder_chains_are_all_returned():
    chains = {"A": BINDER, "B": BINDER, "T": TARGET}
    assert cs.resolve_chains(chains, ["A", "B"], None, [], "c1") == (["A", "B"], "T")


def test_identity_is_over_the_shorter_sequence():
    assert cs.identity("ABC", "ABC") == 1.0
    assert cs.identity("ABC", "XABCX") == 1.0
    assert cs.identity("ABC", "AXC") == pytest.approx(2 / 3)
    assert cs.identity("", "ABC") == 0.0


def test_complex_id_strips_structure_suffixes():
    assert cs.complex_id(Path("d1.pdb")) == "d1"
    assert cs.complex_id(Path("d1.cif.gz")) == "d1"
    assert cs.complex_id(Path("a.b.pdb")) == "a.b"


# --- main(), which is where the record naming lives ---------------------------
#
# A run caught a collision these guard: the derived target record was named after the
# complex it came from, which is also that complex's binder record id, and the pulldown
# stages both chains' files side by side -- "multiple input files for each of the
# following file names" before any folding happened.

DESIGN = Path(__file__).resolve().parents[1] / "data" / "refold" / "design.pdb"


def _run(tmp_path, *args):
    pytest.importorskip("Bio")
    out = tmp_path / "binders.fasta"
    chains = tmp_path / "chains.tsv"
    argv = [
        str(DESIGN),
        "--binder-chains", "B",
        "--binders-out", str(out),
        "--chains-out", str(chains),
        *args,
    ]
    import sys

    old = sys.argv
    sys.argv = ["complex_structures.py", *argv]
    try:
        return cs.main(), out, chains
    finally:
        sys.argv = old


def test_derived_target_does_not_take_the_binder_record_id(tmp_path):
    targets = tmp_path / "targets.fasta"
    rc, binders, _ = _run(tmp_path, "--targets-out", str(targets))
    assert rc == 0
    binder_ids = {rid for rid, _ in cs.read_fasta(binders)}
    target_ids = {rid for rid, _ in cs.read_fasta(targets)}
    assert target_ids == {"target"}
    assert not (binder_ids & target_ids)


def test_supplied_target_colliding_with_a_complex_stem_is_rejected(tmp_path):
    collide = tmp_path / "targets.fasta"
    collide.write_text(f">{cs.complex_id(DESIGN)}\n{TARGET}\n")
    rc, _, _ = _run(tmp_path, "--targets", str(collide))
    assert rc == 1
