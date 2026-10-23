#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "biopython"]
# ///

"""
Tests for bin/refold/prepare_oracle_inputs.py and bin/refold/oracle_summary_to_wide.py.

    uv run --with pytest --with biopython pytest tests/bin/test_refold_oracles.py -q
"""

import csv
import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("Bio")

_BIN = Path(__file__).resolve().parents[2] / "bin" / "refold"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare = _load("prepare_oracle_inputs")
to_wide = _load("oracle_summary_to_wide")


# Target chain A (2 residues), binder chain B (2 residues)
def _pdb(binder_seq: str) -> str:
    binder_3 = {"A": "ALA", "G": "GLY", "S": "SER"}
    lines = [
        "ATOM      1  CA  ALA A  18      32.394 -38.055 -22.849  1.00 37.28           C",
        "ATOM      2  CA  GLY A  19      31.602 -39.292 -22.650  1.00 35.59           C",
    ]
    for i, residue in enumerate(binder_seq):
        lines.append(
            f"ATOM  {3 + i:5d}  CA  {binder_3[residue]} B{1 + i:4d}"
            f"      12.000  -8.000  -{2 + i}.000  1.00 30.00           C"
        )
    return "\n".join(lines) + "\nEND\n"


@pytest.fixture()
def designs(tmp_path: Path) -> Path:
    d = tmp_path / "designs"
    d.mkdir()
    (d / "rfd4_x_0.pdb").write_text(_pdb("AG"))
    (d / "rfd4_x_1.pdb").write_text(_pdb("SA"))
    return d


def test_binder_record_per_design_named_after_the_file(designs: Path, tmp_path: Path):
    sequences = prepare.collect_sequences(sorted(designs.glob("*.pdb")), "B", "A")
    assert set(sequences) == {"rfd4_x_0", "rfd4_x_1"}
    assert sequences["rfd4_x_0"]["binder"] == "AG"
    assert sequences["rfd4_x_1"]["binder"] == "SA"


def test_target_is_shared_across_designs(designs: Path):
    sequences = prepare.collect_sequences(sorted(designs.glob("*.pdb")), "B", "A")
    assert {v["target"] for v in sequences.values()} == {"AG"}
    assert prepare.consensus_target(sequences) == "AG"


def test_consensus_target_takes_the_majority_and_warns(caplog):
    sequences = {
        "a": {"binder": "A", "target": "AG"},
        "b": {"binder": "A", "target": "AG"},
        "c": {"binder": "A", "target": "AGS"},
    }
    assert prepare.consensus_target(sequences) == "AG"
    assert "disagree on the target sequence" in caplog.text


def test_design_id_strips_compound_suffixes():
    assert prepare._design_id(Path("x/rfd4_a_0.cif.gz")) == "rfd4_a_0"
    assert prepare._design_id(Path("x/rfd4_a_0.cif")) == "rfd4_a_0"
    assert prepare._design_id(Path("x/rfd4_a_0.pdb")) == "rfd4_a_0"


def test_a_design_missing_the_binder_chain_is_skipped(tmp_path: Path, caplog):
    only_target = tmp_path / "no_binder.pdb"
    only_target.write_text(
        "ATOM      1  CA  ALA A  18      32.394 -38.055 -22.849  1.00 37.28           C\nEND\n"
    )
    assert prepare.collect_sequences([only_target], "B", "A") == {}
    assert "no chain B" in caplog.text


SUMMARY_ROWS = [
    {
        "target": "target", "binder": "rfd4_x_0", "tool": "boltz", "n": "2",
        "iptm_max": "0.78", "ipsae_max": "0.55",
        "consensus_z": "1.42", "consensus_z_metric": "ipsae", "n_tools": "2",
    },
    {
        "target": "target", "binder": "rfd4_x_0", "tool": "chai1", "n": "2",
        "iptm_max": "0.81", "ipsae_max": "0.62",
        "consensus_z": "1.42", "consensus_z_metric": "ipsae", "n_tools": "2",
    },
    {
        "target": "target", "binder": "rfd4_x_1", "tool": "boltz", "n": "2",
        "iptm_max": "0.40", "ipsae_max": "0.11",
        "consensus_z": "-0.9", "consensus_z_metric": "ipsae", "n_tools": "1",
    },
]


def test_pivot_gives_one_row_per_design_with_tool_prefixed_columns():
    rows, columns = to_wide.pivot(SUMMARY_ROWS, ["iptm_max", "ipsae_max"])
    by_id = {r["id"]: r for r in rows}
    assert set(by_id) == {"rfd4_x_0", "rfd4_x_1"}
    assert by_id["rfd4_x_0"]["chai1_iptm_max"] == "0.81"
    assert by_id["rfd4_x_0"]["boltz_ipsae_max"] == "0.55"
    # A design scored by only one tool leaves the other tool's columns blank rather
    # than carrying the first tool's value across.
    assert by_id["rfd4_x_1"]["chai1_iptm_max"] == ""
    assert columns[0] == "id"


def test_pair_level_columns_are_written_once_not_per_tool():
    rows, columns = to_wide.pivot(SUMMARY_ROWS, ["iptm_max"])
    assert "consensus_z" in columns
    assert not any(c.startswith("boltz_consensus_z") for c in columns)
    assert {r["id"]: r["consensus_z"] for r in rows}["rfd4_x_0"] == "1.42"


def test_several_targets_keep_their_columns_apart():
    rows = SUMMARY_ROWS + [
        {
            "target": "target2", "binder": "rfd4_x_0", "tool": "boltz", "n": "2",
            "iptm_max": "0.20", "ipsae_max": "0.05",
            "consensus_z": "0.1", "consensus_z_metric": "ipsae", "n_tools": "1",
        }
    ]
    pivoted, columns = to_wide.pivot(rows, ["iptm_max"])
    assert "target_boltz_iptm_max" in columns
    assert "target2_boltz_iptm_max" in columns
    row = {r["id"]: r for r in pivoted}["rfd4_x_0"]
    assert row["target_boltz_iptm_max"] == "0.78"
    assert row["target2_boltz_iptm_max"] == "0.20"


def test_empty_summary_still_writes_a_header(tmp_path: Path, monkeypatch):
    summary = tmp_path / "summary.tsv"
    summary.write_text("target\tbinder\ttool\n")
    out = tmp_path / "out.tsv"
    monkeypatch.setattr("sys.argv", ["x", str(summary), "-o", str(out)])
    assert to_wide.main() == 0
    assert out.read_text() == "id\n"


def test_end_to_end_writes_both_fastas(designs: Path, tmp_path: Path, monkeypatch):
    targets = tmp_path / "targets.fasta"
    binders = tmp_path / "binders.fasta"
    monkeypatch.setattr(
        "sys.argv",
        [
            "x", str(designs),
            "--binder-chain", "B", "--target-chain", "A",
            "--targets-out", str(targets), "--binders-out", str(binders),
        ],
    )
    assert prepare.main() == 0
    assert targets.read_text() == ">target\nAG\n"
    assert binders.read_text() == ">rfd4_x_0\nAG\n>rfd4_x_1\nSA\n"


def test_wide_output_is_valid_tsv(tmp_path: Path, monkeypatch):
    summary = tmp_path / "summary.tsv"
    with open(summary, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(SUMMARY_ROWS[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(SUMMARY_ROWS)
    out = tmp_path / "wide.tsv"
    monkeypatch.setattr("sys.argv", ["x", str(summary), "-o", str(out)])
    assert to_wide.main() == 0
    with open(out, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["id"] for r in rows] == ["rfd4_x_0", "rfd4_x_1"]
