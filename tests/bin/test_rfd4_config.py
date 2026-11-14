#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""Tests for bin/rfd4_proteina/rfd4_config.py (stdlib; no RFD4-Proteina container)."""

import importlib.util
import json
from pathlib import Path

import pytest

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "rfd4_proteina" / "rfd4_config.py"
_spec = importlib.util.spec_from_file_location("rfd4_config", _MODPATH)
rfd4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rfd4)


@pytest.mark.parametrize(
    "contig,expected",
    [
        ("50-120, /0, A17-131", ("A", "B")),
        ("A17-131, /0, 50-120", ("A", "B")),
        ("50-120,/0,A17-131", ("A", "B")),
        ("C17-131, /0, 50-120", ("C", "A")),
        ("A1-10, /0, B1-10, /0, 30-40", ("A", "C")),
        # A trailing P/R/D sets the chain type of an expandable segment
        ("20-30P, /0, A5-50", ("A", "B")),
    ],
)
def test_infer_chains(contig, expected) -> None:
    assert rfd4.infer_chains(contig) == expected


def test_infer_chains_named_expandable_joins_its_run() -> None:
    # E1 shares a run with B1-10, so it lands on chain B and there is no designed chain left
    with pytest.raises(SystemExit, match="no designed chain"):
        rfd4.infer_chains("A1-10, /0, E1:1-20, B1-10")


def test_infer_chains_requires_a_preserved_segment() -> None:
    with pytest.raises(SystemExit, match="no preserved segment"):
        rfd4.infer_chains("50-120, /0, 30-40")


def test_target_chain_ids_deduplicates_in_order() -> None:
    assert rfd4.target_chain_ids("B1-5, A1-10, /0, 50-60, /0, A20-30") == ["B", "A"]


def test_parse_hotspots() -> None:
    assert rfd4.parse_hotspots("A56, A115,B7") == [
        "sel('A/*/56')",
        "sel('A/*/115')",
        "sel('B/*/7')",
    ]


def test_parse_hotspots_rejects_a_bare_number() -> None:
    with pytest.raises(SystemExit, match="Could not parse hotspot"):
        rfd4.parse_hotspots("56")


def test_generate_spec_uses_the_default_binder_recipe(tmp_path: Path) -> None:
    out = tmp_path / "spec.json"
    rfd4.generate_spec("pdl1", "/elsewhere/PDL1.cif", "50-120, /0, A17-131", "A56", str(out))
    spec = json.loads(out.read_text())

    assert list(spec) == ["pdl1"]
    design = spec["pdl1"]
    # The input is referenced by basename, so the spec works against Nextflow-staged files
    assert design["input"] == "PDL1.cif"
    assert design["contig"] == "50-120, /0, A17-131"
    assert design["conditions"]["C_CRDGRP"] == {"1": [{"select": "sel('A/*/*/[N,CA,C,O]')"}]}
    assert design["conditions"]["C_SEQ"] == {"true": [{"select": "sel('A')"}]}
    assert design["conditions"]["C_HOT"] == {"true": [{"select": "sel('A/*/56')"}]}


def test_generate_spec_omits_hotspots_when_none_given(tmp_path: Path) -> None:
    out = tmp_path / "spec.json"
    rfd4.generate_spec("d", "t.cif", "50-120, /0, A17-131", None, str(out))
    assert "C_HOT" not in json.loads(out.read_text())["d"]["conditions"]


def test_generate_spec_preserves_annotations_for_partial_diffusion(tmp_path: Path) -> None:
    out = tmp_path / "spec.json"
    rfd4.generate_spec("d", "t.cif", "50-120, /0, A17-131", None, str(out), 0.3)
    assert json.loads(out.read_text())["d"]["input_as_annotated_cif"] is True


def _write_spec(tmp_path: Path, spec: dict) -> Path:
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(spec))
    return path


def test_parse_inputs_resolves_relative_to_the_spec(tmp_path: Path) -> None:
    (tmp_path / "inputs").mkdir()
    path = _write_spec(tmp_path, {"d": {"input": "inputs/target.cif"}})
    assert rfd4.parse_inputs(str(path), None) == [str(tmp_path / "inputs" / "target.cif")]


def test_parse_inputs_covers_global_args_lists_and_duplicates(tmp_path: Path) -> None:
    path = _write_spec(
        tmp_path,
        {
            "global_args": {"input": "shared.cif"},
            "d1": {"input": ["a.cif", "b.cif"]},
            "d2": {"input": "a.cif"},
        },
    )
    names = [Path(p).name for p in rfd4.parse_inputs(str(path), None)]
    assert names == ["shared.cif", "a.cif", "b.cif"]


def test_parse_inputs_keeps_absolute_paths(tmp_path: Path) -> None:
    path = _write_spec(tmp_path, {"d": {"input": "/abs/target.cif"}})
    assert rfd4.parse_inputs(str(path), None) == ["/abs/target.cif"]


def test_stage_rewrites_every_input_to_a_basename(tmp_path: Path) -> None:
    path = _write_spec(
        tmp_path,
        {
            "global_args": {"input": "dir/shared.cif"},
            "d1": {"input": ["x/a.cif", "y/b.cif"], "contig": "50-60, /0, A1-10"},
        },
    )
    out = tmp_path / "staged.json"
    rfd4.stage_spec(str(path), str(out))
    staged = json.loads(out.read_text())

    assert staged["global_args"]["input"] == "shared.cif"
    assert staged["d1"]["input"] == ["a.cif", "b.cif"]
    # Everything else is carried through untouched
    assert staged["d1"]["contig"] == "50-60, /0, A1-10"
