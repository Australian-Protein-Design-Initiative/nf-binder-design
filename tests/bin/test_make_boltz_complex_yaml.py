#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "pyyaml"]
# ///

"""
Unit tests for bin/fold/make_boltz_complex_yaml.py.

    uv run --with pytest --with pyyaml pytest tests/bin/test_make_boltz_complex_yaml.py -q
"""

import importlib.util
from pathlib import Path

import pytest

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "fold" / "make_boltz_complex_yaml.py"
_spec = importlib.util.spec_from_file_location("make_boltz_complex_yaml", _MODPATH)
mbcy = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mbcy)

TARGET = "AFTVTVPKDLYVVEY"
BINDER = "SEELRKKKFESTVLA"


def _write(p: Path, text: str) -> Path:
    p.write_text(text)
    return p


def _csv(p: Path, seq: str) -> Path:
    return _write(p, f"key,sequence\n,{seq}\n")


def test_chain_ordered_msas_accepted(tmp_path):
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    t = _csv(tmp_path / "t.csv", TARGET)
    b = _csv(tmp_path / "b.csv", BINDER)
    data = mbcy.make_boltz_complex_yaml(fasta, msa_paths=[t, b])
    assert [e["protein"]["msa"] for e in data["sequences"]] == [t.name, b.name]


def test_swapped_msas_rejected(tmp_path):
    # Bug class of commit 1ae37a8: a mis-ordered CSV bundle otherwise folds
    # each chain against another chain's MSA with no error from Boltz.
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    t = _csv(tmp_path / "t.csv", TARGET)
    b = _csv(tmp_path / "b.csv", BINDER)
    with pytest.raises(ValueError, match="chain order"):
        mbcy.make_boltz_complex_yaml(fasta, msa_paths=[b, t])


def test_templates_dir_added(tmp_path):
    fasta = _write(tmp_path / "mono.fasta", f">t\n{TARGET}\n")
    templates_dir = tmp_path / "templates"
    templates_dir.mkdir()
    cif = templates_dir / "1abc.cif"
    cif.write_text("data_1abc\n")
    data = mbcy.make_boltz_complex_yaml(fasta, templates_dir=str(templates_dir))
    assert data["templates"] == [{"cif": str(cif)}]


def test_query_only_chains_force_empty_msa_under_msa_server(tmp_path):
    # meta.query_only_chains contract with fold_pulldown: --use_msa_server would
    # otherwise fetch an MSA for every chain, contradicting --create_binder_msa
    # false for a query-only (e.g. binder) chain.
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    data = mbcy.make_boltz_complex_yaml(fasta, use_msa_server=True, query_only_chains=["B"])
    entries = {e["protein"]["id"][0]: e["protein"] for e in data["sequences"]}
    assert entries["A"].get("msa") is None
    assert entries["B"]["msa"] == "empty"


def test_use_msa_server_without_query_only_chains_omits_msa(tmp_path):
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    data = mbcy.make_boltz_complex_yaml(fasta, use_msa_server=True)
    for e in data["sequences"]:
        assert "msa" not in e["protein"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
