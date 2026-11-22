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


def _templates_dir(tmp_path: Path, seq: str) -> Path:
    import hashlib
    import json

    d = tmp_path / "fold_templates"
    d.mkdir()
    (d / "tmpl_001.cif").write_text("data_tmpl_001\n")
    hit = {"cif": "tmpl_001.cif", "template_chain": "A", "query_indices": [0], "template_indices": [0]}
    (d / "index.json").write_text(json.dumps({hashlib.md5(seq.encode()).hexdigest(): [hit]}))
    return d


def test_templates_pinned_to_matching_chain(tmp_path):
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    d = _templates_dir(tmp_path, TARGET)
    data = mbcy.make_boltz_complex_yaml(fasta, templates_dir=str(d))
    assert data["templates"] == [{"cif": str(d / "tmpl_001.cif"), "chain_id": "A", "template_id": "A"}]


def test_template_chains_and_force(tmp_path):
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{TARGET}\n")
    d = _templates_dir(tmp_path, TARGET)
    data = mbcy.make_boltz_complex_yaml(
        fasta, templates_dir=str(d), template_chains=["A"], template_force=True, template_threshold=1.0,
    )
    assert [t["chain_id"] for t in data["templates"]] == ["A"]
    assert data["templates"][0]["force"] is True and data["templates"][0]["threshold"] == 1.0


def test_no_index_no_templates(tmp_path):
    fasta = _write(tmp_path / "mono.fasta", f">t\n{TARGET}\n")
    placeholder = _write(tmp_path / "empty_templates", "placeholder")
    data = mbcy.make_boltz_complex_yaml(fasta, templates_dir=str(placeholder))
    assert "templates" not in data


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


def test_homodimer_shares_first_chain_msa(tmp_path):
    # Boltz groups by exact sequence and rejects different msa paths, even when
    # the two files are the same search under different names.
    fasta = _write(tmp_path / "homo.fasta", f">a\n{TARGET}\n>b\n{TARGET}\n")
    a = _csv(tmp_path / "a.csv", TARGET)
    b = _csv(tmp_path / "b.csv", TARGET)
    data = mbcy.make_boltz_complex_yaml(fasta, msa_paths=[a, b])
    assert [e["protein"]["id"] for e in data["sequences"]] == [["A"], ["B"]]
    assert [e["protein"]["msa"] for e in data["sequences"]] == [a.name, a.name]


def test_homodimer_query_only_copy_keeps_the_real_msa(tmp_path):
    fasta = _write(tmp_path / "homo.fasta", f">a\n{TARGET}\n>b\n{TARGET}\n")
    a = _csv(tmp_path / "a.csv", TARGET)
    b = _csv(tmp_path / "b.csv", TARGET)
    data = mbcy.make_boltz_complex_yaml(fasta, msa_paths=[a, b], query_only_chains=["B"])
    assert [e["protein"]["msa"] for e in data["sequences"]] == [a.name, a.name]


def test_homodimer_msa_server_with_query_only_copy_omits_both(tmp_path):
    fasta = _write(tmp_path / "homo.fasta", f">a\n{TARGET}\n>b\n{TARGET}\n")
    data = mbcy.make_boltz_complex_yaml(fasta, use_msa_server=True, query_only_chains=["B"])
    for e in data["sequences"]:
        assert "msa" not in e["protein"]


def test_homodimer_both_query_only_stays_empty(tmp_path):
    fasta = _write(tmp_path / "homo.fasta", f">a\n{TARGET}\n>b\n{TARGET}\n")
    data = mbcy.make_boltz_complex_yaml(fasta, use_msa_server=True, query_only_chains=["A", "B"])
    assert [e["protein"]["msa"] for e in data["sequences"]] == ["empty", "empty"]


def test_only_identical_chains_share_an_msa(tmp_path):
    fasta = _write(tmp_path / "tri.fasta", f">a\n{TARGET}\n>b\n{TARGET}\n>c\n{BINDER}\n")
    a = _csv(tmp_path / "a.csv", TARGET)
    b = _csv(tmp_path / "b.csv", TARGET)
    c = _csv(tmp_path / "c.csv", BINDER)
    data = mbcy.make_boltz_complex_yaml(fasta, msa_paths=[a, b, c])
    assert [e["protein"]["msa"] for e in data["sequences"]] == [a.name, a.name, c.name]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
