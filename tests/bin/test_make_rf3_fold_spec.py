#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Unit tests for bin/fold/make_rf3_fold_spec.py.

    uv run --with pytest pytest tests/bin/test_make_rf3_fold_spec.py -q
"""

import importlib.util
from pathlib import Path

import pytest

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "fold" / "make_rf3_fold_spec.py"
_spec = importlib.util.spec_from_file_location("make_rf3_fold_spec", _MODPATH)
mrfs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mrfs)

TARGET = "AFTVTVPKDLYVVEY"
BINDER = "SEELRKKKFESTVLA"


def _write(p: Path, text: str) -> Path:
    p.write_text(text)
    return p


def test_chain_ordered_msas_accepted(tmp_path):
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    t = _write(tmp_path / "t.a3m", f">t TaxID=1\n{TARGET}\n")
    b = _write(tmp_path / "b.a3m", f">b TaxID=2\n{BINDER}\n")
    spec = mrfs.make_rf3_fold_spec(fasta, "x", a3m_paths=[t, b])
    assert [c["msa_path"] for c in spec["components"]] == [t.name, b.name]


def test_swapped_msas_rejected(tmp_path):
    # Bug class of commit 1ae37a8: a mis-ordered a3m bundle otherwise folds
    # each chain against another chain's MSA with no error from RF3.
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    t = _write(tmp_path / "t.a3m", f">t TaxID=1\n{TARGET}\n")
    b = _write(tmp_path / "b.a3m", f">b TaxID=2\n{BINDER}\n")
    with pytest.raises(ValueError, match="chain order"):
        mrfs.make_rf3_fold_spec(fasta, "x", a3m_paths=[b, t])


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))


def _pdl1_templates(tmp_path: Path, query: str) -> Path:
    mt_spec = importlib.util.spec_from_file_location("match_templates", _MODPATH.parent / "match_templates.py")
    mt = importlib.util.module_from_spec(mt_spec)
    mt_spec.loader.exec_module(mt)
    pdb = _MODPATH.parents[2] / "examples" / "fold-pulldown" / "input" / "PDL1.pdb"
    accepted, report = mt.match_templates([pdb], [("q", query)], 0.3, 0.3, 40, 4)
    out = tmp_path / "fold_templates"
    mt.write_outputs(accepted, report, out)
    return out


def _pdl1_seq() -> str:
    lines = (_MODPATH.parents[2] / "examples" / "fold-pulldown" / "input" / "targets.fasta").read_text().split(">")[1]
    return "".join(lines.splitlines()[1:])


def test_template_chain_becomes_inplace_structure(tmp_path):
    gemmi = pytest.importorskip("gemmi")
    pytest.importorskip("Bio")
    seq = _pdl1_seq()
    query = "GSHMA" + seq[:40] + "WWW" + seq[40:]
    fasta = _write(tmp_path / "pair.fasta", f">t\n{query}\n>b\n{BINDER}\n")
    t_a3m = _write(tmp_path / "t.a3m", f">q\n{query}\n")
    b_a3m = _write(tmp_path / "b.a3m", f">q\n{BINDER}\n")
    out = tmp_path / "out"
    out.mkdir()
    spec = mrfs.make_rf3_fold_spec(
        fasta, "x", [t_a3m, b_a3m], out_dir=out,
        templates_dir=_pdl1_templates(tmp_path, query), template_chains=["A"],
    )
    target, binder = spec["components"]
    assert target == {"path": f"{mrfs.TEMPLATE_DIR_NAME}/chain_A.cif", "custom_parse_kwargs": {"add_missing_atoms": True}}
    assert binder == {"seq": BINDER, "chain_id": "B", "msa_path": "b.a3m"}
    assert spec["template_selection"] == ["A"]
    assert spec["msa_paths"] == {"A": "t.a3m"}

    st = gemmi.read_structure(str(out / target["path"]))
    st.setup_entities()
    assert gemmi.one_letter_code(st.entities[0].full_sequence) == query
    assert [ch.name for ch in st[0]] == ["A"]
    # The 8 residues the template lacks have no atoms; the rest do.
    assert len(st[0]["A"]) == len(query) - 8


def test_no_templates_unchanged(tmp_path):
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    spec = mrfs.make_rf3_fold_spec(fasta, "x", out_dir=tmp_path)
    assert "template_selection" not in spec and "msa_paths" not in spec
    assert [c["chain_id"] for c in spec["components"]] == ["A", "B"]
