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
