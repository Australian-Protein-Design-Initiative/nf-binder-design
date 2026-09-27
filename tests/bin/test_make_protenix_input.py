#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Unit tests for bin/fold/make_protenix_input.py.

    uv run --with pytest pytest tests/bin/test_make_protenix_input.py -q
"""

import importlib.util
from pathlib import Path

import pytest

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "fold" / "make_protenix_input.py"
_spec = importlib.util.spec_from_file_location("make_protenix_input", _MODPATH)
mpi = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mpi)

TARGET = "AFTVTVPKDLYVVEY"
BINDER = "SEELRKKKFESTVLA"


def _write(p: Path, text: str) -> Path:
    p.write_text(text)
    return p


def test_chain_ordered_msas_accepted(tmp_path):
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    t = _write(tmp_path / "PDL1.protenix_unpaired.a3m", f"#15\t1\n>t\n{TARGET}\n")
    b = _write(tmp_path / "AAA.protenix_unpaired.a3m", f">b\n{BINDER}\n")
    spec = mpi.make_protenix_input(fasta, "x", unpaired_a3m_paths=[t, b])
    assert [s["proteinChain"]["unpairedMsaPath"] for s in spec["sequences"]] == [t.name, b.name]


def test_swapped_msas_rejected(tmp_path):
    # Regression: pulldown bundles were sorted by file name, so a binder id that
    # sorts before the target id gave each chain the other chain's MSA.
    fasta = _write(tmp_path / "pair.fasta", f">t\n{TARGET}\n>b\n{BINDER}\n")
    t = _write(tmp_path / "PDL1.protenix_paired.a3m", f">t\n{TARGET}\n")
    b = _write(tmp_path / "AAA.protenix_paired.a3m", f">b\n{BINDER}\n")
    with pytest.raises(ValueError, match="chain order"):
        mpi.make_protenix_input(fasta, "x", paired_a3m_paths=[b, t])


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
