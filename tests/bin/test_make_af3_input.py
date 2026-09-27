#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Unit tests for bin/fold/make_af3_input.py.

    uv run --with pytest pytest tests/bin/test_make_af3_input.py -q
"""

import importlib.util
import json
from pathlib import Path

import pytest

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "fold" / "make_af3_input.py"
_spec = importlib.util.spec_from_file_location("make_af3_input", _MODPATH)
mai = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mai)

SEQ_A = "MKTAYIAKQR"
SEQ_B = "GSHMLEDPV"


def _write(p: Path, text: str) -> Path:
    p.write_text(text)
    return p


def test_monomer_query_only_when_no_msa(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">x\n{SEQ_A}\n")
    spec = mai.make_af3_input(fasta, "my job/1", tmp_path, seed=7)
    assert spec["name"] == "my_job1"
    assert spec["modelSeeds"] == [7]
    assert spec["dialect"] == "alphafold3" and spec["version"] >= 2
    prot = spec["sequences"][0]["protein"]
    assert prot["id"] == "A" and prot["sequence"] == SEQ_A
    assert prot["templates"] == []
    # AF3 with --run_data_pipeline=false needs both MSAs non-null
    for key in ("unpairedMsaPath", "pairedMsaPath"):
        text = (tmp_path / prot[key]).read_text()
        assert text.splitlines()[1] == SEQ_A


def test_colabfold_header_and_nul_stripped(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">x\n{SEQ_A}\n")
    a3m = _write(tmp_path / "x.a3m", f"#10\t1\n>101\n{SEQ_A}\n>hit\x00\nMKTAYIAKQc-\n")
    spec = mai.make_af3_input(fasta, "x", tmp_path, seed=1, unpaired_a3m_paths=[a3m])
    text = (tmp_path / spec["sequences"][0]["protein"]["unpairedMsaPath"]).read_text()
    assert text.startswith(">101\n" + SEQ_A)
    assert "#" not in text and "\x00" not in text
    assert ">hit\n" in text


def test_query_mismatch_raises(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">x\n{SEQ_A}\n")
    a3m = _write(tmp_path / "x.a3m", f">101\n{SEQ_B}\n")
    with pytest.raises(ValueError, match="query"):
        mai.make_af3_input(fasta, "x", tmp_path, seed=1, unpaired_a3m_paths=[a3m])


def test_multimer_per_chain_files(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">a\n{SEQ_A}\n>b\n{SEQ_B}\n")
    ua = _write(tmp_path / "a.a3m", f">q\n{SEQ_A}\n")
    ub = _write(tmp_path / "b.a3m", f">q\n{SEQ_B}\n")
    pa = _write(tmp_path / "a.p.a3m", f">q\n{SEQ_A}\n>tr|Q8QRZ0|Q8QRZ0_HUMAN\n{SEQ_A}\n")
    pb = _write(tmp_path / "b.p.a3m", f">q\n{SEQ_B}\n>tr|P12345|P12345_HUMAN\n{SEQ_B}\n")
    spec = mai.make_af3_input(
        fasta, "cplx", tmp_path, seed=1, unpaired_a3m_paths=[ua, ub], paired_a3m_paths=[pa, pb]
    )
    ids = [s["protein"]["id"] for s in spec["sequences"]]
    assert ids == ["A", "B"]
    pb_text = (tmp_path / spec["sequences"][1]["protein"]["pairedMsaPath"]).read_text()
    assert "P12345_HUMAN" in pb_text
    json.dumps(spec)


def test_per_chain_count_mismatch_raises(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">a\n{SEQ_A}\n>b\n{SEQ_B}\n")
    ua = _write(tmp_path / "a.a3m", f">q\n{SEQ_A}\n")
    with pytest.raises(ValueError, match="one per chain"):
        mai.make_af3_input(fasta, "x", tmp_path, seed=1, unpaired_a3m_paths=[ua])


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
