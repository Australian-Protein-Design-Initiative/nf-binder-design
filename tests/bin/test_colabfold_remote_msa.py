#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Unit tests for bin/colabfold_remote_msa.py output-path handling.

    uv run --with pytest pytest tests/bin/test_colabfold_remote_msa.py -q

The API call is stubbed out, so these never touch the network.
"""

import importlib.util
from pathlib import Path
from typing import List

import pytest

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "colabfold_remote_msa.py"
_spec = importlib.util.spec_from_file_location("colabfold_remote_msa", _MODPATH)
crm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(crm)


@pytest.fixture
def no_api(monkeypatch):
    """Replace the remote search with a stub, and record whether it ran."""
    calls: List[List[str]] = []

    def fake_run_remote_msa(seqs, prefix, **kwargs):
        calls.append(list(seqs))
        return [f">query\n{s}\n" for s in seqs]

    monkeypatch.setattr(crm, "run_remote_msa", fake_run_remote_msa)
    return calls


def _fasta(tmp_path: Path, *records: str) -> Path:
    p = tmp_path / "in.fasta"
    p.write_text("".join(f">{rid}\n{seq}\n" for rid, seq in (r.split(":") for r in records)))
    return p


def _run(monkeypatch, fasta: Path, out: Path) -> int:
    monkeypatch.setattr(
        crm.sys, "argv", ["colabfold_remote_msa", "--fasta", str(fasta), "-o", str(out)]
    )
    return crm.main()


def test_single_sequence_writes_the_named_file(tmp_path, monkeypatch, no_api):
    fasta = _fasta(tmp_path, "PDL1:AAAA")
    out = tmp_path / "msa" / "PDL1.a3m"
    assert _run(monkeypatch, fasta, out) == 0
    assert out.is_file()


def test_directory_writes_one_a3m_per_record(tmp_path, monkeypatch, no_api):
    """The shape the pipeline uses: MMSEQS_COLABFOLDSEARCH passes `-o result/`."""
    fasta = _fasta(tmp_path, "PDL1:AAAA", "CD47:CCCC")
    out = tmp_path / "result"
    assert _run(monkeypatch, fasta, out) == 0
    assert sorted(p.name for p in out.iterdir()) == ["CD47.a3m", "PDL1.a3m"]


def test_file_output_with_several_sequences_is_rejected(tmp_path, monkeypatch, no_api):
    """Writing siblings beside the requested file would silently ignore the name."""
    fasta = _fasta(tmp_path, "PDL1:AAAA", "CD47:CCCC")
    out = tmp_path / "msa" / "both.a3m"
    assert _run(monkeypatch, fasta, out) == 1
    assert not out.exists()
    assert not (tmp_path / "msa").exists()
    assert no_api == [], "should fail before spending an API call"


def test_existing_directory_with_a_dot_is_not_treated_as_a_file(
    tmp_path, monkeypatch, no_api
):
    fasta = _fasta(tmp_path, "PDL1:AAAA", "CD47:CCCC")
    out = tmp_path / "v1.2"
    out.mkdir()
    assert _run(monkeypatch, fasta, out) == 0
    assert sorted(p.name for p in out.iterdir()) == ["CD47.a3m", "PDL1.a3m"]
