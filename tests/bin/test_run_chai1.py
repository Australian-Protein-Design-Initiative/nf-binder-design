"""Unit tests for bin/fold/run_chai1.py's input preparation (no chai-lab / GPU needed)."""

import gzip
import json
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin" / "fold"
sys.path.insert(0, str(BIN))

import run_chai1 as rc  # noqa: E402


def write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


def test_msa_rows_key_on_taxid_and_label_sources():
    rows = rc.msa_rows(
        ">query\nPEPTIDE\n"
        ">UniRef100_A0A123 n=1 Tax=Homo sapiens TaxID=9606 RepID=A0A123_HUMAN\nPEPTIDA\n"
        ">tr|Q12345|Q12345_MOUSE desc OS=Mus musculus OX=10090\nPEPTIDQ\n"
        ">MGYP000123\nPEPTIDK\n"
    )
    assert [r["source_database"] for r in rows] == ["query", "uniref90", "bfd_uniclust", "bfd_uniclust"]
    assert [r["pairing_key"] for r in rows] == ["", "9606", "10090", ""]


def test_prepare_inputs_writes_one_msa_per_unique_sequence(tmp_path):
    fasta = write(tmp_path, "q.fasta", ">a\nPEPTIDE\n>b\nPEPTIDE\n>c\nMKV\n")
    a = write(tmp_path, "a.a3m", ">q\nPEPTIDE\n>hit TaxID=9606\nPEPTIDA\n")
    c = write(tmp_path, "c.a3m", ">q\nMKV\n")
    out = rc.prepare_inputs(fasta, tmp_path / "w", [a, a, c])
    assert out["fasta"].read_text() == (
        ">protein|name=A\nPEPTIDE\n>protein|name=B\nPEPTIDE\n>protein|name=C\nMKV\n"
    )
    assert set(out["msas"]) == {rc.pqt_basename("PEPTIDE"), rc.pqt_basename("MKV")}
    assert out["m8"] is None


def test_template_rows_stage_gzipped_cifs(tmp_path):
    tdir = tmp_path / "templates"
    tdir.mkdir()
    (tdir / "tmpl001.cif").write_text("data_tmpl001\n")
    seq = "PEPTIDE"
    import hashlib

    index = {hashlib.md5(seq.encode()).hexdigest(): [{
        "cif": "tmpl001.cif", "template_chain": "A", "identity": 0.9,
        "query_indices": [1, 2, 3], "template_indices": [0, 1, 2],
    }]}
    (tdir / "index.json").write_text(json.dumps(index))
    fasta = write(tmp_path, "q.fasta", ">a\nPEPTIDE\n>b\nMKV\n")
    out = rc.prepare_inputs(fasta, tmp_path / "w", None, tdir, ["A"])
    (row,) = [line.split("\t") for line in out["m8"].read_text().splitlines()]
    assert row[:2] == ["A", "tmpl001_A"]
    assert row[6:10] == ["2", "4", "1", "3"]
    with gzip.open(out["cif_folder"] / "TMPL001.cif.gz", "rt") as f:
        assert f.read() == "data_tmpl001\n"
