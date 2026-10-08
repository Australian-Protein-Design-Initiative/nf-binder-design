#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Unit tests for bin/fold/make_openfold3_input.py.

    uv run --with pytest pytest tests/bin/test_make_openfold3_input.py -q
"""

import importlib.util
import sys
from pathlib import Path

import pytest

_BIN = Path(__file__).resolve().parents[2] / "bin" / "fold"
# make_openfold3_input imports its sibling make_af3_input, as it does when run
# as a script from bin/fold.
sys.path.insert(0, str(_BIN))
_spec = importlib.util.spec_from_file_location("make_openfold3_input", _BIN / "make_openfold3_input.py")
moi = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(moi)

SEQ_A = "MKTAYIAKQR"
SEQ_B = "GSHMLEDPV"


def _write(p: Path, text: str) -> Path:
    p.write_text(text)
    return p


def test_monomer_query_only_when_no_msa(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">x\n{SEQ_A}\n")
    spec = moi.make_openfold3_input(fasta, "my job/1", tmp_path)
    (name, query), = spec["queries"].items()
    assert name == "my_job1"
    (chain,) = query["chains"]
    assert chain["molecule_type"] == "protein" and chain["chain_ids"] == ["A"]
    (msa_dir,) = chain["main_msa_file_paths"]
    assert (tmp_path / msa_dir / "colabfold_main.a3m").read_text() == f">query\n{SEQ_A}\n"
    assert not (tmp_path / msa_dir / "uniprot_hits.a3m").exists()


def test_multimer_msas_and_pairing_files(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">a\n{SEQ_A}\n>b\n{SEQ_B}\n")
    a = _write(tmp_path / "a.a3m", f"#10\t1\n>a\n{SEQ_A}\n>hit\nMKTAYIAKQ-\n")
    b = _write(tmp_path / "b.a3m", f">b\n{SEQ_B}\n")
    pa = _write(tmp_path / "a.pair.a3m", f">a\n{SEQ_A}\n>tr|P12345|P12345_HUMAN/1-9\nMKTAYIAKQ-\n")
    pb = _write(tmp_path / "b.pair.a3m", f">b\n{SEQ_B}\n")
    spec = moi.make_openfold3_input(fasta, "cx", tmp_path, [a, b], [pa, pb])
    chains = spec["queries"]["cx"]["chains"]
    dirs = [c["main_msa_file_paths"][0] for c in chains]
    assert len(set(dirs)) == 2
    main_a = (tmp_path / dirs[0] / "colabfold_main.a3m").read_text()
    assert main_a.startswith(f">a\n{SEQ_A}\n") and "#" not in main_a
    assert "P12345_HUMAN" in (tmp_path / dirs[0] / "uniprot_hits.a3m").read_text()


def test_identical_chains_share_msa_dir(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">a\n{SEQ_A}\n>b\n{SEQ_A}\n")
    spec = moi.make_openfold3_input(fasta, "homo", tmp_path)
    dirs = [c["main_msa_file_paths"][0] for c in spec["queries"]["homo"]["chains"]]
    assert dirs[0] == dirs[1]


def test_query_mismatch_rejected(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">a\n{SEQ_A}\n>b\n{SEQ_B}\n")
    a = _write(tmp_path / "a.a3m", f">a\n{SEQ_A}\n")
    b = _write(tmp_path / "b.a3m", f">b\n{SEQ_B}\n")
    with pytest.raises(ValueError, match="query"):
        moi.make_openfold3_input(fasta, "x", tmp_path, [b, a])


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))


def _templates_dir(tmp_path: Path, seq: str) -> Path:
    import hashlib
    import json

    d = tmp_path / "fold_templates"
    d.mkdir()
    (d / "tmpl001.cif").write_text("data_tmpl001\n")
    hit = {"cif": "tmpl001.cif", "template_chain": "A", "query_indices": [0], "template_indices": [0]}
    (d / "index.json").write_text(json.dumps({hashlib.md5(seq.encode()).hexdigest(): [hit]}))
    return d


def test_cif_direct_templates_for_matching_chain(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">a\n{SEQ_A}\n>b\n{SEQ_B}\n")
    out = tmp_path / "out"
    out.mkdir()
    spec = moi.make_openfold3_input(fasta, "x", out, templates_dir=_templates_dir(tmp_path, SEQ_A))
    a, b = next(iter(spec["queries"].values()))["chains"]
    assert a["template_cif_paths"] == [f"{moi.TEMPLATE_DIR_NAME}/tmpl001.cif"]
    assert a["template_cif_chain_ids"] == ["A"]
    assert "template_cif_paths" not in b
    assert (out / moi.TEMPLATE_DIR_NAME / "tmpl001.cif").exists()


def test_template_dir_always_created(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">a\n{SEQ_A}\n")
    out = tmp_path / "out"
    out.mkdir()
    spec = moi.make_openfold3_input(fasta, "x", out)
    assert "template_cif_paths" not in next(iter(spec["queries"].values()))["chains"][0]
    assert (out / moi.TEMPLATE_DIR_NAME).is_dir()


def test_template_chains_restricts(tmp_path):
    fasta = _write(tmp_path / "in.fasta", f">a\n{SEQ_A}\n>b\n{SEQ_A}\n")
    out = tmp_path / "out"
    out.mkdir()
    spec = moi.make_openfold3_input(fasta, "x", out, templates_dir=_templates_dir(tmp_path, SEQ_A), template_chains=["A"])
    a, b = next(iter(spec["queries"].values()))["chains"]
    assert "template_cif_paths" in a and "template_cif_paths" not in b
