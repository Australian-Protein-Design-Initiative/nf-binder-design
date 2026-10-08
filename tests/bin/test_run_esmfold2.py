"""Unit tests for bin/fold/run_esmfold2.py's input building and output shaping.

These cover the parts that run without the esm package / a GPU: FASTA+a3m
pairing, the mmCIF B-factor read-back that feeds ipsae.py's atom_plddts, and the
chain-pair ipTM matrix.
"""

import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin" / "fold"
sys.path.insert(0, str(BIN))

import run_esmfold2 as r2  # noqa: E402


def write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


# --- build_chain_specs ---------------------------------------------------------


def test_single_chain_with_a3m(tmp_path):
    fasta = write(tmp_path, "q.fasta", ">A\nPEPTIDE\n")
    a3m = write(tmp_path, "a.a3m", ">query\nPEPTIDE\n>hit key=9606\nPEPTIDA\n")
    specs = r2.build_chain_specs(fasta, [a3m])
    assert [s.chain_id for s in specs] == ["A"]
    assert specs[0].sequence == "PEPTIDE"
    assert specs[0].a3m == a3m


def test_chain_ids_follow_fasta_record_order(tmp_path):
    fasta = write(tmp_path, "q.fasta", ">target\nPEPTIDE\n>binder\nMKV\n")
    a3m_a = write(tmp_path, "a.a3m", ">q\nPEPTIDE\n")
    a3m_b = write(tmp_path, "b.a3m", ">q\nMKV\n")
    specs = r2.build_chain_specs(fasta, [a3m_a, a3m_b])
    assert [(s.chain_id, s.sequence) for s in specs] == [("A", "PEPTIDE"), ("B", "MKV")]


def test_single_sequence_ignores_a3ms(tmp_path):
    fasta = write(tmp_path, "q.fasta", ">A\nPEPTIDE\n")
    a3m = write(tmp_path, "a.a3m", ">query\nPEPTIDE\n")
    specs = r2.build_chain_specs(fasta, [a3m], single_sequence=True)
    assert specs[0].a3m is None


def test_no_a3m_given_is_single_sequence(tmp_path):
    fasta = write(tmp_path, "q.fasta", ">A\nPEPTIDE\n")
    assert r2.build_chain_specs(fasta, None)[0].a3m is None


def test_query_mismatch_raises(tmp_path):
    fasta = write(tmp_path, "q.fasta", ">A\nPEPTIDE\n")
    a3m = write(tmp_path, "a.a3m", ">wrong\nWWWWWWW\n")
    with pytest.raises(ValueError, match="does not match the FASTA sequence"):
        r2.build_chain_specs(fasta, [a3m])


def test_gapped_query_row_is_accepted(tmp_path):
    """The a3m query may carry gaps; they are stripped before comparison."""
    fasta = write(tmp_path, "q.fasta", ">A\nPEPTIDE\n")
    a3m = write(tmp_path, "a.a3m", ">query\nPEP-TIDE\n")
    assert r2.build_chain_specs(fasta, [a3m])[0].a3m == a3m


def test_wrong_a3m_count_raises(tmp_path):
    fasta = write(tmp_path, "q.fasta", ">A\nPEPTIDE\n>B\nMKV\n")
    a3m = write(tmp_path, "a.a3m", ">q\nPEPTIDE\n")
    with pytest.raises(ValueError, match="one per chain"):
        r2.build_chain_specs(fasta, [a3m])


def test_empty_a3m_falls_back_to_single_sequence(tmp_path):
    fasta = write(tmp_path, "q.fasta", ">A\nPEPTIDE\n")
    a3m = write(tmp_path, "a.a3m", "")
    assert r2.build_chain_specs(fasta, [a3m])[0].a3m is None


def test_empty_fasta_raises(tmp_path):
    with pytest.raises(ValueError, match="No FASTA records"):
        r2.build_chain_specs(write(tmp_path, "q.fasta", ""), None)


# --- atom_plddts_from_cif ------------------------------------------------------


CIF = """data_pred
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N MET A 1 1.0 2.0 3.0 1.0 85.5
ATOM 2 CA MET A 1 1.5 2.5 3.5 1.0 86.0
ATOM 3 C MET A 1 2.0 3.0 4.0 1.0 87.25
#
"""


def test_atom_plddts_in_serial_order():
    assert r2.atom_plddts_from_cif(CIF) == [85.5, 86.0, 87.25]


def test_atom_plddts_indexed_from_lowest_serial():
    """ipsae.py subtracts the lowest serial, so a file starting at 5 must align."""
    shifted = CIF.replace("ATOM 1 N", "ATOM 5 N").replace("ATOM 2 CA", "ATOM 6 CA").replace("ATOM 3 C", "ATOM 7 C")
    assert r2.atom_plddts_from_cif(shifted) == [85.5, 86.0, 87.25]


def test_atom_plddts_without_b_factor_column_raises():
    broken = CIF.replace("_atom_site.B_iso_or_equiv\n", "")
    with pytest.raises(ValueError, match="B_iso_or_equiv"):
        r2.atom_plddts_from_cif(broken)


def test_atom_plddts_without_loop_raises():
    with pytest.raises(ValueError, match="_atom_site loop"):
        r2.atom_plddts_from_cif("data_pred\n#\n")


# --- chain_pair_iptm_matrix ----------------------------------------------------


class FakeTensor:
    def __init__(self, values):
        self._values = values

    def tolist(self):
        return self._values


def test_chain_pair_iptm_matrix_roundtrip():
    matrix = r2.chain_pair_iptm_matrix(FakeTensor([[0.0, 0.8], [0.8, 0.0]]), 2)
    assert matrix == [[0.0, 0.8], [0.8, 0.0]]


def test_chain_pair_iptm_matrix_none():
    assert r2.chain_pair_iptm_matrix(None, 2) is None
