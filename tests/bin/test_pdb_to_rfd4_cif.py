#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "gemmi"]
# ///

"""
Tests for bin/rfd4_proteina/pdb_to_rfd4_cif.py.

    uv run --with pytest --with gemmi pytest tests/bin/test_pdb_to_rfd4_cif.py -q

RFD4-Proteina addresses contigs by label_asym_id and label_seq_id, so those are what these
check. A plain `gemmi convert` gets both wrong for this purpose (see the script's docstring).
"""

import importlib.util
from pathlib import Path
from typing import Dict, List

import pytest

gemmi = pytest.importorskip("gemmi")

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "rfd4_proteina" / "pdb_to_rfd4_cif.py"
_spec = importlib.util.spec_from_file_location("pdb_to_rfd4_cif", _MODPATH)
converter = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(converter)

# Two short chains, numbered from 18 and 200, as an author-numbered PDB would be
PDB = """\
ATOM      1  N   ALA A  18      32.394 -38.055 -22.849  1.00 37.28           N
ATOM      2  CA  ALA A  18      31.602 -39.292 -22.650  1.00 35.59           C
ATOM      3  C   ALA A  18      31.623 -40.168 -23.910  1.00 31.72           C
ATOM      4  O   ALA A  18      30.800 -41.075 -24.010  1.00 30.00           O
ATOM      5  N   GLY A  19      32.500 -39.900 -24.800  1.00 30.00           N
ATOM      6  CA  GLY A  19      32.600 -40.700 -26.000  1.00 30.00           C
ATOM      7  C   GLY A  19      33.800 -40.300 -26.800  1.00 30.00           C
ATOM      8  O   GLY A  19      34.000 -40.900 -27.900  1.00 30.00           O
TER
ATOM      9  N   SER B 200      12.000  -8.000  -2.000  1.00 30.00           N
ATOM     10  CA  SER B 200      11.000  -9.000  -2.500  1.00 30.00           C
ATOM     11  C   SER B 200      11.500 -10.000  -3.500  1.00 30.00           C
ATOM     12  O   SER B 200      10.800 -11.000  -3.800  1.00 30.00           O
TER
END
"""


def _atom_site(cif_path: Path) -> List[Dict[str, str]]:
    doc = gemmi.cif.read(str(cif_path))
    block = doc.sole_block()
    tags = [
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "auth_seq_id",
        "auth_asym_id",
        "label_comp_id",
    ]
    table = block.find("_atom_site.", tags)
    return [{tag: row[i] for i, tag in enumerate(tags)} for row in table]


@pytest.fixture()
def converted(tmp_path: Path) -> Path:
    pdb = tmp_path / "in.pdb"
    pdb.write_text(PDB)
    out = tmp_path / "out.cif"
    converter.convert(pdb, out)
    return out


def test_label_asym_id_is_the_chain_name(converted: Path) -> None:
    rows = _atom_site(converted)
    assert {r["label_asym_id"] for r in rows} == {"A", "B"}
    # and it agrees with the author chain on every atom
    assert all(r["label_asym_id"] == r["auth_asym_id"] for r in rows)


def test_label_seq_id_keeps_the_pdb_residue_numbers(converted: Path) -> None:
    rows = _atom_site(converted)
    assert all(r["label_seq_id"] == r["auth_seq_id"] for r in rows)
    assert {r["label_seq_id"] for r in rows} == {"18", "19", "200"}


def test_label_entity_id_is_numeric(converted: Path) -> None:
    # A non-numeric entity id makes atomworks warn that it cannot read the column
    for row in _atom_site(converted):
        assert row["label_entity_id"].isdigit(), row


def test_residues_and_composition_are_unchanged(converted: Path) -> None:
    rows = _atom_site(converted)
    assert len(rows) == 12
    by_chain = {}
    for row in rows:
        by_chain.setdefault(row["label_asym_id"], set()).add(
            (row["label_seq_id"], row["label_comp_id"])
        )
    assert by_chain["A"] == {("18", "ALA"), ("19", "GLY")}
    assert by_chain["B"] == {("200", "SER")}
