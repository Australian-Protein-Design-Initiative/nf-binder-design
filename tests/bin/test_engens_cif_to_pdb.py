#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "biotite"]
# ///

"""
Regression tests for the CIF->PDB converters: the one in
assets/engens/engens-analysis.qmd and the shared bin/structure_to_pdb.py
that bin/cif_to_pdb.py delegates to.

    uv run --with pytest --with biotite pytest tests/bin/test_engens_cif_to_pdb.py -q

OpenFold3 (and RFdiffusion3) mmCIF writers omit _atom_site.occupancy. The
converter used to go through Bio.PDB.MMCIFParser, which hard-requires that
column, so every such prediction raised KeyError, was caught, warned about and
silently dropped from the ensemble. The function source is pulled out of the
.qmd so this tests the shipped code rather than a copy of it.
"""

import logging
import re
import shutil
import textwrap
from pathlib import Path

import numpy
import pytest

_QMD = Path(__file__).resolve().parents[2] / "assets" / "engens" / "engens-analysis.qmd"


def _load_converter():
    """exec the converter helpers out of the qmd, with their imports."""
    src = _QMD.read_text()
    m = re.search(
        r"^def _drop_unmodelled_atoms\(.*?(?=^def _load_mdtraj)", src, re.S | re.M
    )
    assert m, "could not find the converter helpers in the qmd"
    ns = {"Path": Path, "shutil": shutil, "np": numpy, "log": logging.getLogger("test")}
    exec(
        "import biotite.structure.io.pdb as biotite_pdb\n"
        "import biotite.structure.io.pdbx as biotite_pdbx\n" + m.group(0),
        ns,
    )
    return ns["_structure_to_pdb"]


# Two glycine residues, OpenFold3's _atom_site column set: B_iso_or_equiv and
# pdbx_formal_charge present, occupancy absent.
NO_OCCUPANCY_CIF = textwrap.dedent(
    """\
    data_test
    #
    loop_
    _atom_site.group_PDB
    _atom_site.id
    _atom_site.type_symbol
    _atom_site.label_atom_id
    _atom_site.label_alt_id
    _atom_site.label_comp_id
    _atom_site.label_asym_id
    _atom_site.label_entity_id
    _atom_site.label_seq_id
    _atom_site.pdbx_PDB_ins_code
    _atom_site.Cartn_x
    _atom_site.Cartn_y
    _atom_site.Cartn_z
    _atom_site.B_iso_or_equiv
    _atom_site.pdbx_formal_charge
    _atom_site.auth_seq_id
    _atom_site.auth_comp_id
    _atom_site.auth_asym_id
    _atom_site.auth_atom_id
    _atom_site.pdbx_PDB_model_num
    ATOM 1 N N  . GLY A 1 1 ? 0.000 0.000 0.000 90.0 ? 1 GLY A N  1
    ATOM 2 C CA . GLY A 1 1 ? 1.458 0.000 0.000 90.0 ? 1 GLY A CA 1
    ATOM 3 C C  . GLY A 1 1 ? 2.009 1.420 0.000 90.0 ? 1 GLY A C  1
    ATOM 4 O O  . GLY A 1 1 ? 1.251 2.390 0.000 90.0 ? 1 GLY A O  1
    ATOM 5 N N  . GLY A 1 2 ? 3.332 1.552 0.000 90.0 ? 2 GLY A N  1
    ATOM 6 C CA . GLY A 1 2 ? 3.967 2.859 0.000 90.0 ? 2 GLY A CA 1
    ATOM 7 C C  . GLY A 1 2 ? 5.481 2.730 0.000 90.0 ? 2 GLY A C  1
    ATOM 8 O O  . GLY A 1 2 ? 6.000 1.610 0.000 90.0 ? 2 GLY A O  1
    #
    """
)


def test_converts_cif_without_occupancy(tmp_path):
    """The exact failure: no _atom_site.occupancy must not lose the structure."""
    cif = tmp_path / "openfold3_model.cif"
    cif.write_text(NO_OCCUPANCY_CIF)
    assert "_atom_site.occupancy" not in NO_OCCUPANCY_CIF

    out = _load_converter()(cif, tmp_path / "out.pdb")

    assert out.exists() and out.stat().st_size > 0

    import biotite.structure.io.pdb as biotite_pdb

    arr = biotite_pdb.PDBFile.read(str(out)).get_structure(model=1)
    assert arr.array_length() == 8
    # backbone N, CA, C, O survive for mdtraj's "protein and name N CA C O"
    for resid in (1, 2):
        names = set(arr.atom_name[arr.res_id == resid])
        assert {"N", "CA", "C", "O"} <= names, f"residue {resid}: {names}"


def test_biopython_still_fails_on_this_cif(tmp_path):
    """Pins why biotite is used: Bio.PDB.MMCIFParser cannot read this file."""
    Bio = pytest.importorskip("Bio.PDB")
    cif = tmp_path / "openfold3_model.cif"
    cif.write_text(NO_OCCUPANCY_CIF)
    with pytest.raises(KeyError, match="_atom_site.occupancy"):
        Bio.MMCIFParser(QUIET=True).get_structure("x", str(cif))


def test_pdb_input_passes_through(tmp_path):
    """.pdb inputs are copied, not parsed."""
    convert = _load_converter()
    cif = tmp_path / "m.cif"
    cif.write_text(NO_OCCUPANCY_CIF)
    pdb_in = convert(cif, tmp_path / "in.pdb")
    out = convert(pdb_in, tmp_path / "copied.pdb")
    assert out.read_text() == pdb_in.read_text()


# Same two glycines, plus a terminal OXT the engine never modelled: literal
# 'nan' coordinates, exactly as RFdiffusion3/ProteinMPNN write them.
_OXT_ROW = "ATOM 9 O OXT . GLY A 1 2 ? nan nan nan 90.0 ? 2 GLY A OXT 1\n"
assert NO_OCCUPANCY_CIF.endswith("#\n")
NAN_OXT_CIF = NO_OCCUPANCY_CIF[: -len("#\n")] + _OXT_ROW + "#\n"


def test_unmodelled_atom_does_not_lose_the_structure(tmp_path):
    """A NaN placeholder atom must cost that atom, not the whole structure."""
    cif = tmp_path / "rfd3_model.cif"
    cif.write_text(NAN_OXT_CIF)

    out = _load_converter()(cif, tmp_path / "out.pdb")

    import biotite.structure.io.pdb as biotite_pdb

    arr = biotite_pdb.PDBFile.read(str(out)).get_structure(model=1)
    assert "OXT" not in set(arr.atom_name)          # the unplaced atom is gone
    assert arr.array_length() == 8                  # everything else survived
    assert not numpy.isnan(arr.coord).any()
    # the residue OXT belonged to is still present, described by its other atoms
    assert {"N", "CA", "C", "O"} <= set(arr.atom_name[arr.res_id == 2])


def test_nan_fixture_really_is_nan_and_really_breaks_the_writer():
    """Guard the guard: the NaN fixture must actually reach biotite as NaN,
    and must be what makes an unfiltered write fail. Without this, the test
    above can pass for the wrong reason."""
    import biotite.structure.io.pdb as biotite_pdb
    import biotite.structure.io.pdbx as biotite_pdbx
    from biotite.structure import BadStructureError
    import tempfile, os

    with tempfile.TemporaryDirectory() as d:
        cif = os.path.join(d, "m.cif")
        Path(cif).write_text(NAN_OXT_CIF)
        arr = biotite_pdbx.get_structure(biotite_pdbx.CIFFile.read(cif), model=1)
        assert arr.array_length() == 9, "OXT row did not make it into the fixture"
        assert numpy.isnan(arr.coord).any(axis=-1).sum() == 1
        with pytest.raises(BadStructureError):
            f = biotite_pdb.PDBFile()
            f.set_structure(arr)          # unfiltered: this is the failure we fix


# --- bin/structure_to_pdb.py: same two defects, the tracked implementation ---

def _load_structure_to_pdb():
    import importlib.util

    path = Path(__file__).resolve().parents[2] / "bin" / "structure_to_pdb.py"
    spec = importlib.util.spec_from_file_location("structure_to_pdb", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_structure_to_pdb_reads_cif_without_occupancy(tmp_path):
    mod = _load_structure_to_pdb()
    cif = tmp_path / "m.cif"
    cif.write_text(NO_OCCUPANCY_CIF)
    arr = mod.read_structure(cif)
    assert arr.array_length() == 8


def test_structure_to_pdb_drops_unmodelled_atoms(tmp_path):
    """The refold scoring path must not die on a placeholder OXT."""
    mod = _load_structure_to_pdb()
    cif = tmp_path / "m.cif"
    cif.write_text(NAN_OXT_CIF)

    arr = mod.drop_unmodelled_atoms(mod.read_structure(cif), str(cif))

    assert arr.array_length() == 8
    assert "OXT" not in set(arr.atom_name)
    assert not numpy.isnan(arr.coord).any()
    # and the result is now writable as PDB, which is the point
    import biotite.structure.io.pdb as biotite_pdb

    f = biotite_pdb.PDBFile()
    f.set_structure(arr)
    f.write(str(tmp_path / "out.pdb"))
    assert (tmp_path / "out.pdb").stat().st_size > 0


def test_structure_to_pdb_leaves_clean_structures_alone(tmp_path):
    mod = _load_structure_to_pdb()
    cif = tmp_path / "m.cif"
    cif.write_text(NO_OCCUPANCY_CIF)
    arr = mod.read_structure(cif)
    assert mod.drop_unmodelled_atoms(arr, "x") is arr
