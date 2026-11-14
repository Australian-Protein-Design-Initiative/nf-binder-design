#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "biotite", "numpy"]
# ///

"""
Tests for bin/rmsd4all.py's chain selection.

    uv run --with pytest --with biotite --with numpy pytest tests/bin/test_rmsd4all.py -q
"""

import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("biotite")

import biotite.structure as struc  # noqa: E402
import numpy as np  # noqa: E402

_BIN = Path(__file__).resolve().parents[2] / "bin"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


rmsd4all = _load("rmsd4all")


def _structure(chain_order):
    """One CA atom per residue, chains laid out in the order given."""
    atoms = []
    for chain, n_residues in chain_order:
        for i in range(n_residues):
            atom = struc.Atom(
                [float(i), 0.0, 0.0],
                chain_id=chain,
                res_id=i + 1,
                res_name="GLY",
                atom_name="CA",
                element="C",
            )
            atoms.append(atom)
    return struc.array(atoms)


def test_chains_come_back_in_the_order_asked_for_not_the_file_order():
    # Written binder-first, as RFD4-Proteina writes its designs.
    structure = _structure([("B", 3), ("A", 5)])
    selected = rmsd4all._filter_structure_by_chains(structure, ["A", "B"])
    assert [str(c) for c in np.unique(selected.chain_id)] == ["A", "B"]
    # The A atoms come first because A was asked for first.
    assert str(selected.chain_id[0]) == "A"
    assert len(selected) == 8


def test_two_structures_written_in_opposite_chain_order_line_up():
    # rmsd_all pairs atoms by position, so the two selections have to agree
    # chain-for-chain however each file happens to be written.
    design = _structure([("B", 3), ("A", 5)])
    prediction = _structure([("A", 5), ("B", 3)])
    fixed = rmsd4all._filter_structure_by_chains(design, ["B", "A"])
    mobile = rmsd4all._filter_structure_by_chains(prediction, ["B", "A"])
    assert [str(c) for c in fixed.chain_id] == [str(c) for c in mobile.chain_id]


def test_a_chain_that_is_not_there_is_skipped_rather_than_failing():
    structure = _structure([("A", 4)])
    selected = rmsd4all._filter_structure_by_chains(structure, ["A", "Z"])
    assert len(selected) == 4


def test_asking_for_no_chain_that_exists_gives_an_empty_selection():
    structure = _structure([("A", 4)])
    selected = rmsd4all._filter_structure_by_chains(structure, ["Z"])
    assert len(selected) == 0


def test_no_chain_list_leaves_the_structure_alone():
    structure = _structure([("B", 3), ("A", 5)])
    assert rmsd4all._filter_structure_by_chains(structure, None) is structure
