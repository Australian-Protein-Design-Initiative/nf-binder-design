#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "biotite<=1.6", "numpy", "peppr"]
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


def _helix(chain, n_residues, origin, backbone_only=False):
    """A helical chain, so the residues are not collinear.

    The geometry matters: superimpose_homologs keeps every anchor on a straight
    poly-alanine rod, but prunes a real fold down to a core, which is what made it
    the wrong basis for atom correspondence.
    """
    names = ["N", "CA", "C", "O"] if backbone_only else ["N", "CA", "C", "O", "CB"]
    atoms = []
    for i in range(n_residues):
        angle = i * np.deg2rad(100.0)
        centre = np.array(
            [2.3 * np.cos(angle), 2.3 * np.sin(angle), 1.5 * i]
        ) + np.asarray(origin, dtype=float)
        for j, name in enumerate(names):
            atoms.append(
                struc.Atom(
                    centre + np.array([0.0, 0.0, 0.3 * j]),
                    chain_id=chain,
                    res_id=i + 1,
                    res_name="ALA",
                    atom_name=name,
                    element=name[0],
                )
            )
    return atoms


def _two_chain_complex(backbone_only_chain=None, shift=(0.0, 0.0, 0.0)):
    """Two helices side by side, close enough to make a real interface."""
    ligand_origin = np.array([7.0, 0.0, 0.0]) + np.asarray(shift, dtype=float)
    return struc.array(
        _helix("A", 24, [0.0, 0.0, 0.0], backbone_only_chain == "A")
        + _helix("B", 24, ligand_origin, backbone_only_chain == "B")
    )


def test_matching_keeps_every_aligned_residue_not_just_superposition_anchors():
    # superimpose_homologs prunes to a non-outlier core - 15 of 24 residues for this
    # helix, and 23 of 117 for a real designed chain. Atom correspondence has to keep
    # all of them, or fnat counts contacts over a fraction of the interface.
    chain_a = rmsd4all._filter_structure_by_chains(_two_chain_complex(), ["A"])
    anchors = struc.superimpose_homologs(
        chain_a[chain_a.atom_name == "CA"], chain_a[chain_a.atom_name == "CA"]
    )[2]
    assert len(anchors) < 24, "fixture no longer exercises the pruning it guards"

    fixed, mobile = rmsd4all._match_atoms(chain_a, chain_a)
    assert len(fixed) == len(chain_a)
    assert len(mobile) == len(chain_a)


def test_matching_falls_back_to_the_atoms_both_structures_have():
    # A backbone-only design against an all-atom prediction: the shared atom names
    # are the backbone, and both sides must come back the same length.
    design = rmsd4all._filter_structure_by_chains(
        _two_chain_complex(backbone_only_chain="A"), ["A"]
    )
    prediction = rmsd4all._filter_structure_by_chains(_two_chain_complex(), ["A"])
    fixed, mobile = rmsd4all._match_atoms(design, prediction)
    assert len(fixed) == len(mobile)
    assert set(fixed.atom_name.tolist()) == {"N", "CA", "C", "O"}


def test_dockq_is_one_for_a_pose_identical_to_the_reference():
    pytest.importorskip("peppr")
    structure = _two_chain_complex()
    scores = rmsd4all.dockq_scores(structure, structure, ["A"], ["A"], ["B"], ["B"])
    assert scores["dockq"] == pytest.approx(1.0, abs=1e-4)
    assert scores["dockq_fnat"] == pytest.approx(1.0, abs=1e-4)
    assert scores["dockq_irmsd"] == pytest.approx(0.0, abs=1e-4)


def test_dockq_drops_when_the_ligand_is_in_the_wrong_place():
    pytest.importorskip("peppr")
    scores = rmsd4all.dockq_scores(
        _two_chain_complex(),
        _two_chain_complex(shift=(6.0, 4.0, 0.0)),
        ["A"],
        ["A"],
        ["B"],
        ["B"],
    )
    assert scores["dockq"] < 0.5
    assert scores["dockq_lrmsd"] > 1.0
