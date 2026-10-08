#!/usr/bin/env python
# /// script
# requires-python = ">=3.10"
# dependencies = ["pytest", "numpy"]
# ///

"""Tests for bin/ipsae.py helpers (structure extension detection, AF2 list-wrapped
PAE JSON, Protenix / OpenFold3 key aliases and OpenFold3 summaries)."""

import importlib.util
from pathlib import Path

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "ipsae.py"
_spec = importlib.util.spec_from_file_location("ipsae", _MODPATH)
ipsae = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ipsae)


def test_unwrap_af2_list_wrapped_pae():
    inner = {"predicted_aligned_error": [[0.1, 0.2], [0.2, 0.1]]}
    assert ipsae.unwrap_json_object([inner]) == inner
    assert ipsae.unwrap_json_object(inner) is inner


def test_normalize_protenix_full_data_keys():
    raw = {"atom_plddt": [0.9, 0.8], "token_pair_pae": [[0.1, 1.0], [1.0, 0.1]]}
    out = ipsae.normalize_token_pae_json(raw)
    assert out["atom_plddts"] == [0.9, 0.8]
    assert out["pae"] == [[0.1, 1.0], [1.0, 0.1]]
    # AF3/RF3 names are left alone
    already = {"atom_plddts": [50.0], "pae": [[0.0]]}
    assert ipsae.normalize_token_pae_json(already)["atom_plddts"] == [50.0]


def test_summary_confidences_path_protenix(tmp_path):
    pae = tmp_path / "cx_full_data_sample_0.json"
    summary = tmp_path / "cx_summary_confidence_sample_0.json"
    pae.write_text("{}")
    summary.write_text("{}")
    got = ipsae.summary_confidences_path(str(pae))
    assert Path(got).name == "cx_summary_confidence_sample_0.json"


def test_normalize_openfold3_plddt_key():
    raw = {"plddt": [91.0, 88.0], "pde": [[0.1]], "pae": [[0.2]]}
    assert ipsae.normalize_token_pae_json(raw)["atom_plddts"] == [91.0, 88.0]


def test_openfold3_format_detected_from_aggregated_sibling(tmp_path):
    pae = tmp_path / "cx_seed_42_sample_1_confidences.json"
    agg = tmp_path / "cx_seed_42_sample_1_confidences_aggregated.json"
    pae.write_text("{}")
    assert ipsae.resolve_input_format("auto", "cx_seed_42_sample_1_model.cif", str(pae)) == "af3"
    agg.write_text("{}")
    assert ipsae.openfold3_aggregated_path(str(pae)) == str(agg)
    assert ipsae.resolve_input_format("auto", "cx_seed_42_sample_1_model.cif", str(pae)) == "openfold3"


def test_openfold3_chain_pair_iptm_dict():
    summary = {"chain_pair_iptm": {"(A, B)": 0.92, "(A, C)": 0.4}}
    iptm = ipsae.openfold3_chain_pair_iptm(summary, ["A", "B", "C"])
    assert iptm["A"]["B"] == iptm["B"]["A"] == 0.92
    assert iptm["C"]["A"] == 0.4 and iptm["B"]["C"] == 0.0


def test_split_structure_name_uses_real_extension():
    assert ipsae.split_structure_name("model_0.cif") == ("model_0", True)
    assert ipsae.split_structure_name("model_0.pdb") == ("model_0", False)
    assert ipsae.split_structure_name("model_0.npz") is None


def test_split_structure_name_ignores_embedded_extension():
    # rfd3 -> Boltz names carry the RFdiffusion3 backbone filename, extension and
    # all, in the middle of the design id. This must still read as a PDB.
    name = "epea_tipfix_0_model_0.cif_b0_d1_model_0.pdb"
    assert ipsae.split_structure_name(name) == (
        "epea_tipfix_0_model_0.cif_b0_d1_model_0",
        False,
    )


def test_resolve_input_format_embedded_cif_in_pdb_name():
    name = "epea_tipfix_0_model_0.cif_b0_d1_model_0.pdb"
    assert ipsae.resolve_input_format("auto", name, "whatever_pae.json") == "af2"
