#!/usr/bin/env python
# /// script
# requires-python = ">=3.10"
# dependencies = ["pytest", "numpy"]
# ///

"""Tests for bin/ipsae.py helpers (AF2 list-wrapped PAE JSON, Protenix key aliases)."""

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
    already = {"atom_plddts": [50.0], "pae": [[0.0]]}
    assert ipsae.normalize_token_pae_json(already)["atom_plddts"] == [50.0]


def test_summary_confidences_path_protenix(tmp_path):
    pae = tmp_path / "cx_full_data_sample_0.json"
    summary = tmp_path / "cx_summary_confidence_sample_0.json"
    pae.write_text("{}")
    summary.write_text("{}")
    got = ipsae.summary_confidences_path(str(pae))
    assert Path(got).name == "cx_summary_confidence_sample_0.json"
