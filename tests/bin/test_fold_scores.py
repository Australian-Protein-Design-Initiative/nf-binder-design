#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "numpy"]
# ///

"""
Unit tests for the fold.nf score-TSV pipeline:
  bin/fold/parse_fold_confidence.py  (rf3 / protenix / af2 / af2_mono -> normalized row)
  bin/fold/merge_fold_scores.py      (per-tool TSVs -> master, boltz mapped)

Run from the repo root (host python has no pytest):

    uv run --with pytest --with numpy pytest tests/bin/test_fold_scores.py -q

Fixtures are inline (the example results/ tree is gitignored), mirroring each
engine's real confidence-JSON schema. Assertions lock the normalized schema:
plddt is rescaled to 0-1, equivalent scores share a column name, and asymmetric
per-chain-pair values are dropped.
"""

import csv
import io
import json
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np

BIN = Path(__file__).resolve().parents[2] / "bin" / "fold"
CANON = [
    "tool", "id", "model", "batch", "msa_depth", "original_file", "predictions_file",
    "ranking_score", "ptm", "iptm", "plddt", "pae", "pde", "has_clash",
    "ipsae", "ipsae_d0chn", "ipsae_d0dom", "pdockq", "pdockq2", "lis",
]


def _run(args, **kw):
    return subprocess.run([sys.executable, *args], capture_output=True, text=True, check=True, **kw)


def _rows(tsv_text):
    return list(csv.DictReader(io.StringIO(tsv_text), delimiter="\t"))


def _parse(tmp_path, tool, payload, **flags):
    j = tmp_path / "conf.json"
    j.write_text(json.dumps(payload))
    args = [str(BIN / "parse_fold_confidence.py"), "--tool", tool, "--id", "cx",
            "--model", "m0", "--original-file", "s.cif", "--predictions-file",
            f"{tool}_s.cif", "--json", str(j)]
    for k, v in flags.items():
        args += [f"--{k}", str(v)]
    return _run(args).stdout


def test_rf3_normalized(tmp_path):
    payload = {"ranking_score": 0.80, "ptm": 0.78, "iptm": 0.81,
               "overall_plddt": 0.819, "overall_pae": 9.4, "overall_pde": 2.5,
               "has_clash": False, "chain_ptm": [0.83, 0.8],
               "chain_pair_pae": [[None, 11.9], [None, None]]}
    rows = _rows(_parse(tmp_path, "rf3", payload))
    assert list(rows[0].keys()) == CANON            # exact canonical schema
    r = rows[0]
    assert r["tool"] == "rf3" and r["plddt"] == "0.819" and r["pae"] == "9.4"
    assert r["has_clash"] == "false"
    assert "chain_ptm" not in r and "chain_pair_pae" not in r   # asymmetric dropped


def test_protenix_plddt_rescaled(tmp_path):
    payload = {"ranking_score": 0.87, "ptm": 0.87, "iptm": 0.87, "plddt": 88.15,
               "gpde": 0.43, "has_clash": False,
               "chain_pair_iptm": [[0.0, 0.87], [0.87, 0.0]]}
    r = _rows(_parse(tmp_path, "protenix", payload))[0]
    assert abs(float(r["plddt"]) - 0.8815) < 1e-9   # 0-100 -> 0-1
    assert r["pde"] == "0.43" and "chain_pair_iptm" not in r


def test_af3_summary_and_full_json(tmp_path):
    payload = {"ranking_score": 0.71, "ptm": 0.66, "iptm": 0.74, "has_clash": 0.0,
               "fraction_disordered": 0.02, "chain_pair_iptm": [[0.8, 0.74], [0.74, 0.7]]}
    full = tmp_path / "full.json"
    full.write_text(json.dumps({"atom_plddts": [80.0, 90.0], "pae": [[1.0, 3.0], [5.0, 7.0]]}))
    r = _rows(_parse(tmp_path, "af3", payload, **{"full-json": full}))[0]
    assert r["tool"] == "af3" and r["iptm"] == "0.74"
    assert abs(float(r["plddt"]) - 0.85) < 1e-9     # 0-100 -> 0-1
    assert abs(float(r["pae"]) - 4.0) < 1e-9
    assert r["has_clash"] == "false"


def test_af3_without_full_json_leaves_plddt_blank(tmp_path):
    r = _rows(_parse(tmp_path, "af3", {"ranking_score": 0.5, "ptm": 0.5, "iptm": 0.4, "has_clash": 1.0}))[0]
    assert r["plddt"] == "" and r["has_clash"] == "true"


def test_openfold3_aggregated_and_full_json(tmp_path):
    payload = {"avg_plddt": 82.5, "gpde": 1.2, "iptm": 0.61, "ptm": 0.7, "disorder": 0.1,
               "has_clash": 0.0, "sample_ranking_score": 0.64,
               "chain_ptm": {"A": 0.7}, "chain_pair_iptm": {"(A, B)": 0.61}}
    full = tmp_path / "full.json"
    full.write_text(json.dumps({"plddt": [80.0, 85.0], "pde": [[0.5]], "pae": [[2.0, 4.0], [6.0, 8.0]]}))
    r = _rows(_parse(tmp_path, "openfold3", payload, **{"full-json": full}))[0]
    assert r["tool"] == "openfold3" and r["ranking_score"] == "0.64" and r["iptm"] == "0.61"
    assert abs(float(r["plddt"]) - 0.825) < 1e-9
    assert r["pde"] == "1.2" and abs(float(r["pae"]) - 5.0) < 1e-9
    assert r["has_clash"] == "false"


def test_rf3_ipsae_tsv_merged(tmp_path):
    payload = {"ranking_score": 0.80, "ptm": 0.78, "iptm": 0.81,
               "overall_plddt": 0.819, "overall_pae": 9.4, "overall_pde": 2.5,
               "has_clash": False}
    ipsae = tmp_path / "model_10_10_ipsae.tsv"
    ipsae.write_text(
        "\nChn1 Chn2 PAE Dist Type ipSAE ipSAE_d0chn ipSAE_d0dom pDockQ pDockQ2 LIS\n"
        "A B 10 10 min 0.42 0.41 0.40 0.3 0.2 0.1\n"
    )
    r = _rows(_parse(tmp_path, "rf3", payload, **{"ipsae-tsv": ipsae}))[0]
    assert r["ipsae"] == "0.42" and r["ipsae_d0chn"] == "0.41"
    assert r["pdockq2"] == "0.2" and r["lis"] == "0.1"
    # canonical af2/rf3 table
    canon = tmp_path / "rf3_fold_scores.tsv"
    with canon.open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(CANON)
        w.writerow(["rf3", "cx", "s0", "1", "", "o.cif", "rf3_o.cif", "0.8", "0.78",
                    "0.81", "0.819", "9.4", "2.5", "false", "", "", "", "", "", ""])
    # native boltz table (batch/msa_depth optional columns present here)
    boltz = tmp_path / "boltz_fold_scores.tsv"
    with boltz.open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "model", "batch", "msa_depth", "original_file", "predictions_file",
                    "confidence_score", "ptm", "iptm", "complex_plddt",
                    "complex_pde", "ipsae_min", "ipSAE_d0chn", "ipSAE_d0dom",
                    "pDockQ", "pDockQ2", "LIS", "pair_chains_iptm_0_1"])
        w.writerow(["cx", "0", "2", "512", "cm0.cif", "boltz_cm0.cif", "0.875", "0.82",
                    "0.88", "0.874", "0.48", "0.63", "0.60", "0.59", "0.31", "0.22", "0.15", "0.72"])
    out = tmp_path / "master.tsv"
    _run([str(BIN / "merge_fold_scores.py"), "--input", str(canon),
          "--input", str(boltz), "-o", str(out)])
    rows = _rows(out.read_text())
    assert [r["tool"] for r in rows] == ["boltz", "rf3"]     # sorted by tool
    b = next(r for r in rows if r["tool"] == "boltz")
    assert b["ranking_score"] == "0.875" and b["plddt"] == "0.874"
    assert b["pde"] == "0.48" and b["ipsae"] == "0.63"
    assert b["ipsae_d0chn"] == "0.60" and b["ipsae_d0dom"] == "0.59"
    assert b["pdockq"] == "0.31" and b["pdockq2"] == "0.22" and b["lis"] == "0.15"
    assert b["batch"] == "2" and b["msa_depth"] == "512"
    assert b["predictions_file"] == "boltz_cm0.cif"
    r_rf3 = next(r for r in rows if r["tool"] == "rf3")
    assert r_rf3["batch"] == "1" and r_rf3["msa_depth"] == ""
    assert list(rows[0].keys()) == CANON                    # master is canonical


def test_batch_and_msa_depth_columns(tmp_path):
    payload = {"ranking_score": 0.5, "ptm": 0.5, "iptm": 0.5}
    r = _rows(_parse(tmp_path, "rf3", payload, **{"batch": "2", "msa-depth": "512"}))[0]
    assert r["batch"] == "2" and r["msa_depth"] == "512"
    # unset -> blank, not omitted (fixed column position in every row)
    r2 = _rows(_parse(tmp_path, "rf3", payload))[0]
    assert r2["batch"] == "" and r2["msa_depth"] == ""


def test_af2_mono_uses_af2_parser(tmp_path):
    pkl = tmp_path / "result_model_1.pkl"
    with pkl.open("wb") as f:
        pickle.dump({
            "ranking_confidence": 0.921,
            "ptm": 0.81,
            "plddt": np.array([90.0, 92.0]),
        }, f)
    out = _run([
        str(BIN / "parse_fold_confidence.py"),
        "--tool", "af2_mono", "--id", "cx", "--model", "1",
        "--original-file", "relaxed_model_1_ptm_pred_0.cif",
        "--predictions-file", "af2_mono_cx_run1_relaxed_model_1_ptm_pred_0.cif",
        "--pkl", str(pkl),
    ]).stdout
    r = _rows(out)[0]
    assert r["tool"] == "af2_mono"
    assert r["ranking_score"] == "0.921"
    assert r["ptm"] == "0.81"
    assert r["iptm"] == ""  # monomer_ptm pickle has no iptm
    assert abs(float(r["plddt"]) - 0.91) < 1e-9  # mean([90, 92]) / 100


def test_esmfold2_summary_and_full_json(tmp_path):
    """run_esmfold2.py writes pLDDT already on 0-1 and reports no ranking score."""
    payload = {"ptm": 0.72, "iptm": 0.64, "plddt": 0.813,
               "chain_ids": ["A", "B"], "seed": 42, "sample": 0,
               "chain_pair_iptm": [[0.0, 0.64], [0.64, 0.0]]}
    full = tmp_path / "full.json"
    full.write_text(json.dumps({"pae": [[2.0, 4.0], [6.0, 8.0]], "atom_plddts": [81.0, 82.0]}))
    r = _rows(_parse(tmp_path, "esmfold2", payload, **{"full-json": full}))[0]
    assert r["tool"] == "esmfold2" and r["ptm"] == "0.72" and r["iptm"] == "0.64"
    assert abs(float(r["plddt"]) - 0.813) < 1e-9
    assert abs(float(r["pae"]) - 5.0) < 1e-9
    assert r["ranking_score"] == "" and r["pde"] == "" and r["has_clash"] == ""


def test_esmfold2_fast_uses_esmfold2_parser_with_own_tool_tag(tmp_path):
    r = _rows(_parse(tmp_path, "esmfold2_fast", {"ptm": 0.8, "iptm": 0.5, "plddt": 0.7}))[0]
    assert r["tool"] == "esmfold2_fast" and r["iptm"] == "0.5" and r["ranking_score"] == ""


def test_esmfold2_monomer_leaves_iptm_blank(tmp_path):
    r = _rows(_parse(tmp_path, "esmfold2", {"ptm": 0.9, "iptm": None, "plddt": 0.88}))[0]
    assert r["iptm"] == "" and r["ptm"] == "0.9" and r["pae"] == ""
