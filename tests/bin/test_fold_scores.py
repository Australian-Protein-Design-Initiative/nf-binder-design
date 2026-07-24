#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "numpy"]
# ///

"""
Unit tests for the fold.nf score-TSV pipeline:
  bin/parse_fold_confidence.py  (rf3 / protenix / af2 -> normalized row)
  bin/merge_fold_scores.py      (per-tool TSVs -> master, boltz mapped)

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
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
CANON = [
    "tool", "id", "model", "original_file", "predictions_file",
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


def test_merge_maps_boltz_and_concats(tmp_path):
    # canonical af2/rf3 table
    canon = tmp_path / "rf3_fold_scores.tsv"
    with canon.open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(CANON)
        w.writerow(["rf3", "cx", "s0", "o.cif", "rf3_o.cif", "0.8", "0.78",
                    "0.81", "0.819", "9.4", "2.5", "false", "", "", "", "", "", ""])
    # native boltz table
    boltz = tmp_path / "boltz_fold_scores.tsv"
    with boltz.open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "model", "original_file", "predictions_file",
                    "confidence_score", "ptm", "iptm", "complex_plddt",
                    "complex_pde", "ipsae_min", "pair_chains_iptm_0_1"])
        w.writerow(["cx", "0", "cm0.cif", "boltz_cm0.cif", "0.875", "0.82",
                    "0.88", "0.874", "0.48", "0.63", "0.72"])
    out = tmp_path / "master.tsv"
    _run([str(BIN / "merge_fold_scores.py"), "--input", str(canon),
          "--input", str(boltz), "-o", str(out)])
    rows = _rows(out.read_text())
    assert [r["tool"] for r in rows] == ["boltz", "rf3"]     # sorted by tool
    b = next(r for r in rows if r["tool"] == "boltz")
    assert b["ranking_score"] == "0.875" and b["plddt"] == "0.874"
    assert b["pde"] == "0.48" and b["ipsae"] == "0.63"
    assert b["predictions_file"] == "boltz_cm0.cif"
    assert list(rows[0].keys()) == CANON                    # master is canonical
