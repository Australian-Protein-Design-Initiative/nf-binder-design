#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "pandas"]
# ///

"""
Unit tests for bin/parse_boltz_confidence.py (shared between fold.nf's
FOLD_PARSE_BOLTZ_CONFIDENCE and the Boltz refold comparison modules).

    uv run --with pytest --with pandas pytest tests/bin/test_parse_boltz_confidence.py -q
"""

import csv
import io
import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"


def _run(args):
    return subprocess.run([sys.executable, str(BIN / "parse_boltz_confidence.py"), *args],
                           capture_output=True, text=True, check=True)


def _rows(tsv_text):
    return list(csv.DictReader(io.StringIO(tsv_text), delimiter="\t"))


def test_batch_and_msa_depth_columns(tmp_path):
    j = tmp_path / "conf.json"
    j.write_text(json.dumps({"confidence_score": 0.9, "ptm": 0.8}))
    out = _run(["--json", str(j), "--id", "cx", "--model", "0",
                "--batch", "2", "--msa-depth", "512"]).stdout
    r = _rows(out)[0]
    assert r["batch"] == "2" and r["msa_depth"] == "512"
    assert list(r.keys()).index("batch") == list(r.keys()).index("model") + 1


def test_batch_and_msa_depth_omitted_when_unset(tmp_path):
    j = tmp_path / "conf.json"
    j.write_text(json.dumps({"confidence_score": 0.9}))
    out = _run(["--json", str(j), "--id", "cx", "--model", "0"]).stdout
    r = _rows(out)[0]
    assert "batch" not in r and "msa_depth" not in r


def test_merge_ipsae_fills_min_row_columns(tmp_path):
    j = tmp_path / "conf.json"
    j.write_text(json.dumps({"confidence_score": 0.9, "ptm": 0.8}))
    ipsae = tmp_path / "model_10_10_ipsae.tsv"
    ipsae.write_text(
        "Chn1\tChn2\tPAE\tDist\tType\tipSAE\tipSAE_d0chn\tipSAE_d0dom\tipTM_af\tipTM_d0chn"
        "\tpDockQ\tpDockQ2\tLIS\tn0res\tn0chn\tn0dom\td0res\td0chn\td0dom"
        "\tnres1\tnres2\tdist1\tdist2\tModel\n"
        "A\tB\t10\t10\tmin\t0.42\t0.41\t0.40\t0.5\t0.5\t0.3\t0.2\t0.1"
        "\t10\t10\t10\t5\t5\t5\t10\t10\t5\t5\tmodel\n"
    )
    out = _run(["--json", str(j), "--id", "cx", "--merge-ipsae", str(ipsae)]).stdout
    r = _rows(out)[0]
    assert r["ipSAE_min"] == "0.42"
    assert r["ipSAE_d0chn"] == "0.41" and r["ipSAE_d0dom"] == "0.4"
    assert r["pDockQ"] == "0.3" and r["pDockQ2"] == "0.2" and r["LIS"] == "0.1"


def test_merge_ipsae_no_min_row_is_harmless(tmp_path):
    j = tmp_path / "conf.json"
    j.write_text(json.dumps({"confidence_score": 0.9}))
    ipsae = tmp_path / "empty_ipsae.tsv"
    ipsae.write_text(
        "Chn1\tChn2\tPAE\tDist\tType\tipSAE\tipSAE_d0chn\tipSAE_d0dom\tipTM_af\tipTM_d0chn"
        "\tpDockQ\tpDockQ2\tLIS\tn0res\tn0chn\tn0dom\td0res\td0chn\td0dom"
        "\tnres1\tnres2\tdist1\tdist2\tModel\n"
    )
    out = _run(["--json", str(j), "--id", "cx", "--merge-ipsae", str(ipsae)]).stdout
    r = _rows(out)[0]
    assert "ipSAE_min" not in r


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
