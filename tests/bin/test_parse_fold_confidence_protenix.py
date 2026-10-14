#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "numpy"]
# ///

"""
Unit tests for the protenix branch of bin/fold/parse_fold_confidence.py.

Regression cover for an empty `pae` column on every protenix fold: the parser
did not open the full-data JSON, and protenix names its NxN matrix
`token_pair_pae` rather than `pae` as af3/openfold3/esmfold2 do. The result was
a blank column, which makes Germinal's `external_pae` gate unevaluable -- and an
unevaluable gate is not a pass, so accept rates collapsed to zero for reasons
unrelated to the structures.

    uv run --with pytest --with numpy pytest tests/bin/test_parse_fold_confidence_protenix.py -q
"""

import csv
import io
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "fold" / "parse_fold_confidence.py"

# A protenix summary_confidence JSON carries no whole-structure PAE, only the
# per-chain-pair means. This is the shape the real files have.
SUMMARY = {
    "plddt": 92.0,            # protenix reports pLDDT on 0-100
    "gpde": 0.4,
    "ptm": 0.9,
    "iptm": 0.87,
    "ranking_score": 0.88,
    "chain_pair_pae_mean": [[2.0, 5.0], [6.0, 3.0]],
    "has_clash": False,
}


def _run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args],
                          capture_output=True, text=True, check=True).stdout


def _row(tsv_text):
    return list(csv.DictReader(io.StringIO(tsv_text), delimiter="\t"))[0]


def _write(tmp_path, summary=SUMMARY, full=None):
    s = tmp_path / "x_summary_confidence_sample_0.json"
    s.write_text(json.dumps(summary))
    paths = ["--tool", "protenix", "--id", "x", "--model", "sample_0",
             "--json", str(s)]
    if full is not None:
        f = tmp_path / "x_full_data_sample_0.json"
        f.write_text(json.dumps(full))
        paths += ["--full-json", str(f)]
    return paths


def test_pae_is_mean_of_the_whole_token_matrix(tmp_path):
    # Deliberately asymmetric: a mean over the whole matrix is 2.5, while a mean
    # of only the upper triangle or only the diagonal would not be.
    full = {"token_pair_pae": [[1.0, 2.0], [3.0, 4.0]]}
    row = _row(_run(*_write(tmp_path, full=full)))
    assert float(row["pae"]) == 2.5


def test_pae_absent_without_full_json(tmp_path):
    """The flag stays optional, so callers that do not stage the file still work."""
    row = _row(_run(*_write(tmp_path)))
    assert row["pae"] == ""


def test_summary_chain_pair_pae_is_not_used_as_pae(tmp_path):
    """chain_pair_pae_mean is per chain pair, not the whole-structure PAE.

    Germinal's external_pae is a mean over the entire NxN matrix, so silently
    substituting the summary's chain-pair block would report a different
    quantity under the same column name.
    """
    row = _row(_run(*_write(tmp_path)))
    assert row["pae"] == ""


def test_other_protenix_fields_still_parse(tmp_path):
    row = _row(_run(*_write(tmp_path, full={"token_pair_pae": [[1.0]]})))
    assert float(row["plddt"]) == 0.92      # rescaled from 0-100 to 0-1
    assert float(row["ptm"]) == 0.9
    assert float(row["iptm"]) == 0.87
    assert float(row["pde"]) == 0.4
