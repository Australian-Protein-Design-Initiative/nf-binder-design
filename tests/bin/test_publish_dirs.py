#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
--refold_publish_dir and --fold_publish_dir must be derived from --method in
nextflow.config, never assigned in a workflow script.

Nextflow 26 evaluates the top-level `params.x = ...` lines of every included
script, and main.nf includes all of them, so a workflow-script assignment is
applied whatever --method asks for: the first include (workflows/rfd.nf) sent
every method's fold outputs to results/rfd/refold.

    uv run --with pytest pytest tests/bin/test_publish_dirs.py -q
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PARAMS = ["fold_publish_dir", "refold_publish_dir"]


def test_no_workflow_script_assigns_a_publish_dir() -> None:
    offenders = [
        f"{path.relative_to(ROOT)}: {param}"
        for path in sorted(ROOT.glob("workflows/*.nf"))
        for param in PARAMS
        if re.search(rf"^\s*params\.{param}\s*=", path.read_text(), re.M)
    ]
    assert not offenders, f"assign these in nextflow.config instead: {offenders}"


def test_nextflow_config_derives_both_from_method() -> None:
    text = (ROOT / "nextflow.config").read_text()
    for param in PARAMS:
        assert re.search(rf"^\s*{param}\s*=", text, re.M), f"{param} is not set in nextflow.config"
