#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""Tests for bin/fold/assemble_af2_multimer_msas.py (stdlib; no AF2 container)."""

import importlib.util
from pathlib import Path

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "fold" / "assemble_af2_multimer_msas.py"
_spec = importlib.util.spec_from_file_location("assemble_af2_multimer_msas", _MODPATH)
asm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(asm)


def test_a3m_to_stockholm_uniquifies_duplicate_ids() -> None:
    records = [
        ("101", "ACDE"),
        ("hitA", "AC-E"),
        ("101", "A-DE"),
        ("hitA extra", "ACDE"),
    ]
    sto = asm.a3m_to_stockholm(records)
    names = [
        line.split()[0]
        for line in sto.splitlines()
        if line and not line.startswith(("#", "//"))
    ]
    assert names == ["101", "hitA", "101_2", "hitA_2"]
    assert len(names) == len(set(names))
