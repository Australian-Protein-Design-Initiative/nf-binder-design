#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Unit tests for bin/fold/stockholm_gs_de.py - recovering `#=GS <name> DE ...`
descriptions that alphafold.data.parsers.parse_stockholm drops.

Run from the repo root with uv (host python has no pytest):

    uv run --with pytest pytest tests/bin/test_stockholm_gs_de.py -q
"""

import importlib.util
from pathlib import Path

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "fold" / "stockholm_gs_de.py"
_spec = importlib.util.spec_from_file_location("stockholm_gs_de", _MODPATH)
gsde = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gsde)


# A truncated, real excerpt of examples/fold-multimer/work/.../uniref90_hits.sto
# (query row + two hit rows) - see the module docstring in af2_msas_to_a3m.py.
STOCKHOLM_EXCERPT = """\
# STOCKHOLM 1.0

#=GS UniRef90_UPI0009E0E323/2-223     DE [subseq from] Programmed cell death 1 ligand 1 n=4 Tax=Homo sapiens TaxID=9606 RepID=UPI0009E0E323
#=GS UniRef90_Q9NZQ7/18-239           DE [subseq from] Programmed cell death 1 ligand 1 n=25 Tax=Catarrhini TaxID=9526 RepID=PD1L1_HUMAN

complex.chain_B                          AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQF
UniRef90_UPI0009E0E323/2-223             AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQF
UniRef90_Q9NZQ7/18-239                   AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQF
//
"""


def test_parse_gs_de_extracts_taxid_and_repid():
    gs_de = gsde.parse_gs_de(STOCKHOLM_EXCERPT)
    assert gs_de["UniRef90_UPI0009E0E323/2-223"] == (
        "[subseq from] Programmed cell death 1 ligand 1 n=4 Tax=Homo sapiens "
        "TaxID=9606 RepID=UPI0009E0E323"
    )
    assert "TaxID=9526" in gs_de["UniRef90_Q9NZQ7/18-239"]
    assert "RepID=PD1L1_HUMAN" in gs_de["UniRef90_Q9NZQ7/18-239"]


def test_parse_gs_de_ignores_non_gs_lines():
    gs_de = gsde.parse_gs_de(STOCKHOLM_EXCERPT)
    assert "complex.chain_B" not in gs_de
    assert len(gs_de) == 2


def test_parse_gs_de_empty_text():
    assert gsde.parse_gs_de("") == {}


def test_parse_gs_de_first_de_line_wins():
    text = "#=GS foo DE first\n#=GS foo DE second\n"
    assert gsde.parse_gs_de(text) == {"foo": "first"}


def test_merge_descriptions_appends_de_text_when_present():
    names = ["complex.chain_B", "UniRef90_Q9NZQ7/18-239"]
    gs_de = gsde.parse_gs_de(STOCKHOLM_EXCERPT)
    merged = gsde.merge_descriptions(names, gs_de)
    assert merged[0] == "complex.chain_B"  # query row has no GS DE line -> unchanged
    assert merged[1].startswith("UniRef90_Q9NZQ7/18-239 [subseq from]")
    assert "TaxID=9526" in merged[1]


def test_merge_descriptions_passes_through_names_without_gs_de():
    assert gsde.merge_descriptions(["no_such_name"], {}) == ["no_such_name"]
