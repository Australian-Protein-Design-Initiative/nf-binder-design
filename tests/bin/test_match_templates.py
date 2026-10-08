#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "gemmi", "biopython"]
# ///

"""
Unit tests for bin/fold/match_templates.py.

    uv run --with pytest --with gemmi --with biopython pytest tests/bin/test_match_templates.py -q
"""

import importlib.util
import json
from pathlib import Path

import gemmi
import pytest

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("match_templates", ROOT / "bin" / "fold" / "match_templates.py")
mt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mt)

INPUT = ROOT / "examples" / "fold-pulldown" / "input"
PDL1_PDB = INPUT / "PDL1.pdb"
IL7RA_PDB = INPUT / "IL7RA.pdb"


def _targets():
    return dict(mt.parse_fasta(INPUT / "targets.fasta"))


def _run(tmp_path, templates, queries, **kw):
    opts = dict(min_identity=0.3, min_coverage=0.3, min_aligned=40, max_per_chain=4)
    opts.update(kw)
    accepted, report = mt.match_templates(templates, queries, **opts)
    out = tmp_path / "out"
    mt.write_outputs(accepted, report, out)
    return json.loads((out / "index.json").read_text()), report, out


def test_each_target_matches_its_own_structure(tmp_path):
    t = _targets()
    index, _, _ = _run(tmp_path, [PDL1_PDB, IL7RA_PDB], list(t.items()))
    pdl1 = index[mt.seq_key(t["PDL1"])]
    il7ra = index[mt.seq_key(t["IL7RA"])]
    assert [h["source"] for h in pdl1] == ["PDL1.pdb:A"]
    assert [h["source"] for h in il7ra] == ["IL7RA.pdb:A"]
    assert pdl1[0]["identity"] == 1.0 and pdl1[0]["coverage"] == 1.0


def test_mapping_indexes_matching_residues(tmp_path):
    seq = _targets()["PDL1"]
    query = seq[:20] + seq[26:]  # a deletion relative to the template
    index, _, out = _run(tmp_path, [PDL1_PDB], [("del", query)])
    hit = index[mt.seq_key(query)][0]
    tmpl_seq = gemmi.one_letter_code([r.name for r in gemmi.read_structure(str(out / hit["cif"]))[0]["A"]])
    pairs = list(zip(hit["query_indices"], hit["template_indices"]))
    assert len(pairs) == len(query)
    assert all(query[q] == tmpl_seq[t] for q, t in pairs)


def test_thresholds_reject_weak_matches(tmp_path):
    t = _targets()
    index, report, _ = _run(tmp_path, [PDL1_PDB], [("IL7RA", t["IL7RA"])], min_identity=0.5)
    assert index == {}
    assert report[0][6] == "rejected"


def test_short_alignment_rejected(tmp_path):
    """A partial match passing identity and coverage is still too short to be a fold."""
    seq = _targets()["PDL1"]
    query = seq[:30] + "WWEWWKWWDWWRWWNWWQWWHWWYWWCWWM"
    index, report, _ = _run(tmp_path, [PDL1_PDB], [("frag", query)])
    assert index == {}
    assert report[0][6] == "rejected" and report[0][7].startswith("aligned <")


def test_short_query_matched_end_to_end(tmp_path):
    """Below min_aligned the whole query must align, so peptides stay templatable."""
    seq = _targets()["PDL1"]
    query = seq[:25]
    index, _, _ = _run(tmp_path, [PDL1_PDB], [("peptide", query)])
    assert len(index[mt.seq_key(query)]) == 1


def test_max_per_chain_keeps_best(tmp_path):
    copy = tmp_path / "PDL1_copy.pdb"
    copy.write_bytes(PDL1_PDB.read_bytes())
    t = _targets()
    index, report, _ = _run(tmp_path, [PDL1_PDB, copy], [("PDL1", t["PDL1"])], max_per_chain=1)
    assert len(index[mt.seq_key(t["PDL1"])]) == 1
    assert any(r[7] == "rank > 1" for r in report)


def test_complex_template_split_by_chain(tmp_path):
    st = gemmi.read_structure(str(PDL1_PDB))
    other = gemmi.read_structure(str(IL7RA_PDB))[0][0]
    other.name = "B"
    st[0].add_chain(other)
    cplx = tmp_path / "complex.pdb"
    st.write_pdb(str(cplx))
    t = _targets()
    index, _, _ = _run(tmp_path, [cplx], list(t.items()))
    assert index[mt.seq_key(t["PDL1"])][0]["source"] == "complex.pdb:A"
    assert index[mt.seq_key(t["IL7RA"])][0]["source"] == "complex.pdb:B"


def test_written_cif_is_single_chain_with_revision_date(tmp_path):
    t = _targets()
    index, _, out = _run(tmp_path, [PDL1_PDB], [("PDL1", t["PDL1"])])
    cif = out / index[mt.seq_key(t["PDL1"])][0]["cif"]
    block = gemmi.cif.read(str(cif)).sole_block()
    assert block.find_value("_pdbx_audit_revision_history.revision_date") == mt.REVISION_DATE
    assert set(block.find_loop("_atom_site.label_asym_id")) == {"A"}
    assert set(block.find_loop("_atom_site.auth_asym_id")) == {"A"}
    assert len(block.find_loop("_entity_poly_seq.mon_id")) == len(t["PDL1"])
    # OpenFold3's CIF-direct parser needs these two
    assert block.find_value("_entity_poly.pdbx_seq_one_letter_code_can") == t["PDL1"]
    assert list(block.find_loop("_pdbx_poly_seq_scheme.asym_id")) == ["A"] * len(t["PDL1"])


@pytest.mark.parametrize("name", ["notes.txt", "model.json"])
def test_non_structure_files_skipped(tmp_path, name):
    junk = tmp_path / name
    junk.write_text("x")
    t = _targets()
    index, _, _ = _run(tmp_path, [junk, PDL1_PDB], [("PDL1", t["PDL1"])])
    assert index[mt.seq_key(t["PDL1"])][0]["source"] == "PDL1.pdb:A"
