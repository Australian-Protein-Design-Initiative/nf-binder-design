#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Pure-python unit tests for bin/msa_taxonomy.py - the single source of taxonomy
truth for multimer paired MSAs (fold.nf Phase 2).

Run from the repo root with uv (host python has no pytest):

    uv run --with pytest pytest tests/bin/test_msa_taxonomy.py -q

The renderer assertions are checked against the ACTUAL regexes each downstream
engine uses to pair chains (copied verbatim below from the container ground
truth in plans/fold-nf-multimer-paired-msa.md §0), so a rendering regression
that would silently break pairing fails here.
"""

import importlib.util
import io
import re
from pathlib import Path

import pytest

# Load bin/msa_taxonomy.py by path (it is a standalone PEP-723 script, not a package).
_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "msa_taxonomy.py"
_spec = importlib.util.spec_from_file_location("msa_taxonomy", _MODPATH)
mt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mt)


# --- Downstream engines' real pairing regexes (ground truth, §0) ----------------

# RF3 / atomworks _msa_loading_utils.extract_tax_id
RF3_TAXID_RE = re.compile(r"TaxID=(\d+)")
# Protenix msa_utils.MSAPairingEngine.get_species_ids
PROTENIX_UNIPROT_RE = re.compile(
    r"(?:tr|sp)\|[A-Z0-9]{6,10}(?:_\d+)?\|[A-Z0-9]{1,10}_(?P<SpeciesId>[A-Z0-9]{1,5})"
)
PROTENIX_UNIREF_RE = re.compile(r"^UniRef100_[^_]+_([^_/]+)")


# --- Realistic header samples ---------------------------------------------------

QUERY_HEADER = "101"
QUERY_SEQ = "MKTAYIAKQR"

# jackhmmer_af2 route: rich UniProt/UniRef descriptions (taxonomy recoverable).
UNIPROT_TR = "tr|Q8QRZ0|Q8QRZ0_9BETA Membrane glycoprotein UL119 OS=Panine betaherpesvirus 2 OX=1608321 GN=UL119 PE=4 SV=1"
UNIREF_FULL = "UniRef100_A0A7T7DGJ2 Ba152/Ba151 n=1 Tax=Baboon cytomegalovirus TaxID=120505 RepID=A0A7T7DGJ2_9BETA"
UNIPROT_SP = "sp|P12345|XPROT_HUMAN Example protein OS=Homo sapiens OX=9606 GN=X PE=1 SV=2"
# mmseqs2_colabfold route: bare UniRef100 id + tab-separated stats (NO taxonomy).
COLABFOLD_BARE = "UniRef100_A4GW17\t197\t0.959\t5.349E-53\t0\t221\t222\t17\t238\t290"


def build_a3m(*hits):
    lines = [f">{QUERY_HEADER}", QUERY_SEQ]
    for i, h in enumerate(hits):
        lines += [f">{h}", "MKTAYIAKQ" + "ACDEFG"[i % 6]]
    return "\n".join(lines) + "\n"


# --- annotate_header ------------------------------------------------------------


@pytest.mark.parametrize(
    "header,acc,tax,species",
    [
        (UNIPROT_TR, "Q8QRZ0", "1608321", "9BETA"),      # OX= -> taxid; tr| -> mnemonic
        (UNIREF_FULL, "A0A7T7DGJ2", "120505", "9BETA"),  # TaxID= ; RepID= -> mnemonic
        (UNIPROT_SP, "P12345", "9606", "HUMAN"),         # OX= ; sp| -> mnemonic
        (COLABFOLD_BARE, "A4GW17", None, None),          # bare colabfold: no taxonomy
    ],
)
def test_annotate_header(header, acc, tax, species):
    got_acc, got_tax, got_species = mt.annotate_header(header)
    assert got_acc == acc
    assert got_tax == tax
    assert got_species == species


def test_parse_a3m_flags_query_and_keeps_sequences():
    recs = mt.parse_a3m(build_a3m(UNIPROT_TR, COLABFOLD_BARE))
    assert len(recs) == 3
    assert recs[0].is_query and recs[0].sequence == QUERY_SEQ
    assert not recs[1].is_query
    assert recs[1].tax_id == "1608321" and recs[1].species == "9BETA"
    assert recs[2].tax_id is None and recs[2].species is None


# --- RF3 renderer: every taxonomic hit must carry a TaxID= token ----------------


def test_render_rf3_appends_taxid_and_rf3_can_parse_it():
    recs = mt.parse_a3m(build_a3m(UNIPROT_TR, UNIREF_FULL, UNIPROT_SP, COLABFOLD_BARE))
    out = mt.render_rf3_a3m(recs)
    headers = [ln[1:] for ln in out.splitlines() if ln.startswith(">")]
    # query header untouched
    assert headers[0] == QUERY_HEADER
    # tr| hit (taxid from OX=) gets an appended TaxID= that atomworks will parse
    assert RF3_TAXID_RE.search(headers[1]).group(1) == "1608321"
    # UniRef hit already had TaxID=120505 - not duplicated
    assert headers[2].count("TaxID=") == 1
    assert RF3_TAXID_RE.search(headers[2]).group(1) == "120505"
    # sp| hit (OX=9606)
    assert RF3_TAXID_RE.search(headers[3]).group(1) == "9606"
    # bare ColabFold hit stays taxonomy-less (unpaired tail)
    assert "TaxID=" not in headers[4]


def test_render_rf3_preserves_all_rows():
    recs = mt.parse_a3m(build_a3m(UNIPROT_TR, COLABFOLD_BARE))
    out = mt.render_rf3_a3m(recs)
    assert sum(1 for ln in out.splitlines() if ln.startswith(">")) == 3


# --- Protenix renderer: paired headers must match Protenix's species regexes ----


def test_render_protenix_paired_matches_protenix_uniref_regex():
    recs = mt.parse_a3m(build_a3m(UNIPROT_TR, UNIREF_FULL, UNIPROT_SP, COLABFOLD_BARE))
    paired = mt.render_protenix_paired_a3m(recs)
    headers = [ln[1:] for ln in paired.splitlines() if ln.startswith(">")]
    # query first, then only species-bearing hits (3 of 4; bare ColabFold dropped)
    assert headers[0] == QUERY_HEADER
    hit_headers = headers[1:]
    assert len(hit_headers) == 3
    species = [PROTENIX_UNIREF_RE.match(h).group(1) for h in hit_headers]
    assert species == ["9BETA", "9BETA", "HUMAN"]


def test_render_protenix_unpaired_keeps_everything_verbatim():
    recs = mt.parse_a3m(build_a3m(UNIPROT_TR, COLABFOLD_BARE))
    unpaired = mt.render_protenix_unpaired_a3m(recs)
    headers = [ln[1:] for ln in unpaired.splitlines() if ln.startswith(">")]
    assert headers == [QUERY_HEADER, UNIPROT_TR, COLABFOLD_BARE]


def test_protenix_paired_accession_has_no_underscore():
    # _UNIREF_REGEX = ^UniRef100_[^_]+_(species): the accession segment must not
    # contain an underscore or the species group would capture the wrong token.
    recs = mt.parse_a3m(build_a3m(UNIREF_FULL))
    paired = mt.render_protenix_paired_a3m(recs)
    hit = [ln[1:] for ln in paired.splitlines() if ln.startswith(">")][1]
    acc_seg = hit.split("_")[1]  # UniRef100_<acc>_<species>
    assert "_" not in acc_seg


# --- Boltz renderer: key,sequence CSV keyed on taxid ----------------------------


def test_render_boltz_csv_keys_on_taxid():
    recs = mt.parse_a3m(build_a3m(UNIPROT_TR, UNIREF_FULL, COLABFOLD_BARE))
    buf = io.StringIO()
    mt.render_boltz_csv(recs, buf)
    rows = [r for r in buf.getvalue().splitlines()]
    assert rows[0] == "key,sequence"
    # query row: blank key
    assert rows[1].startswith(",")
    # taxonomic hits carry their numeric key
    assert rows[2].startswith("1608321,")
    assert rows[3].startswith("120505,")
    # bare ColabFold hit: blank key (Boltz -> unpaired)
    assert rows[4].startswith(",")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
