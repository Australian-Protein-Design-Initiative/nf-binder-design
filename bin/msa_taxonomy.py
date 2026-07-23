#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Canonical taxonomy annotation for multimer paired MSAs (fold.nf Phase 2).

This is the single place taxonomy logic lives. It parses each a3m hit header
ONCE into a canonical record (sequence, accession, tax_id, species_mnemonic)
and then RENDERS the per-tool file each engine needs, because the three offline
engines key cross-chain pairing on *different* header fields (ground truth,
plans/fold-nf-multimer-paired-msa.md §0):

  - RF3 (atomworks parse_a3m):  re.search(r"TaxID=(\\d+)", header) - numeric taxid.
  - Protenix (MSAPairingEngine.get_species_ids): the species *mnemonic*
      (_HUMAN, _9BETA) via
        _UNIPROT_REGEX = (?:tr|sp)\\|[A-Z0-9]{6,10}(?:_\\d+)?\\|[A-Z0-9]{1,10}_(?P<SpeciesId>[A-Z0-9]{1,5})
        _UNIREF_REGEX  = ^UniRef100_[^_]+_([^_/]+)
      NOT the numeric TaxID=.
  - Boltz (data/parse/csv.py): a `key,sequence` CSV where `key` == taxid;
      empty/nan -> -1 (unpaired). Cross-chain rows sharing a key are paired.

Because RF3 keys on numeric taxid and Protenix on the species mnemonic, no single
a3m header form feeds both at full coverage - hence one canonical parse + N
renderers rather than one universal a3m.

Header coverage note: the `--msa_method jackhmmer_af2` route preserves rich
UniProt/UniRef descriptions (TaxID=, OX=, RepID=ACC_SPECIES, tr|ACC|X_SPECIES),
so taxonomy is recoverable. The `--msa_method mmseqs2_colabfold` route emits bare
`>UniRef100_ACC<TAB>stats` headers with NO taxonomy - those records degrade to
unpaired here, and ColabFold multimer pairing should instead use --use_msa_server
(or the ticket/pair depth booster, deferred). Every render logs the resolved
paired depth per chain so a thin/absent pairing is never silent.
"""

import argparse
import csv
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, TextIO, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


# --- Header-parsing patterns (order = precedence for extraction) ----------------

# Numeric taxonomy id. UniRef descriptions use `TaxID=<n>`; UniProt (sp|/tr|)
# descriptions use `OX=<n>` (organism taxonomy id). Either resolves the numeric id
# RF3/Boltz pair on.
_TAXID_RE = re.compile(r"(?:TaxID|OX)[=:](\d+)")

# Species mnemonic (HUMAN, 9BETA, ...):
#  - sp|ACC|ENTRYNAME_SPECIES  /  tr|ACC|ENTRYNAME_SPECIES
_SPTR_SPECIES_RE = re.compile(r"\b(?:sp|tr)\|[A-Za-z0-9]+\|[A-Za-z0-9]+_([A-Za-z0-9]{1,10})")
#  - RepID=ACC_SPECIES  (UniRef description tail)
_REPID_SPECIES_RE = re.compile(r"RepID=[A-Za-z0-9]+_([A-Za-z0-9]{1,10})")
#  - compact UniRef token `UniRef100_ACC_SPECIES` (acc carries no underscore)
_UNIREF_COMPACT_RE = re.compile(r"^UniRef\d+_[A-Za-z0-9]+_([A-Za-z0-9]{1,10})(?:\s|$)")

# Accession, for synthesising Protenix-compatible compact headers.
_SPTR_ACC_RE = re.compile(r"\b(?:sp|tr)\|([A-Za-z0-9]+)\|")
_UNIREF_ACC_RE = re.compile(r"^UniRef\d+_([A-Za-z0-9]+)")
_REPID_ACC_RE = re.compile(r"RepID=([A-Za-z0-9]+)_")


@dataclass
class MsaRecord:
    """One canonical a3m record. `is_query` marks the first (target) sequence."""

    header: str  # original header text (without the leading '>')
    sequence: str  # aligned sequence, verbatim (a3m may carry lower-case inserts)
    accession: Optional[str] = None
    tax_id: Optional[str] = None
    species: Optional[str] = None
    is_query: bool = False


def _first_token(header: str) -> str:
    # a3m/ColabFold headers separate the id from stats/description by space OR tab.
    return re.split(r"[\s\t]", header.strip(), maxsplit=1)[0]


def _extract_tax_id(header: str) -> Optional[str]:
    m = _TAXID_RE.search(header)
    return m.group(1) if m else None


def _extract_species(header: str) -> Optional[str]:
    for rx in (_SPTR_SPECIES_RE, _REPID_SPECIES_RE, _UNIREF_COMPACT_RE):
        m = rx.search(header)
        if m:
            return m.group(1)
    return None


def _extract_accession(header: str) -> Optional[str]:
    for rx in (_SPTR_ACC_RE, _UNIREF_ACC_RE, _REPID_ACC_RE):
        m = rx.search(header)
        if m:
            return m.group(1)
    # Fall back to the first whitespace/tab-delimited token, minus any UniRefNN_ prefix.
    tok = _first_token(header)
    tok = re.sub(r"^UniRef\d+_", "", tok)
    tok = re.sub(r"[^A-Za-z0-9]", "", tok)
    return tok or None


def annotate_header(header: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Return (accession, tax_id, species_mnemonic) parsed from one header."""
    return _extract_accession(header), _extract_tax_id(header), _extract_species(header)


def parse_a3m(text: str) -> List[MsaRecord]:
    """Parse a3m text into canonical records. First record is flagged is_query."""
    records: List[MsaRecord] = []
    header: Optional[str] = None
    seq_parts: List[str] = []

    def flush() -> None:
        if header is None:
            return
        seq = "".join(seq_parts)
        acc, tax_id, species = annotate_header(header)
        records.append(
            MsaRecord(
                header=header,
                sequence=seq,
                accession=acc,
                tax_id=tax_id,
                species=species,
                is_query=(len(records) == 0),
            )
        )

    for line in text.splitlines():
        if line.startswith(">"):
            flush()
            header = line[1:].rstrip("\n")
            seq_parts = []
        elif header is not None:
            seq_parts.append(line.strip())
    flush()
    return records


# --- Renderers ------------------------------------------------------------------


def render_rf3_a3m(records: List[MsaRecord]) -> str:
    """a3m whose hit headers carry `TaxID=<n>` so atomworks pairs by numeric taxid.

    The query (row 0) is kept verbatim and first. Hits with a resolved tax_id get a
    canonical `TaxID=<n>` token appended if one is not already present (atomworks
    only needs the token to exist in the header line). Hits with no tax_id keep their
    original header and remain in the file as unpaired-tail rows.
    """
    lines: List[str] = []
    for rec in records:
        if rec.is_query:
            lines += [f">{rec.header}", rec.sequence]
            continue
        header = rec.header
        if rec.tax_id and not re.search(r"TaxID=\d+", header):
            header = f"{header} TaxID={rec.tax_id}"
        lines += [f">{header}", rec.sequence]
    return "\n".join(lines) + "\n"


def render_protenix_paired_a3m(records: List[MsaRecord]) -> str:
    """a3m whose hit headers are Protenix-parseable compact `UniRef100_ACC_SPECIES`.

    Only records with a resolved species mnemonic are emitted as hits (Protenix pairs
    on species; species-less rows contribute nothing to pairing). The query is kept
    first, verbatim.
    """
    lines: List[str] = []
    for i, rec in enumerate(records):
        if rec.is_query:
            lines += [f">{rec.header}", rec.sequence]
            continue
        if not rec.species:
            continue
        acc = (rec.accession or f"seq{i}").replace("_", "")
        lines += [f">UniRef100_{acc}_{rec.species}", rec.sequence]
    return "\n".join(lines) + "\n"


def render_protenix_unpaired_a3m(records: List[MsaRecord]) -> str:
    """The plain merged a3m (all records, original headers) for `unpairedMsaPath`."""
    lines: List[str] = []
    for rec in records:
        lines += [f">{rec.header}", rec.sequence]
    return "\n".join(lines) + "\n"


def render_boltz_csv(records: List[MsaRecord], out: TextIO) -> None:
    """Boltz `key,sequence` CSV; key=<tax_id> (blank -> Boltz treats as unpaired).

    The query row (row 0) is emitted with a blank key: Boltz always treats each
    chain MSA's first row as the query and pairs the queries positionally.
    """
    w = csv.writer(out)
    w.writerow(["key", "sequence"])
    for rec in records:
        key = "" if rec.is_query else (rec.tax_id or "")
        w.writerow([key, rec.sequence])


def _count_paired(records: List[MsaRecord], attr: str) -> int:
    """Count non-query hits with a resolvable pairing key for a given attribute."""
    return sum(1 for r in records if not r.is_query and getattr(r, attr))


# --- CLI ------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--a3m", required=True, help="Input per-chain a3m (query first).")
    parser.add_argument(
        "--tool", required=True, choices=["rf3", "protenix", "boltz"], help="Target engine."
    )
    parser.add_argument("--out", help="Output path (rf3 a3m, or boltz csv).")
    parser.add_argument("--paired-out", help="Protenix pairedMsaPath a3m output.")
    parser.add_argument("--unpaired-out", help="Protenix unpairedMsaPath a3m output.")
    parser.add_argument(
        "--chain-id", default="?", help="Chain id, for the paired-depth log line only."
    )
    args = parser.parse_args()

    records = parse_a3m(Path(args.a3m).read_text())
    if not records:
        log.error("No records parsed from %s", args.a3m)
        return 1

    if args.tool == "rf3":
        if not args.out:
            parser.error("--tool rf3 needs --out")
        Path(args.out).write_text(render_rf3_a3m(records))
        depth = _count_paired(records, "tax_id")
        log.info("chain %s: RF3 a3m %d rows, %d with TaxID= (pairable)", args.chain_id, len(records), depth)
    elif args.tool == "protenix":
        if not (args.paired_out and args.unpaired_out):
            parser.error("--tool protenix needs --paired-out and --unpaired-out")
        Path(args.paired_out).write_text(render_protenix_paired_a3m(records))
        Path(args.unpaired_out).write_text(render_protenix_unpaired_a3m(records))
        depth = _count_paired(records, "species")
        log.info(
            "chain %s: Protenix unpaired %d rows, paired %d species-tagged rows",
            args.chain_id, len(records), depth,
        )
    elif args.tool == "boltz":
        if not args.out:
            parser.error("--tool boltz needs --out")
        with open(args.out, "w", newline="") as f:
            render_boltz_csv(records, f)
        depth = _count_paired(records, "tax_id")
        log.info("chain %s: Boltz csv %d rows, %d keyed by taxid (pairable)", args.chain_id, len(records), depth)

    return 0


if __name__ == "__main__":
    sys.exit(main())
