#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["gemmi", "biopython"]
# ///

"""
Match user-supplied template structures to the chains being folded (fold.nf --templates).

Every protein chain in every template file (.pdb/.cif, single chains or whole
complexes) is aligned against every unique query sequence (local alignment,
BLOSUM62), much like AF2's template search with the user's files as the
database. A template chain is accepted for a query sequence when it passes
--min-identity (over aligned positions), --min-coverage (fraction of the query
covered) and --min-aligned (aligned residues, relaxed to the query length for
queries shorter than that); the best --max-per-chain by identity x coverage are
kept.

Outputs, in --outdir:
  - tmplNNN.cif: one normalised single-chain mmCIF per accepted template chain
    (chain A, resolved residues only, renumbered 1..n, with _entity_poly_seq and
    a _pdbx_audit_revision_history date, which AF3 and others require).
  - index.json: {md5(query sequence): [{cif, template_chain, source, identity,
    coverage, query_indices, template_indices}, ...]}, best first. Indices are
    0-based and index tmplNNN.cif's residues, which are all resolved.
  - templates_matched.tsv: every query x template-chain pair considered.
Engine input generators look each chain up by the md5 of its sequence.
"""

import argparse
import hashlib
import json
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import gemmi
from Bio import Align
from Bio.Align import substitution_matrices

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

TEMPLATE_SUFFIXES = (".pdb", ".ent", ".cif", ".mmcif", ".pdb.gz", ".ent.gz", ".cif.gz", ".mmcif.gz")
MIN_CHAIN_LENGTH = 10
# A fixed, old date: engines that filter templates by release date must keep these.
REVISION_DATE = "2000-01-01"


def seq_key(seq: str) -> str:
    return hashlib.md5(seq.upper().encode()).hexdigest()


@dataclass
class TemplateChain:
    source: str
    chain_id: str
    residues: List[gemmi.Residue] = field(default_factory=list)
    sequence: str = ""


@dataclass
class Match:
    query_key: str
    query_ids: List[str]
    template: TemplateChain
    identity: float
    coverage: float
    query_indices: List[int]
    template_indices: List[int]

    @property
    def score(self) -> float:
        return self.identity * self.coverage


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header = line[1:].split()[0] if line[1:].split() else ""
            seq = []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def is_template_file(path: Path) -> bool:
    return path.name.lower().endswith(TEMPLATE_SUFFIXES)


def read_template_chains(path: Path) -> List[TemplateChain]:
    """Protein chains of the first model, keeping amino-acid residues that have a CA."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    if len(st) == 0:
        log.warning("%s: no models - skipped", path.name)
        return []
    chains: List[TemplateChain] = []
    for chain in st[0]:
        tc = TemplateChain(source=path.name, chain_id=chain.name)
        letters: List[str] = []
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid() or res.find_atom("CA", "*") is None:
                continue
            code = info.one_letter_code.upper()
            letters.append(code if code.isalpha() and code != " " else "X")
            tc.residues.append(res)
        tc.sequence = "".join(letters)
        if len(tc.sequence) >= MIN_CHAIN_LENGTH:
            chains.append(tc)
    if not chains:
        log.warning("%s: no protein chains of >= %d resolved residues", path.name, MIN_CHAIN_LENGTH)
    return chains


def make_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def align(aligner: Align.PairwiseAligner, query: str, template: str) -> Tuple[float, float, List[int], List[int]]:
    """identity over aligned pairs, query coverage, and the 0-based index mapping."""
    q = query.upper()
    alignment = aligner.align(q, template.upper())[0]
    q_idx: List[int] = []
    t_idx: List[int] = []
    for (qs, qe), (ts, te) in zip(*alignment.aligned):
        q_idx.extend(range(qs, qe))
        t_idx.extend(range(ts, te))
    if not q_idx:
        return 0.0, 0.0, [], []
    matches = sum(1 for i, j in zip(q_idx, t_idx) if q[i] == template[j].upper())
    return matches / len(q_idx), len(q_idx) / len(query), q_idx, t_idx


def write_template_cif(tc: TemplateChain, name: str, out_path: Path) -> None:
    """Single chain 'A', resolved residues only, renumbered, with entity sequence and a revision date."""
    st = gemmi.Structure()
    st.name = name
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for i, res in enumerate(tc.residues, start=1):
        new = gemmi.Residue()
        new.name = res.name
        new.seqid = gemmi.SeqId(i, " ")
        new.entity_type = gemmi.EntityType.Polymer
        new.het_flag = "A"
        # Boltz names chains by label_asym_id, so keep it equal to the auth chain id
        # (gemmi would otherwise assign "Axp").
        new.subchain = "A"
        for atom in res:
            if atom.altloc not in ("\0", "A"):
                continue
            a = atom.clone()
            a.altloc = "\0"
            new.add_atom(a)
        chain.add_residue(new)
    model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    for ent in st.entities:
        if ent.entity_type == gemmi.EntityType.Polymer:
            ent.name = "1"
            ent.full_sequence = [r.name for r in tc.residues]
            ent.polymer_type = gemmi.PolymerType.PeptideL
    st.assign_label_seq_id()
    doc = st.make_mmcif_document()
    block = doc.sole_block()
    block.name = name
    # OpenFold3's CIF-direct parser reads the canonical one-letter sequence, which
    # gemmi does not write.
    block.set_mmcif_category("_entity_poly.", {
        "entity_id": ["1"],
        "type": ["polypeptide(L)"],
        "pdbx_strand_id": ["A"],
        "pdbx_seq_one_letter_code": [tc.sequence],
        "pdbx_seq_one_letter_code_can": [tc.sequence],
    })
    n = len(tc.residues)
    nums = [str(i) for i in range(1, n + 1)]
    names = [r.name for r in tc.residues]
    block.set_mmcif_category("_pdbx_poly_seq_scheme.", {
        "asym_id": ["A"] * n, "entity_id": ["1"] * n, "seq_id": nums, "mon_id": names,
        "ndb_seq_num": nums, "pdb_seq_num": nums, "auth_seq_num": nums,
        "pdb_mon_id": names, "auth_mon_id": names, "pdb_strand_id": ["A"] * n,
        "pdb_ins_code": ["."] * n, "hetero": ["n"] * n,
    })
    chem_comp = block.find_mmcif_category("_chem_comp.")
    for row in chem_comp:
        row[1] = "'PEPTIDE LINKING'" if row[0] == "GLY" else "'L-PEPTIDE LINKING'"
    loop = block.init_loop("_pdbx_audit_revision_history.", ["ordinal", "data_content_type", "major_revision", "minor_revision", "revision_date"])
    loop.add_row(["1", "'Structure model'", "1", "0", REVISION_DATE])
    doc.write_file(str(out_path))


def match_templates(
    template_paths: List[Path],
    query_records: List[Tuple[str, str]],
    min_identity: float,
    min_coverage: float,
    min_aligned: int,
    max_per_chain: int,
) -> Tuple[Dict[str, List[Match]], List[List[str]]]:
    queries: Dict[str, Tuple[str, List[str]]] = {}
    for qid, seq in query_records:
        if not seq:
            continue
        key = seq_key(seq)
        queries.setdefault(key, (seq.upper(), []))[1].append(qid)

    template_chains: List[TemplateChain] = []
    for p in sorted(template_paths, key=lambda x: x.name):
        if not is_template_file(p):
            log.warning("%s: not a .pdb/.cif file - skipped", p.name)
            continue
        template_chains.extend(read_template_chains(p))
    log.info("%d template chain(s) from %d file(s); %d unique query sequence(s)",
             len(template_chains), len(template_paths), len(queries))

    aligner = make_aligner()
    accepted: Dict[str, List[Match]] = {}
    report: List[List[str]] = []
    for key, (seq, ids) in queries.items():
        # Queries shorter than min_aligned must be covered end to end, so a short
        # peptide can still be templated without opening the door to fragments.
        min_len = min(min_aligned, len(seq))
        candidates: List[Match] = []
        for tc in template_chains:
            identity, coverage, q_idx, t_idx = align(aligner, seq, tc.sequence)
            m = Match(key, ids, tc, identity, coverage, q_idx, t_idx)
            if identity < min_identity:
                report.append(_report_row(m, "rejected", f"identity < {min_identity}"))
            elif coverage < min_coverage:
                report.append(_report_row(m, "rejected", f"coverage < {min_coverage}"))
            elif len(q_idx) < min_len:
                report.append(_report_row(m, "rejected", f"aligned < {min_len}"))
            else:
                candidates.append(m)
        candidates.sort(key=lambda m: m.score, reverse=True)
        for rank, m in enumerate(candidates, start=1):
            if rank <= max_per_chain:
                accepted.setdefault(key, []).append(m)
                report.append(_report_row(m, "accepted", "", rank))
            else:
                report.append(_report_row(m, "rejected", f"rank > {max_per_chain}", rank))
        if key not in accepted:
            log.warning("No template matched %s", ",".join(ids))
    return accepted, report


def _report_row(m: Match, status: str, reason: str, rank: Optional[int] = None) -> List[str]:
    return [
        ",".join(m.query_ids), f"{m.template.source}:{m.template.chain_id}",
        f"{m.identity:.3f}", f"{m.coverage:.3f}", str(len(m.query_indices)),
        str(rank) if rank else "", status, reason,
    ]


def write_outputs(accepted: Dict[str, List[Match]], report: List[List[str]], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    cif_names: Dict[Tuple[str, str], str] = {}
    index: Dict[str, List[dict]] = {}
    for key, matches in accepted.items():
        for m in matches:
            tkey = (m.template.source, m.template.chain_id)
            if tkey not in cif_names:
                # No underscore: OpenFold3 splits template ids on "_".
                cif_names[tkey] = f"tmpl{len(cif_names) + 1:03d}.cif"
                write_template_cif(m.template, cif_names[tkey][:-4], outdir / cif_names[tkey])
            index.setdefault(key, []).append({
                "cif": cif_names[tkey],
                "template_chain": "A",
                "source": f"{m.template.source}:{m.template.chain_id}",
                "identity": round(m.identity, 4),
                "coverage": round(m.coverage, 4),
                "query_indices": m.query_indices,
                "template_indices": m.template_indices,
            })
    (outdir / "index.json").write_text(json.dumps(index))
    header = ["query_ids", "template", "identity", "coverage", "aligned", "rank", "status", "reason"]
    lines = ["\t".join(header)] + ["\t".join(r) for r in report]
    (outdir / "templates_matched.tsv").write_text("\n".join(lines) + "\n")
    log.info("Wrote %d template mmCIF(s) for %d query sequence(s) to %s", len(cif_names), len(index), outdir)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--templates", nargs="+", required=True, help="Template structure files (.pdb/.cif)")
    parser.add_argument("--queries", nargs="+", required=True, help="FASTA file(s) of the chains that may be templated")
    parser.add_argument("--min-identity", type=float, default=0.3, help="Minimum identity over aligned positions (default: 0.3)")
    parser.add_argument("--min-coverage", type=float, default=0.3, help="Minimum fraction of the query covered (default: 0.3)")
    parser.add_argument("--min-aligned", type=int, default=40,
                        help="Minimum aligned residues, or the query length if shorter (default: 40)")
    parser.add_argument("--max-per-chain", type=int, default=4, help="Templates kept per query sequence (default: 4)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    records: List[Tuple[str, str]] = []
    for q in args.queries:
        records.extend(parse_fasta(Path(q)))
    accepted, report = match_templates(
        [Path(p) for p in args.templates], records,
        args.min_identity, args.min_coverage, args.min_aligned, args.max_per_chain,
    )
    write_outputs(accepted, report, Path(args.outdir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
