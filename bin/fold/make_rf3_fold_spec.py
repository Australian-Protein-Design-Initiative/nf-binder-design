#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Generic RF3 (RosettaFold3) input JSON for the standalone fold.nf
ROSETTAFOLD3_FOLD subworkflow.

NOT to be confused with bin/rfd3/make_rf3_input_spec.py, which is hard-wired
to a 2-body target+binder complex for the binder-design pipeline (rfd3.nf) and
always emits target+binder roles plus optional template handling.

Builds one {seq, chain_id, msa_path?} component per input FASTA record (chain
IDs A, B, C, ... in file order), with no target/binder roles, no template, and
no binder postprocessing - this is exactly what rf3's generic
InferenceInput.from_json_dict() / components_to_atom_array() accepts (ground
truth confirmed by inspecting rc-foundry:0.2.0-weights' rf3/utils/inference.py
on 2026-07-17: it builds an atom array generically from `components`, with no
target/binder-specific logic at all).

Multimer: pass one --a3m per chain in record order (chain A, B, C, ...); each
becomes that component's msa_path. RF3 (atomworks) pairs chains internally by
matching TaxID=<n> parsed from the a3m hit headers, so the per-chain a3m must
be TaxID=-annotated (bin/fold/msa_taxonomy.py --tool rf3 does this). A single --a3m
with a single-record FASTA is the monomer case (unchanged).

Templates: RF3 has no separate template input - a templated chain is given as a
structure component ({"path": ...}) and named in template_selection, so the
file's sequence becomes the chain's sequence. With --templates-dir
(bin/fold/match_templates.py output), each chain in --template-chains (default:
all) with a match gets an "in-place" mmCIF built from its best template: the
FASTA sequence, with template coordinates on aligned residues (backbone only
where the residue differs) and no atoms elsewhere; atomworks then fills those
in as unresolved from the entity sequence. Its MSA moves to the
top-level msa_paths. RF3 uses one template per chain.
"""

import argparse
import hashlib
import json
import logging
import string
import sys
from pathlib import Path
from typing import List, Optional

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)


def parse_fasta_records(fasta_path: Path) -> List[str]:
    """Return one sequence string per FASTA record, in file order (no external deps)."""
    records: List[str] = []
    current: List[str] = []
    for line in fasta_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                records.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        records.append("".join(current))
    return records


def a3m_query(a3m_path: Path) -> Optional[str]:
    """First sequence in an a3m (skipping ColabFold '#' lines), gaps removed."""
    seq: List[str] = []
    seen_header = False
    for raw in a3m_path.read_text().replace("\x00", "").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if seen_header:
                break
            seen_header = True
        elif seen_header:
            seq.append(line)
    return "".join(seq).replace("-", "").upper() if seen_header else None


def _check_query(path: Path, seq: str, chain_index: int) -> None:
    # A mis-ordered bundle otherwise folds each chain with another chain's MSA
    # without any error from RF3 (this is the bug class fixed in commit 1ae37a8).
    query = a3m_query(path)
    if query is not None and query != seq.upper():
        raise ValueError(
            f"chain {chain_index}: first sequence in {path.name} does not match the FASTA "
            f"record - MSA files are not in chain order"
        )


TEMPLATE_DIR_NAME = "rf3_templates"
BACKBONE = ("N", "CA", "C", "O")


def best_template(templates_dir: Optional[Path], seq: str) -> Optional[dict]:
    if templates_dir is None or not (templates_dir / "index.json").exists():
        return None
    index = json.loads((templates_dir / "index.json").read_text())
    hits = index.get(hashlib.md5(seq.upper().encode()).hexdigest(), [])
    return hits[0] if hits else None


def write_inplace_cif(template_cif: Path, hit: dict, seq: str, chain_id: str, out_path: Path) -> None:
    """The query chain with the template's coordinates on its aligned residues."""
    import gemmi

    tmpl = gemmi.read_structure(str(template_cif))[0][hit["template_chain"]]
    mapping = dict(zip(hit["query_indices"], hit["template_indices"]))
    names = [gemmi.expand_one_letter(c, gemmi.ResidueKind.AA) or "UNK" for c in seq.upper()]

    st = gemmi.Structure()
    st.name = out_path.stem
    model = gemmi.Model("1")
    chain = gemmi.Chain(chain_id)
    for qi, name in enumerate(names):
        if qi not in mapping:
            continue
        src = tmpl[mapping[qi]]
        res = gemmi.Residue()
        res.name = name
        res.seqid = gemmi.SeqId(qi + 1, " ")
        res.entity_type = gemmi.EntityType.Polymer
        res.het_flag = "A"
        res.subchain = chain_id
        same = src.name == name
        for atom in src:
            if same or atom.name in BACKBONE:
                res.add_atom(atom.clone())
        chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    for ent in st.entities:
        ent.name = "1"
        ent.full_sequence = names
        ent.polymer_type = gemmi.PolymerType.PeptideL
    st.assign_label_seq_id()
    doc = st.make_mmcif_document()
    doc.sole_block().name = out_path.stem
    doc.write_file(str(out_path))


def make_rf3_fold_spec(
    fasta_path: Path,
    name: str,
    a3m_paths: Optional[List[Path]] = None,
    out_dir: Path = Path("."),
    templates_dir: Optional[Path] = None,
    template_chains: Optional[List[str]] = None,
) -> dict:
    sequences = parse_fasta_records(fasta_path)
    if not sequences:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    chain_ids = list(string.ascii_uppercase)
    if len(sequences) > len(chain_ids):
        raise ValueError(f"{fasta_path} has more than {len(chain_ids)} chains; not supported")

    # a3m paths are matched to chains by position (record order): exactly one
    # per chain (a monomer is just the n==1 case). bin/fold/msa_taxonomy.py renders
    # the per-chain TaxID=-annotated a3m.
    if a3m_paths and len(a3m_paths) != len(sequences):
        raise ValueError(
            f"--a3m got {len(a3m_paths)} file(s) for {len(sequences)} chain(s); "
            f"pass exactly one a3m per chain in record order"
        )

    (out_dir / TEMPLATE_DIR_NAME).mkdir(parents=True, exist_ok=True)
    components = []
    templated: List[str] = []
    msa_paths = {}
    for i, (chain_id, seq) in enumerate(zip(chain_ids, sequences)):
        comp = {"seq": seq, "chain_id": chain_id}
        if a3m_paths:
            _check_query(a3m_paths[i], seq, i)
            # Basename so fold.nf can overwrite the staged a3m in-task
            # (MSA subsample) without rewriting this JSON.
            comp["msa_path"] = a3m_paths[i].name
        hit = best_template(templates_dir, seq) if (not template_chains or chain_id in template_chains) else None
        if hit is not None:
            rel = f"{TEMPLATE_DIR_NAME}/chain_{chain_id}.cif"
            write_inplace_cif(templates_dir / hit["cif"], hit, seq, chain_id, out_dir / rel)
            if "msa_path" in comp:
                msa_paths[chain_id] = comp["msa_path"]
            # add_missing_atoms keeps residues the template lacks (as unresolved,
            # so untemplated) rather than dropping them from the chain.
            comp = {"path": rel, "custom_parse_kwargs": {"add_missing_atoms": True}}
            templated.append(chain_id)
            log.info("chain %s: template %s", chain_id, hit["source"])
        components.append(comp)

    spec = {"name": name, "components": components}
    if templated:
        spec["template_selection"] = templated
    if msa_paths:
        spec["msa_paths"] = msa_paths
    return spec


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="FASTA file (one record per chain)")
    parser.add_argument("--name", required=True, help="RF3 'name' field")
    parser.add_argument(
        "--a3m",
        nargs="+",
        default=None,
        help="Optional a3m MSA(s): one per chain in record order (or a single a3m for a monomer)",
    )
    parser.add_argument("--templates-dir", dest="templates_dir", default=None, help="Matched-templates directory (bin/fold/match_templates.py)")
    parser.add_argument("--template-chains", dest="template_chains", nargs="+", default=None, help="Chain ID(s) that may be templated (default: all)")
    parser.add_argument("-o", "--output", required=True, help="Output JSON path")
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    spec = make_rf3_fold_spec(
        fasta_path=Path(args.fasta),
        name=args.name,
        a3m_paths=[Path(p) for p in args.a3m] if args.a3m else None,
        out_dir=out_path.parent,
        templates_dir=Path(args.templates_dir) if args.templates_dir else None,
        template_chains=args.template_chains,
    )
    # rf3 fold's `inputs=` argument is a JSON list of examples (see
    # bin/rfd3/make_rf3_input_spec.py's json-batch mode for the same
    # convention), even for a single example.
    with open(out_path, "w") as f:
        json.dump([spec], f, indent=2)
    log.info("Wrote RF3 fold input JSON (1 example, %d chain(s)) to %s", len(spec["components"]), out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
