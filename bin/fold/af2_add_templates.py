#!/usr/bin/env python
# /// script
# requires-python = ">=3.8"
# dependencies = ["numpy"]
# ///

"""
Add matched user templates (bin/fold/match_templates.py) to an AF2 features.pkl,
ahead of any templates AF2's own search found.

AF2 runs with --use_precomputed_msas and loads features.pkl directly, so templates
have to be in it. Features are built with AF2's own
templates._extract_template_features from our explicit residue mapping. That
bypasses HhsearchHitFeaturizer, whose duplicate filter would reject a template
identical to the query (the usual case for a known target structure).

Handles both layouts:
  - monomer features (monomer presets, af2_mono chain-break mode, subsampled
    features): one-hot template_aatype, per-template names and sum_probs. With
    --chain-layout (af2_mono), chains are the concatenated FASTA records.
  - multimer features: integer template_aatype in AF2's own order, template_all_atom_mask,
    chains given by asym_id.
User templates go first, and at most 4 templates are kept per row set, the most
any AF2 model uses. Monomer models 1 and 2 use templates; models 3-5 ignore them.

Must run inside the AF2 container (/app/alphafold on sys.path).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import pickle
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

sys.path.insert(0, "/app/alphafold")

import numpy as np  # noqa: E402
from alphafold.common import residue_constants  # noqa: E402
from alphafold.data import mmcif_parsing, templates  # noqa: E402

MAX_TEMPLATES = 4
KALIGN = "/usr/bin/kalign"


def parse_fasta(path: Path) -> List[str]:
    seqs: List[str] = []
    cur: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line)
    if cur:
        seqs.append("".join(cur))
    return seqs


def chain_ranges(feats: dict, seqs: List[str], multimer: bool, layout: Optional[Path]) -> List[Tuple[int, int]]:
    if multimer:
        asym = np.asarray(feats["asym_id"]).astype(int)
        return [(int(np.argmax(asym == a)), int(np.argmax(asym == a)) + int((asym == a).sum())) for a in range(1, len(seqs) + 1)]
    if layout is not None:
        lengths = [c["length"] for c in json.loads(layout.read_text())["chains"]]
    else:
        lengths = [len(s) for s in seqs]
    starts = np.cumsum([0] + lengths[:-1])
    return [(int(s), int(s) + n) for s, n in zip(starts, lengths)]


def chain_hits(templates_dir: Path, seq: str) -> List[dict]:
    index = json.loads((templates_dir / "index.json").read_text())
    return index.get(hashlib.md5(seq.upper().encode()).hexdigest(), [])[:MAX_TEMPLATES]


def hit_features(templates_dir: Path, hit: dict, seq: str) -> Dict[str, np.ndarray]:
    path = templates_dir / hit["cif"]
    file_id = path.stem
    parsed = mmcif_parsing.parse(file_id=file_id, mmcif_string=path.read_text())
    if parsed.mmcif_object is None:
        raise ValueError(f"{path.name}: could not parse: {parsed.errors}")
    tmpl_seq = parsed.mmcif_object.chain_to_seqres[hit["template_chain"]]
    feats, warning = templates._extract_template_features(
        mmcif_object=parsed.mmcif_object,
        pdb_id=file_id,
        mapping=dict(zip(hit["query_indices"], hit["template_indices"])),
        template_sequence=tmpl_seq,
        query_sequence=seq,
        template_chain_id=hit["template_chain"],
        kalign_binary_path=KALIGN,
    )
    if warning:
        log.warning("%s: %s", hit["source"], warning)
    return feats


def add_templates(feats: dict, seqs: List[str], templates_dir: Path,
                  template_chains: Optional[List[str]], layout: Optional[Path]) -> int:
    multimer = "template_all_atom_mask" in feats
    n_res = feats["template_all_atom_positions"].shape[1]
    ranges = chain_ranges(feats, seqs, multimer, layout)
    chain_ids = [chr(ord("A") + i) for i in range(len(seqs))]

    per_chain = {}
    for cid, seq, rng in zip(chain_ids, seqs, ranges):
        if template_chains and cid not in template_chains:
            continue
        hits = chain_hits(templates_dir, seq)
        if hits:
            per_chain[cid] = (rng, [hit_features(templates_dir, h, seq) for h in hits], hits)
            log.info("chain %s: %d template(s): %s", cid, len(hits), ", ".join(h["source"] for h in hits))
    n_new = max((len(v[1]) for v in per_chain.values()), default=0)
    if n_new == 0:
        return 0

    gap = residue_constants.HHBLITS_AA_TO_ID["-"]
    aatype = np.zeros((n_new, n_res, 22), np.float32)
    aatype[..., gap] = 1.0
    masks = np.zeros((n_new, n_res, residue_constants.atom_type_num), np.float32)
    positions = np.zeros((n_new, n_res, residue_constants.atom_type_num, 3), np.float32)
    names: List[bytes] = []
    for t in range(n_new):
        names.append(";".join(h[t]["source"] for _, _, h in per_chain.values() if t < len(h)).encode())
        for (start, end), hit_feats, _ in per_chain.values():
            if t < len(hit_feats):
                aatype[t, start:end] = hit_feats[t]["template_aatype"]
                masks[t, start:end] = hit_feats[t]["template_all_atom_masks"]
                positions[t, start:end] = hit_feats[t]["template_all_atom_positions"]

    def prepend(key: str, new: np.ndarray) -> None:
        old = feats.get(key)
        merged = new if old is None or len(old) == 0 else np.concatenate([new, old.astype(new.dtype)], axis=0)
        feats[key] = merged[:MAX_TEMPLATES]

    if multimer:
        order = np.asarray(residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE)
        prepend("template_aatype", np.take(order, np.argmax(aatype, axis=-1)).astype(feats["template_aatype"].dtype))
        prepend("template_all_atom_mask", masks)
        prepend("template_all_atom_positions", positions)
        feats["num_templates"] = np.asarray(len(feats["template_aatype"]), dtype=np.int32)
    else:
        prepend("template_aatype", aatype)
        prepend("template_all_atom_masks", masks)
        prepend("template_all_atom_positions", positions)
        prepend("template_domain_names", np.asarray(names, dtype=np.object_))
        prepend("template_sequence", np.asarray([b""] * n_new, dtype=np.object_))
        prepend("template_sum_probs", np.ones((n_new, 1), np.float32))
    return n_new


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", required=True, help="features.pkl to update in place")
    ap.add_argument("--fasta", required=True, help="Query FASTA (one record per chain)")
    ap.add_argument("--templates-dir", required=True, help="Matched-templates directory (bin/fold/match_templates.py)")
    ap.add_argument("--template-chains", nargs="+", default=None, help="Chain ID(s) that may be templated (default: all)")
    ap.add_argument("--chain-layout", default=None, help="af2_mono chain layout JSON (af2_monomer_features_from_msas.py)")
    args = ap.parse_args()

    if not (Path(args.templates_dir) / "index.json").exists():
        log.info("No template index in %s; features unchanged", args.templates_dir)
        return 0
    path = Path(args.features)
    if not path.exists():
        # Shallow --msa_subsample jobs leave features.pkl to AF2's own pipeline.
        log.warning("%s not found; this job folds without user templates", path)
        return 0
    with open(path, "rb") as fh:
        feats = pickle.load(fh)
    n = add_templates(
        feats, parse_fasta(Path(args.fasta)), Path(args.templates_dir), args.template_chains,
        Path(args.chain_layout) if args.chain_layout else None,
    )
    if n:
        with open(path, "wb") as fh:
            pickle.dump(feats, fh, protocol=4)
        log.info("Added %d template row(s) to %s", n, path)
    else:
        log.info("No templates matched; features unchanged")
    return 0


if __name__ == "__main__":
    sys.exit(main())
