#!/usr/bin/env python
# /// script
# requires-python = ">=3.10"
# ///

"""
Chai-1 structure prediction driver for fold.nf's CHAI1_FOLD.

The `chai-lab fold` CLI keeps the PAE matrix in memory only, so this script calls
chai_lab.chai1.run_inference() itself and writes the PAE alongside each structure
for ipSAE. Inputs are prepared from the pipeline's shared MSA stage and matched
templates; Chai-1 is never pointed at its MSA or template servers.

Conventions this script has to honour (chai-lab 0.6.1):
  - MSAs are read from <msa_directory>/<sha256(SEQUENCE)>.aligned.pqt, a parquet
    with columns sequence, source_database, pairing_key, comment. Row 0 is the
    query (source_database "query").
  - MULTIMER PAIRING IS DONE BY THE LIBRARY: rows of different chains that share a
    pairing_key are paired (data/dataset/msas/preprocess.py), and rows with an
    empty key stay as that chain's unpaired tail. The key is an arbitrary string,
    so we use the numeric taxid parsed by msa_taxonomy.py (Chai's own jackhmmer
    route keys UniProt hits on the OS= species name instead).
  - source_database is a model feature. Hits get "uniref90" when the header
    starts with UniRef and "bfd_uniclust" otherwise, as in Chai's own ColabFold
    route (data/dataset/msas/colabfold.py). Chai caps uniref90 rows at 10,000.
  - Templates come from an m8 hit table whose query_id is the FASTA entity name
    and whose subject_id is <ID>_<chain>. Chai loads <ID upper>.cif.gz from
    CHAI_TEMPLATE_CIF_FOLDER and only downloads from RCSB when it is missing, so
    the matched templates are placed there first. Chai realigns each hit with
    kalign and uses at most 4 templates per chain.

Outputs per diffusion sample, in the same shape as run_esmfold2.py so the shared
FOLD_PARSE_CONFIDENCE / ipsae.py af3 path consumes them unchanged:

    <name>_seed_<S>_sample_<i>_model.cif
    <name>_seed_<S>_sample_<i>_confidences.json          (atom_plddts, pae)
    <name>_seed_<S>_sample_<i>_summary_confidences.json  (ranking_score, ptm, iptm, chain_pair_iptm)
"""

import argparse
import gzip
import hashlib
import json
import logging
import os
import random
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional

from make_af3_input import CHAIN_IDS, _match_per_chain, clean_a3m, parse_fasta_records
from msa_taxonomy import parse_a3m
from run_esmfold2 import atom_plddts_from_cif

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

MAX_TEMPLATES = 4


def pqt_basename(seq: str) -> str:
    """Chai-1's expected_basename(): sha256 of the upper-cased query sequence."""
    return hashlib.sha256(seq.upper().encode()).hexdigest() + ".aligned.pqt"


def msa_rows(a3m_text: str) -> List[Dict[str, str]]:
    """Rows of a Chai-1 aligned.pqt, query first, pairing_key = taxid (or blank)."""
    rows = []
    for rec in parse_a3m(a3m_text):
        if rec.is_query:
            src, key = "query", ""
        else:
            src = "uniref90" if rec.header.startswith("UniRef") else "bfd_uniclust"
            key = rec.tax_id or ""
        rows.append({"sequence": rec.sequence, "source_database": src, "pairing_key": key, "comment": rec.header})
    return rows


def template_rows(templates_dir: Path, seq: str, chain_id: str, cif_folder: Path) -> List[list]:
    """m8 rows for one chain's matched templates; copies each CIF in as <ID>.cif.gz."""
    index_path = templates_dir / "index.json"
    if not index_path.exists():
        return []
    hits = json.loads(index_path.read_text()).get(hashlib.md5(seq.upper().encode()).hexdigest(), [])
    rows = []
    for rank, h in enumerate(hits[:MAX_TEMPLATES]):
        stem = Path(h["cif"]).stem
        if "_" in stem:
            raise ValueError(f"template file {h['cif']}: Chai-1 splits template ids on '_'")
        dest = cif_folder / f"{stem.upper()}.cif.gz"
        if not dest.exists():
            with open(templates_dir / h["cif"], "rb") as src, gzip.open(dest, "wb") as out:
                shutil.copyfileobj(src, out)
        q, t = h["query_indices"], h["template_indices"]
        rows.append([
            chain_id, f"{stem}_{h['template_chain']}", round(100 * h.get("identity", 0), 1), len(q), 0, 0,
            min(q) + 1, max(q) + 1, min(t) + 1, max(t) + 1, rank, 0, "",
        ])
    return rows


def prepare_inputs(
    fasta_path: Path,
    work_dir: Path,
    a3m_paths: Optional[List[Path]] = None,
    templates_dir: Optional[Path] = None,
    template_chains: Optional[List[str]] = None,
) -> Dict[str, object]:
    """Write Chai-1's FASTA, per-sequence MSA rows and template table under work_dir."""
    sequences = parse_fasta_records(fasta_path)
    if not sequences:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    if len(sequences) > len(CHAIN_IDS):
        raise ValueError(f"{len(sequences)} chains exceeds the {len(CHAIN_IDS)} supported")
    per_chain = _match_per_chain(a3m_paths, len(sequences), "--a3m")

    work_dir.mkdir(parents=True, exist_ok=True)
    cif_folder = work_dir / "template_cifs"
    cif_folder.mkdir(exist_ok=True)
    fasta_out = work_dir / "chai_input.fasta"
    fasta_out.write_text("".join(f">protein|name={CHAIN_IDS[i]}\n{s}\n" for i, s in enumerate(sequences)))

    msas: Dict[str, List[Dict[str, str]]] = {}
    m8: List[list] = []
    for i, seq in enumerate(sequences):
        cid = CHAIN_IDS[i]
        if per_chain[i] is not None and pqt_basename(seq) not in msas:
            rows = msa_rows(clean_a3m(Path(per_chain[i]), seq, cid))
            msas[pqt_basename(seq)] = rows
            log.info("chain %s: %d MSA rows, %d keyed by taxid", cid, len(rows), sum(1 for r in rows if r["pairing_key"]))
        if templates_dir is not None and (not template_chains or cid in template_chains):
            rows = template_rows(templates_dir, seq, cid, cif_folder)
            if rows:
                log.info("chain %s: %d template(s)", cid, len(rows))
            m8 += rows

    m8_path = None
    if m8:
        m8_path = work_dir / "templates.m8"
        m8_path.write_text("".join("\t".join(str(v) for v in r) + "\n" for r in m8))
    return {"fasta": fasta_out, "msas": msas, "m8": m8_path, "cif_folder": cif_folder, "n_chains": len(sequences)}


def summary_json(ranking, plddt_mean: float, n_chains: int, seed: int, sample: int) -> Dict[str, object]:
    """AF3-shaped summary from a chai SampleRanking (all tensors have a leading batch dim of 1)."""
    pair = ranking.ptm_scores.per_chain_pair_iptm[0].tolist()
    return {
        "ranking_score": float(ranking.aggregate_score.item()),
        "ptm": float(ranking.ptm_scores.complex_ptm.item()),
        "iptm": float(ranking.ptm_scores.interface_ptm.item()),
        "plddt": plddt_mean,
        "has_clash": bool(ranking.clash_scores.has_inter_chain_clashes.any().item()),
        "chain_pair_iptm": [[float(v) for v in row] for row in pair],
        "chain_ids": list(CHAIN_IDS[:n_chains]),
        "seed": seed,
        "sample": sample,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True, type=Path, help="Input FASTA; one record per chain.")
    p.add_argument("--name", required=True, help="Output file prefix (fold.nf meta.id).")
    p.add_argument("--output-dir", type=Path, default=Path("output"))
    p.add_argument("--a3m", type=Path, nargs="*", help="One a3m per FASTA record, in record order (none: single sequence).")
    p.add_argument("--templates-dir", type=Path, help="Matched-templates directory (bin/fold/match_templates.py)")
    p.add_argument("--template-chains", nargs="+", help="Chain ID(s) that may be templated (default: all)")
    p.add_argument("--seed", type=int, help="Model seed (default: random, recorded in the file names)")
    p.add_argument("--num-diffusion-samples", type=int, default=5)
    p.add_argument("--num-trunk-recycles", type=int, default=3)
    p.add_argument("--num-diffn-timesteps", type=int, default=200)
    p.add_argument("--no-esm-embeddings", action="store_true", help="Disable the ESM2 language-model embeddings")
    p.add_argument("--work-dir", type=Path, default=Path("chai_inputs"))
    p.add_argument("--device", default=None)
    args = p.parse_args()

    seed = args.seed if args.seed is not None else random.SystemRandom().randrange(2**31)
    inputs = prepare_inputs(args.fasta, args.work_dir, args.a3m, args.templates_dir, args.template_chains)
    os.environ["CHAI_TEMPLATE_CIF_FOLDER"] = str(inputs["cif_folder"].resolve())

    import pandas as pd
    import torch
    from chai_lab.chai1 import run_inference

    msa_dir = args.work_dir / "msas"
    msa_dir.mkdir(exist_ok=True)
    for basename, rows in inputs["msas"].items():
        pd.DataFrame.from_records(rows).to_parquet(msa_dir / basename)

    log.info(
        "Chai-1 %s: %d chain(s), %d MSA(s), templates %s, %d sample(s), seed %d",
        args.name, inputs["n_chains"], len(inputs["msas"]), "yes" if inputs["m8"] else "no",
        args.num_diffusion_samples, seed,
    )
    raw_dir = args.work_dir / "raw"
    candidates = run_inference(
        fasta_file=inputs["fasta"],
        output_dir=raw_dir,
        use_esm_embeddings=not args.no_esm_embeddings,
        msa_directory=msa_dir if inputs["msas"] else None,
        template_hits_path=inputs["m8"],
        num_trunk_recycles=args.num_trunk_recycles,
        num_diffn_timesteps=args.num_diffn_timesteps,
        num_diffn_samples=args.num_diffusion_samples,
        seed=seed,
        device=args.device,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for i, (cif, ranking) in enumerate(zip(candidates.cif_paths, candidates.ranking_data, strict=True)):
        prefix = args.output_dir / f"{args.name}_seed_{seed}_sample_{i}"
        cif_text = Path(cif).read_text()
        Path(f"{prefix}_model.cif").write_text(cif_text)
        token_plddts = candidates.plddt[i].tolist()
        confidences = {
            "atom_plddts": atom_plddts_from_cif(cif_text),
            "pae": candidates.pae[i].tolist(),
            "token_plddts": token_plddts,
        }
        Path(f"{prefix}_confidences.json").write_text(json.dumps(confidences))
        summary = summary_json(ranking, float(torch.tensor(token_plddts).mean()), inputs["n_chains"], seed, i)
        Path(f"{prefix}_summary_confidences.json").write_text(json.dumps(summary, indent=2))
        log.info("%s: score %.3f pTM %.3f ipTM %.3f", prefix.name, summary["ranking_score"], summary["ptm"], summary["iptm"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
