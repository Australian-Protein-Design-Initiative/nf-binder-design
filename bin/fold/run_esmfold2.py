#!/usr/bin/env python
# /// script
# requires-python = ">=3.12"
# ///

"""
ESMFold2 structure prediction driver for fold.nf's ESMFOLD2_FOLD.

The `esm` package (3.4.1.post1) ships NO command-line entry point - no
`[project.scripts]`, no `__main__.py` - so unlike every other fold.nf engine this
one is driven by our own script around the documented Python API:

    model = EsmFold2Model.from_pretrained("biohub/ESMFold2", device="cuda").eval()
    result = ESMFold2InputBuilder().fold(model, spi, num_diffusion_samples=N, seed=S)

Conventions this script has to honour (esm 3.4.1.post1):
  - MULTIMER PAIRING IS DONE BY THE LIBRARY. We pass one a3m per chain and it
    builds the paired + block-diagonal layout itself, reading `key=<taxid>` from
    each hit header (models/esmfold2/paired_msa.py). Render those headers with
    `bin/fold/msa_taxonomy.py --tool esmfold2`; rows with no `key=` stay in the
    file as that chain's unpaired tail.
  - `fold()` returns a bare MolecularComplexResult for num_diffusion_samples=1
    and a LIST otherwise, so the result is always normalised to a list here.
  - result.plddt is per token on a 0-1 scale (not 0-100), result.pae is a token
    x token matrix in Angstroms, and result.pde / .num_tokens are never
    populated by decode().
  - The `-Fast` checkpoint has no MSA encoder and silently IGNORES MSAs, so MSA
    runs need the full `biohub/ESMFold2` (see --model).

Outputs per diffusion sample, named so the shared FOLD_PARSE_CONFIDENCE /
ipsae.py path can consume them unchanged:

    <name>_seed_<S>_sample_<i>_model.cif
    <name>_seed_<S>_sample_<i>_confidences.json          (atom_plddts, pae)
    <name>_seed_<S>_sample_<i>_summary_confidences.json  (ptm, iptm, chain_pair_iptm)

Those two JSONs deliberately use AlphaFold3's key names and shapes: ipsae.py's
af3 branch reads `atom_plddts` indexed by mmCIF atom serial plus a sibling
*_summary_confidences.json with a chain-ordered `chain_pair_iptm` matrix, so
emitting that shape needs no new ipsae format. `atom_plddts` is read back out of
the mmCIF B-factor column we just wrote (to_mmcif() puts pLDDT there), which
guarantees it is indexed exactly like the atoms ipsae.py parses from that same
file rather than relying on an assumed token-to-atom expansion.

ESMFold2 reports no ranking score, so that column is left blank downstream -
rank on iptm / plddt / ipsae instead.
"""

import argparse
import json
import logging
import string
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from make_af3_input import _match_per_chain, parse_fasta_records, read_a3m_records

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

CHAIN_IDS = string.ascii_uppercase


@dataclass
class ChainSpec:
    """One protein chain: its ESMFold2 chain id, sequence and optional a3m."""

    chain_id: str
    sequence: str
    a3m: Optional[Path] = None


def build_chain_specs(
    fasta_path: Path,
    a3m_paths: Optional[List[Path]] = None,
    single_sequence: bool = False,
) -> List[ChainSpec]:
    """Pair FASTA records with per-chain a3ms, in record order.

    The a3m's first row must be the chain's own query sequence - ESMFold2 treats
    row 0 as the query, so a mismatch silently folds the wrong alignment.
    """
    sequences = parse_fasta_records(fasta_path)
    if not sequences:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    if len(sequences) > len(CHAIN_IDS):
        raise ValueError(f"{len(sequences)} chains exceeds the {len(CHAIN_IDS)} supported")

    per_chain = [None] * len(sequences) if single_sequence \
        else _match_per_chain(a3m_paths, len(sequences), "--a3m")

    specs: List[ChainSpec] = []
    for i, seq in enumerate(sequences):
        cid = CHAIN_IDS[i]
        a3m = per_chain[i]
        if a3m is not None:
            records = read_a3m_records(Path(a3m).read_text())
            if not records:
                log.warning("chain %s: %s has no records - folding it single-sequence", cid, Path(a3m).name)
                a3m = None
            else:
                first = records[0][1].replace("-", "")
                if first.upper() != seq.upper():
                    raise ValueError(
                        f"chain {cid}: first sequence in {Path(a3m).name} does not match the "
                        f"FASTA sequence (the query must be the first MSA row)"
                    )
        specs.append(ChainSpec(chain_id=cid, sequence=seq, a3m=Path(a3m) if a3m is not None else None))
    return specs


def atom_plddts_from_cif(cif_text: str) -> List[float]:
    """Per-atom B-factors from an mmCIF, indexed by (atom serial - lowest serial).

    to_mmcif() writes pLDDT into the B-factor column, and ipsae.py indexes its
    atom_plddts array by `_atom_site.id` minus the lowest serial in the file, so
    reading the same file back reproduces exactly that indexing.
    """
    tags: List[str] = []
    rows: List[List[str]] = []
    in_loop = False
    for line in cif_text.splitlines():
        s = line.strip()
        if s.startswith("_atom_site."):
            tags.append(s.split(".", 1)[1].split()[0])
            in_loop = True
            continue
        if in_loop:
            if not s or s.startswith("#") or s.startswith("_") or s == "loop_":
                break
            rows.append(s.split())

    if not tags or not rows:
        raise ValueError("mmCIF has no _atom_site loop to read B-factors from")
    for required in ("id", "B_iso_or_equiv"):
        if required not in tags:
            raise ValueError(f"mmCIF _atom_site loop has no {required} column")
    i_id = tags.index("id")
    i_b = tags.index("B_iso_or_equiv")

    serials: List[int] = []
    values: List[float] = []
    for row in rows:
        if len(row) <= max(i_id, i_b):
            continue
        try:
            serials.append(int(row[i_id]))
            values.append(float(row[i_b]))
        except ValueError:
            continue
    if not serials:
        raise ValueError("mmCIF _atom_site loop has no parseable atom rows")

    base = min(serials)
    plddts = [0.0] * (max(serials) - base + 1)
    for serial, value in zip(serials, values):
        plddts[serial - base] = value
    return plddts


def chain_pair_iptm_matrix(pair_chains_iptm, n_chains: int) -> Optional[List[List[float]]]:
    """ESMFold2's [n_chains, n_chains] ipTM as a plain nested list, or None.

    ipsae.py's af3 branch indexes this by chain letter position (A->0, B->1),
    which is the order build_chain_specs assigns, so no remapping is needed.
    """
    if pair_chains_iptm is None:
        return None
    matrix = [[float(v) for v in row] for row in pair_chains_iptm.tolist()]
    if len(matrix) != n_chains:
        log.warning(
            "pair_chains_iptm is %dx%d for %d chains - writing it through unchanged",
            len(matrix), len(matrix[0]) if matrix else 0, n_chains,
        )
    return matrix


def sample_outputs(result, specs: Sequence[ChainSpec], seed: int, sample: int) -> Dict[str, object]:
    """Build the mmCIF text plus the two confidence JSONs for one sample."""
    cif_text = result.complex.to_mmcif()
    plddt = result.plddt
    mean_plddt = float(plddt.mean()) if plddt is not None else None

    confidences: Dict[str, object] = {"atom_plddts": atom_plddts_from_cif(cif_text)}
    if result.pae is not None:
        confidences["pae"] = result.pae.tolist()
    if plddt is not None:
        # Named token_plddts, not plddt: ipsae.py's normalize_token_pae_json()
        # aliases a bare "plddt" key onto atom_plddts, which this is not.
        confidences["token_plddts"] = [float(v) for v in plddt.tolist()]

    summary: Dict[str, object] = {
        "ptm": None if result.ptm is None else float(result.ptm),
        "iptm": None if result.iptm is None else float(result.iptm),
        "plddt": mean_plddt,
        "chain_ids": [s.chain_id for s in specs],
        "seed": seed,
        "sample": sample,
    }
    matrix = chain_pair_iptm_matrix(result.pair_chains_iptm, len(specs))
    if matrix is not None:
        summary["chain_pair_iptm"] = matrix
    return {"cif": cif_text, "confidences": confidences, "summary": summary}


def write_sample(out_dir: Path, name: str, seed: int, sample: int, outputs: Dict[str, object]) -> str:
    prefix = f"{name}_seed_{seed}_sample_{sample}"
    (out_dir / f"{prefix}_model.cif").write_text(outputs["cif"])
    (out_dir / f"{prefix}_confidences.json").write_text(json.dumps(outputs["confidences"]))
    (out_dir / f"{prefix}_summary_confidences.json").write_text(json.dumps(outputs["summary"], indent=2))
    return prefix


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--fasta", required=True, type=Path, help="Input FASTA; one record per chain.")
    p.add_argument("--name", required=True, help="Output file prefix (fold.nf meta.id).")
    p.add_argument("--output-dir", type=Path, default=Path("output"), help="Directory for the outputs.")
    p.add_argument(
        "--a3m", type=Path, nargs="*",
        help="One a3m per FASTA record, in record order. Multimer a3ms should carry "
             "key=<taxid> headers (msa_taxonomy.py --tool esmfold2).",
    )
    p.add_argument(
        "--single-sequence", action="store_true",
        help="Ignore --a3m and fold from sequence alone (ESMFold2's single-sequence mode).",
    )
    p.add_argument("--model", default="biohub/ESMFold2", help="HuggingFace repo id or local snapshot dir.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--num-diffusion-samples", type=int, default=1)
    p.add_argument("--num-loops", type=int, help="esm default: 20")
    p.add_argument("--num-sampling-steps", type=int, help="esm default: 200")
    p.add_argument("--msa-max-depth", type=int, help="Rows kept per loop; esm default: 1024")
    p.add_argument(
        "--kernel-backend", default="fused", choices=["fused", "cuequivariance", "none"],
        help="model.set_kernel_backend(); 'fused' falls back to the reference path without Triton.",
    )
    p.add_argument("--esmc-precision", default="bf16", choices=["bf16", "fp32", "fp8"])
    p.add_argument("--chunk-size", type=int, help="model.set_chunk_size(); esm default 64, 0 disables chunking.")
    p.add_argument("--ccd-cache", type=Path, help="Directory holding ccd.pkl (else it is fetched from HuggingFace).")
    p.add_argument("--local-files-only", action="store_true", help="Never reach HuggingFace for weights.")
    p.add_argument("--device", default="cuda")
    args = p.parse_args()

    specs = build_chain_specs(args.fasta, args.a3m, args.single_sequence)
    n_with_msa = sum(1 for s in specs if s.a3m is not None)
    log.info(
        "ESMFold2 %s: %d chain(s), %d with an MSA, %d diffusion sample(s), seed %d",
        args.name, len(specs), n_with_msa, args.num_diffusion_samples, args.seed,
    )
    if n_with_msa and "fast" in str(args.model).lower():
        log.warning(
            "--model %s looks like the -Fast checkpoint, which has no MSA encoder and "
            "IGNORES MSAs silently. Use biohub/ESMFold2 for MSA runs.", args.model
        )

    import torch
    from esm.models.esmfold2 import (
        ESMFold2InputBuilder,
        EsmFold2Model,
        ProteinInput,
        StructurePredictionInput,
    )
    from esm.utils.msa import MSA

    model = EsmFold2Model.from_pretrained(
        args.model,
        device=args.device,
        esmc_precision=args.esmc_precision,
        local_files_only=args.local_files_only,
    ).eval()
    model.set_kernel_backend(None if args.kernel_backend == "none" else args.kernel_backend)
    if args.chunk_size is not None:
        model.set_chunk_size(None if args.chunk_size == 0 else args.chunk_size)

    builder = ESMFold2InputBuilder(ccd_cache=args.ccd_cache)
    spi = StructurePredictionInput(
        sequences=[
            ProteinInput(
                id=s.chain_id,
                sequence=s.sequence,
                msa=MSA.from_a3m(str(s.a3m)) if s.a3m is not None else None,
            )
            for s in specs
        ]
    )

    fold_kwargs = {"num_diffusion_samples": args.num_diffusion_samples, "seed": args.seed}
    for key, value in (
        ("num_loops", args.num_loops),
        ("num_sampling_steps", args.num_sampling_steps),
        ("msa_max_depth", args.msa_max_depth),
    ):
        if value is not None:
            fold_kwargs[key] = value

    with torch.inference_mode():
        results = builder.fold(model, spi, **fold_kwargs)
    if not isinstance(results, list):
        results = [results]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for i, result in enumerate(results):
        prefix = write_sample(args.output_dir, args.name, args.seed, i, sample_outputs(result, specs, args.seed, i))
        log.info(
            "%s: pLDDT %.3f pTM %s ipTM %s",
            prefix,
            float(result.plddt.mean()) if result.plddt is not None else float("nan"),
            f"{result.ptm:.3f}" if result.ptm is not None else "n/a",
            f"{result.iptm:.3f}" if result.iptm is not None else "n/a",
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
