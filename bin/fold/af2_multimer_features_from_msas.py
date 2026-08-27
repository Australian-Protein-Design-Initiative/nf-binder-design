#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///
"""
Build a multimer (or monomer) features.pkl from an AF2 precomputed-MSA directory.

The custom AF2 container's predict_structure() loads features.pkl directly when
--use_precomputed_msas=true and never re-reads msas/*.sto|*.a3m|*.hhr. fold_pulldown
assembles those per-chain files but cannot run AF2's data pipeline (template search
is not gated by use_precomputed_msas). This script builds features.pkl the same way
native pipeline_multimer.DataPipeline.process does after MSA search: per-chain
make_sequence_features / make_msa_features, convert_monomer_features,
add_assembly_features, pair_and_merge, pad_msa.

Empty templates are used (no mmCIF lookup). Ranked zero-hit arrays are required:
np.array([]) makes convert_monomer_features argmax crash.

Must run inside the AF2 container (/app/alphafold on sys.path).
"""

from __future__ import annotations

import argparse
import logging
import pickle
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

sys.path.insert(0, "/app/alphafold")

import numpy as np  # noqa: E402
from alphafold.common import protein  # noqa: E402
from alphafold.data import feature_processing, msa_pairing, parsers, pipeline  # noqa: E402
from alphafold.data.pipeline_multimer import (  # noqa: E402
    add_assembly_features,
    convert_monomer_features,
    pad_msa,
)

STO_NAMES = ("uniref90_hits.sto", "mgnify_hits.sto")
A3M_CANDIDATES = ("bfd_uniref_hits.a3m", "bfd_uniclust_hits.a3m")
UNIPROT_STO = "uniprot_hits.sto"
N_ATOM_TYPES = 37
N_TEMPLATE_AATYPE = 22


def empty_templates(n_res: int) -> Dict[str, np.ndarray]:
    """Zero-hit templates with rank so convert_monomer_features can argmax."""
    return {
        "template_aatype": np.zeros((0, n_res, N_TEMPLATE_AATYPE), dtype=np.float32),
        "template_all_atom_masks": np.zeros((0, n_res, N_ATOM_TYPES), dtype=np.float32),
        "template_all_atom_positions": np.zeros(
            (0, n_res, N_ATOM_TYPES, 3), dtype=np.float32
        ),
        "template_domain_names": np.array([], dtype=object),
        "template_sequence": np.array([], dtype=object),
        "template_sum_probs": np.array([], dtype=np.float32),
    }


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    seqs, descs = parsers.parse_fasta(path.read_text())
    return list(zip(descs, seqs))


def chain_msa_dir(msas_root: Path, chain_id: str) -> Path:
    for candidate in (msas_root / "msas" / chain_id, msas_root / chain_id):
        if candidate.is_dir() and any(candidate.iterdir()):
            return candidate
    if chain_id == "A":
        flat = msas_root / "msas"
        if flat.is_dir() and (
            any(flat.glob("*.sto")) or any(flat.glob("*.a3m"))
        ):
            files_only = [p for p in flat.iterdir() if p.is_file()]
            subdirs = [p for p in flat.iterdir() if p.is_dir()]
            if files_only and not subdirs:
                return flat
    raise FileNotFoundError(f"no MSA directory for chain {chain_id} under {msas_root}")


def _parse_stockholm(path: Path) -> Optional[parsers.Msa]:
    try:
        msa = parsers.parse_stockholm(path.read_text())
    except (IndexError, ValueError) as exc:
        # Duplicate hit IDs in a3m-converted STO files make AF2 concatenate
        # rows by name and then index past the end of shorter sequences.
        log.warning("skipping unparseable Stockholm %s: %s", path, exc)
        return None
    if not msa.sequences:
        log.warning("empty Stockholm %s", path)
        return None
    return msa


def load_chain_msas(chain_dir: Path) -> List[parsers.Msa]:
    msas: List[parsers.Msa] = []
    for name in STO_NAMES:
        path = chain_dir / name
        if path.is_file():
            msa = _parse_stockholm(path)
            if msa is not None:
                msas.append(msa)
    a3m_path: Optional[Path] = None
    for name in A3M_CANDIDATES:
        candidate = chain_dir / name
        if candidate.is_file():
            a3m_path = candidate
            break
    if a3m_path is not None:
        msa = parsers.parse_a3m(a3m_path.read_text())
        if msa.sequences:
            msas.append(msa)
        else:
            log.warning("empty a3m %s", a3m_path)
    if not msas:
        raise ValueError(f"no usable MSA files in {chain_dir}")
    return msas


def load_uniprot_msa(chain_dir: Path, fallback: parsers.Msa) -> parsers.Msa:
    path = chain_dir / UNIPROT_STO
    if path.is_file():
        msa = _parse_stockholm(path)
        if msa is not None:
            return msa
        log.warning("unusable uniprot Stockholm %s; using main MSA for pairing", path)
    return fallback


def all_seq_features(msa: parsers.Msa) -> Dict[str, np.ndarray]:
    feats = pipeline.make_msa_features([msa])
    valid = msa_pairing.MSA_FEATURES + ("msa_species_identifiers",)
    return {f"{k}_all_seq": v for k, v in feats.items() if k in valid}


def monomer_features(
    sequence: str,
    description: str,
    msas: Sequence[parsers.Msa],
    pairing_msa: parsers.Msa,
    *,
    add_pairing: bool,
) -> pipeline.FeatureDict:
    seq_feats = pipeline.make_sequence_features(
        sequence=sequence, description=description, num_res=len(sequence)
    )
    msa_feats = pipeline.make_msa_features(tuple(msas))
    feats: pipeline.FeatureDict = {
        **seq_feats,
        **msa_feats,
        **empty_templates(len(sequence)),
    }
    if add_pairing:
        feats.update(all_seq_features(pairing_msa))
    return feats


def build_features(fasta_path: Path, msas_root: Path) -> pipeline.FeatureDict:
    records = parse_fasta(fasta_path)
    if not records:
        raise ValueError(f"{fasta_path} has no FASTA records")
    if len(records) > protein.PDB_MAX_CHAINS:
        raise ValueError(
            f"{fasta_path} has {len(records)} chains; PDB format supports "
            f"at most {protein.PDB_MAX_CHAINS}"
        )

    unique_seqs = {seq for _desc, seq in records}
    is_homomer_or_monomer = len(unique_seqs) == 1
    add_pairing = not is_homomer_or_monomer

    all_chain_features: Dict[str, pipeline.FeatureDict] = {}
    by_sequence: Dict[str, pipeline.FeatureDict] = {}
    for chain_id, (description, sequence) in zip(protein.PDB_CHAIN_IDS, records):
        if sequence in by_sequence:
            all_chain_features[chain_id] = {
                k: (v.copy() if hasattr(v, "copy") else v)
                for k, v in by_sequence[sequence].items()
            }
            continue
        chain_dir = chain_msa_dir(msas_root, chain_id)
        msas = load_chain_msas(chain_dir)
        pairing = load_uniprot_msa(chain_dir, msas[0])
        chain_feats = monomer_features(
            sequence, description, msas, pairing, add_pairing=add_pairing
        )
        chain_feats = convert_monomer_features(chain_feats, chain_id=chain_id)
        all_chain_features[chain_id] = chain_feats
        by_sequence[sequence] = chain_feats
        log.info(
            "chain %s (%s): %d MSA rows from %s",
            chain_id,
            description.split()[0] if description else chain_id,
            int(chain_feats["num_alignments"]),
            chain_dir,
        )

    all_chain_features = add_assembly_features(all_chain_features)
    merged = feature_processing.pair_and_merge(all_chain_features=all_chain_features)
    return pad_msa(merged, 512)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--fasta", required=True, help="Query FASTA (one record per chain)")
    parser.add_argument(
        "--msas-dir",
        required=True,
        help="AF2 per-target directory (contains msas/A, msas/B, ... or flat msas/)",
    )
    args = parser.parse_args()

    msas_root = Path(args.msas_dir)
    features = build_features(Path(args.fasta), msas_root)
    features_path = msas_root / "features.pkl"
    with open(features_path, "wb") as f:
        pickle.dump(features, f, protocol=4)
    log.info(
        "Wrote %s (msa %s, aatype %s)",
        features_path,
        features["msa"].shape,
        features["aatype"].shape,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
