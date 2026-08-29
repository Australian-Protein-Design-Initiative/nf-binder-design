#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///
"""
Build a *monomer* features.pkl for a multi-chain complex, using the chain-break trick.

Why this exists
---------------
AF2-multimer needs a paired MSA to reason about which homolog of chain A goes with
which homolog of chain B. For a de novo designed binder there are no homologs at all,
so `uniprot_hits` is depth 1 and no pairing happens: the merged MSA is purely block
diagonal, deep on the target and a single row on the binder. That imbalance is not an
error, but it is not what the multimer head was trained to exploit either.

The monomer models never had a pairing concept to begin with. The standard way to fold
a complex with them - used by dl_binder_design's af2_initial_guess and by ColabFold's
early complex mode - is to concatenate the chains into one sequence and jump the
residue index by a large offset at each chain break. AF2's relative-position encoding
clips at 32, so any offset above that reads as "not covalently connected". A
block-diagonal MSA (deep target block, single-sequence binder block) then carries
exactly the information we actually have, and nothing is pretending to be paired.

This gives a second AF2 opinion that does not depend on the multimer head. It is NOT an
independent engine: monomer_ptm and multimer share an architecture family and a
training corpus, so it must not be counted as another vote alongside Boltz / RF3 /
Protenix in a cross-engine consensus. Its value is as a controlled contrast - same
weights lineage, different complex-assembly assumption.

Method
------
Native monomer pipeline output, assembled by hand:
  * sequence features over the concatenated sequence, with `residue_index` carrying a
    +`--chain-break-offset` jump at every chain boundary;
  * one MSA whose row 0 is the full concatenated query, and whose remaining rows are
    each chain's own hits padded with gaps outside that chain's columns;
  * zero-hit template placeholders (no mmCIF lookup).

Must run inside the AF2 container (/app/alphafold on sys.path).
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

sys.path.insert(0, "/app/alphafold")

import numpy as np  # noqa: E402
from alphafold.common import protein  # noqa: E402
from alphafold.data import parsers, pipeline  # noqa: E402

# Reuse the loaders from the multimer builder so both stay in step.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from af2_multimer_features_from_msas import (  # noqa: E402
    chain_msa_dir,
    empty_templates,
    load_chain_msas,
    parse_fasta,
)


def block_diagonal_msa(
    sequences: Sequence[str], per_chain: Sequence[List[parsers.Msa]]
) -> parsers.Msa:
    """Pad each chain's MSA rows with gaps outside that chain's columns.

    Row 0 is the concatenated query. AF2 requires the first MSA row to be the query
    itself; every other row is one chain's hit, gapped everywhere else.
    """
    lengths = [len(s) for s in sequences]
    starts = np.cumsum([0] + lengths[:-1])
    total = sum(lengths)

    rows = ["".join(sequences)]
    deletions = [[0] * total]
    descs = ["query"]

    for idx, (start, length, msas) in enumerate(zip(starts, lengths, per_chain)):
        left, right = "-" * int(start), "-" * (total - int(start) - length)
        zeros_l, zeros_r = [0] * int(start), [0] * (total - int(start) - length)
        for msa in msas:
            for seq, dels, desc in zip(msa.sequences, msa.deletion_matrix, msa.descriptions):
                if len(seq) != length:
                    raise ValueError(
                        f"chain {idx} MSA row '{desc}' has length {len(seq)}, "
                        f"expected {length} - the MSA does not match the query"
                    )
                rows.append(left + seq + right)
                deletions.append(zeros_l + list(dels) + zeros_r)
                descs.append(desc)

    return parsers.Msa(sequences=rows, deletion_matrix=deletions, descriptions=descs)


def build_features(
    fasta_path: Path, msas_root: Path, chain_break_offset: int
) -> Tuple[pipeline.FeatureDict, List[Tuple[str, int]]]:
    records = parse_fasta(fasta_path)
    if not records:
        raise ValueError(f"{fasta_path} has no FASTA records")

    sequences = [seq for _d, seq in records]
    descriptions = [d for d, _s in records]
    concat = "".join(sequences)

    per_chain: List[List[parsers.Msa]] = []
    for chain_id, seq, desc in zip(protein.PDB_CHAIN_IDS, sequences, descriptions):
        chain_dir = chain_msa_dir(msas_root, chain_id)
        msas = load_chain_msas(chain_dir)
        depth = sum(len(m.sequences) for m in msas)
        log.info("chain %s (%s): %d residues, %d MSA rows from %s",
                 chain_id, desc.split()[0] if desc else chain_id, len(seq), depth, chain_dir)
        per_chain.append(msas)

    feats = pipeline.make_sequence_features(
        sequence=concat, description=fasta_path.stem, num_res=len(concat)
    )

    # Chain break: +offset at every boundary. AF2 clips relative positions at 32, so
    # anything above that reads as a discontinuity; 200 is the dl_binder_design value.
    residue_index = []
    cursor = 0
    for i, seq in enumerate(sequences):
        if i:
            cursor += chain_break_offset
        residue_index.extend(range(cursor, cursor + len(seq)))
        cursor += len(seq)
    feats["residue_index"] = np.array(residue_index, dtype=np.int32)

    merged = block_diagonal_msa(sequences, per_chain)
    feats.update(pipeline.make_msa_features([merged]))
    feats.update(empty_templates(len(concat)))

    layout = list(zip(list(protein.PDB_CHAIN_IDS)[: len(sequences)],
                      [len(s) for s in sequences]))
    return feats, layout


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="Query FASTA (one record per chain)")
    ap.add_argument("--msas-dir", required=True,
                    help="AF2 per-target directory (contains msas/A, msas/B, ...)")
    ap.add_argument("--chain-break-offset", type=int, default=200,
                    help="residue_index jump at each chain break (default 200)")
    ap.add_argument("--layout-out", help="write chain layout JSON for the splitter")
    args = ap.parse_args()

    msas_root = Path(args.msas_dir)
    feats, layout = build_features(Path(args.fasta), msas_root, args.chain_break_offset)

    out = msas_root / "features.pkl"
    with open(out, "wb") as fh:
        pickle.dump(feats, fh, protocol=4)
    log.info("Wrote %s (msa %s, aatype %s, residue_index %d..%d)",
             out, feats["msa"].shape, feats["aatype"].shape,
             int(feats["residue_index"][0]), int(feats["residue_index"][-1]))

    if args.layout_out:
        Path(args.layout_out).write_text(json.dumps({
            "chains": [{"chain_id": c, "length": n} for c, n in layout],
            "chain_break_offset": args.chain_break_offset,
        }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
