#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""
Bridge a ColabFold-produced a3m alignment into an AlphaFold2 "precomputed MSA"
per-target directory that modules/local/af2/alphafold2.nf (the ALPHAFOLD2 GPU
predict process) can consume via --use_precomputed_msas=true.

Ground truth (probed against the AF2 CUDA-12 custom container's
alphafold/data/pipeline.py and run_alphafold.py, 2026-07-17 - see
plans/fold-nf-multi-method-folding.md): when --use_precomputed_msas is set,
predict_structure() in run_alphafold.py loads features.pkl directly and never
re-reads msas/*.sto|*.a3m|*.hhr - so the only load-bearing artifact for the
predict stage is a valid features.pkl, not the raw per-source MSA files the
CPU MSA stage (ALPHAFOLD2_JACKHMMER_MSA) writes. This is a deliberate
deviation from the plan's original "write fake per-source .sto/.hhr files"
sketch: building features.pkl directly is simpler and more robust (it also
sidesteps the fact that AF2's monomer pipeline.process() runs template search
unconditionally - it is not gated by use_precomputed_msas at all - so faking
an empty pdb_hits.hhr would not actually have skipped it).

This script builds that features.pkl directly from the ColabFold a3m using
AlphaFold's own feature-building functions (alphafold.data.pipeline.
make_sequence_features / make_msa_features), with empty (zero-hit) template
features - ColabFold's local/remote search does not do template search in
this pipeline, and empty template arrays are valid AF2 monomer input (see
alphafold/data/templates.py: HhsearchHitFeaturizer.get_templates() returns
zero-length arrays, of the correct dtype, when given zero hits).

Caveat: a ColabFold MSA (MMseqs2 uniref30 + colabfold_envdb) is not
equivalent in depth/composition to AlphaFold's native
jackhmmer(uniref90+mgnify)+hhblits(bfd+uniref30) blend, so predictions fed
through this bridge may differ (typically modestly) from the native
jackhmmer_af2 route - see plans/fold-nf-multi-method-folding.md §3.

Phase 1 scope: monomer only (a single FASTA record / single a3m block).
"""

import argparse
import logging
import pickle
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

# The AF2 container only puts /app/alphafold on sys.path for
# /app/alphafold/run_alphafold.py itself (via Python's script-relative
# sys.path[0] insertion) - PYTHONPATH is otherwise empty, so scripts run from
# elsewhere (like this one, staged from bin/) must add it explicitly.
sys.path.insert(0, "/app/alphafold")

import numpy as np  # noqa: E402
from alphafold.data import parsers, pipeline  # noqa: E402
from alphafold.data.templates import TEMPLATE_FEATURES  # noqa: E402


def build_features(fasta_path: Path, a3m_path: Path) -> dict:
    """Build an AF2 monomer feature dict from a query FASTA + ColabFold a3m."""
    input_seqs, input_descs = parsers.parse_fasta(fasta_path.read_text())
    if len(input_seqs) != 1:
        raise ValueError(
            f"{fasta_path} has {len(input_seqs)} records; the a3m->AF2 bridge only "
            "supports monomer input in this phase (see plans/"
            "fold-nf-multi-method-folding.md Phase 1 scope)."
        )
    sequence = input_seqs[0]
    description = input_descs[0]

    msa = parsers.parse_a3m(a3m_path.read_text())
    if not msa.sequences:
        raise ValueError(f"{a3m_path} contains no sequences")
    if msa.sequences[0].replace("-", "") != sequence:
        log.warning(
            "First a3m record does not exactly match the FASTA query sequence "
            "(gaps aside) in %s - continuing, but double-check the ColabFold search "
            "was run against this same FASTA.",
            a3m_path,
        )

    sequence_features = pipeline.make_sequence_features(
        sequence=sequence, description=description, num_res=len(sequence)
    )
    msa_features = pipeline.make_msa_features((msa,))

    # See module docstring: no template search is performed for
    # ColabFold-sourced MSAs in this bridge.
    template_features = {name: np.array([], dtype=dtype) for name, dtype in TEMPLATE_FEATURES.items()}

    return {**sequence_features, **msa_features, **template_features}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="Query FASTA (single record)")
    parser.add_argument("--a3m", required=True, help="ColabFold a3m alignment (query-first)")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="AF2 per-target directory to write features.pkl into (must be named meta.id)",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    features = build_features(Path(args.fasta), Path(args.a3m))

    features_path = out_dir / "features.pkl"
    with open(features_path, "wb") as f:
        pickle.dump(features, f, protocol=4)
    log.info(
        "Wrote %s (%d MSA rows, 0 templates)",
        features_path,
        int(features["num_alignments"][0]),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
