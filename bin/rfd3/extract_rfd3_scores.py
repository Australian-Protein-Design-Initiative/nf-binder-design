#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pandas",
# ]
# ///

"""
Extract metrics from RFDiffusion3 and RosettaFold3 output JSONs to TSV.
"""

import argparse
import json
import logging
import os
import re
import sys
from pathlib import Path
from typing import Any, Optional

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

RFD3_METRICS_SKIP = {"num_residues", "num_residues_in", "diffused_com", "fixed_com"}
RFD3_METRICS_RENAME = {
    "n_clashing.interresidue_clashes_w_sidechain": "n_clashing_w_sidechain",
    "n_clashing.interresidue_clashes_w_backbone": "n_clashing_w_backbone",
}
# RFD3 CIF names use ".cif_b0_d0"; RF3 summary stems often use "_cif_b0_d0" after an infix (e.g. _rf3_).
BACKBONE_SUFFIX_RE = re.compile(r"(?:\.|_)cif_b\d+_d\d+$")


def _first_off_diagonal_non_null(matrix: list[list[Any]]) -> Optional[float]:
    """Return first non-null off-diagonal value from 2x2 (or larger) matrix."""
    for i, row in enumerate(matrix):
        for j, v in enumerate(row):
            if i != j and v is not None:
                return float(v)
    return None


def extract_rfdiffusion3(json_path: Path) -> dict[str, Any]:
    with open(json_path) as f:
        data = json.load(f)

    stem = json_path.stem
    spec = data.get("specification") or {}
    metrics = data.get("metrics") or {}

    row: dict[str, Any] = {"backbone_id": stem}
    row["rfd3_input"] = os.path.basename(spec.get("input", ""))
    row["rfd3_contig"] = spec.get("contig", "")

    hotspots = spec.get("select_hotspots") or {}
    row["rfd3_hotspots"] = " ".join(f"{k}:{v}" for k, v in sorted(hotspots.items()))

    for key, val in metrics.items():
        if key in RFD3_METRICS_SKIP or isinstance(val, (list, dict)):
            continue
        out_key = RFD3_METRICS_RENAME.get(key, key)
        row[f"rfd3_{out_key}"] = val

    return row


def _chain_ptm_indices_for_target_binder(
    data: dict[str, Any],
    chain_ptm: list,
    target_chain: Optional[str],
    binder_chain: Optional[str],
) -> tuple[Optional[int], Optional[int]]:
    """Map target/binder chain IDs to chain_ptm indices (two-chain models)."""
    if not target_chain or not binder_chain or len(chain_ptm) < 2:
        return None, None
    tc = target_chain.strip()
    bc = binder_chain.strip()
    cids = data.get("chain_ids")
    if isinstance(cids, list) and len(cids) == len(chain_ptm):

        def find_idx(ch: str) -> Optional[int]:
            for i, x in enumerate(cids):
                if str(x).strip() == ch:
                    return i
            return None

        it, ib = find_idx(tc), find_idx(bc)
        if it is not None and ib is not None:
            return it, ib
    order = sorted([tc, bc])
    try:
        return order.index(tc), order.index(bc)
    except ValueError:
        return None, None


def extract_rf3(
    json_path: Path,
    target_chain: Optional[str] = None,
    binder_chain: Optional[str] = None,
) -> dict[str, Any]:
    with open(json_path) as f:
        data = json.load(f)

    name = json_path.name
    if name.endswith("_summary_confidences.json"):
        stem = name[: -len("_summary_confidences.json")]
    else:
        stem = json_path.stem

    # Match RFD3 metrics JSON stem (no .cif_b*_d*): strip outer RF3 summary suffix, then diffusion replica.
    backbone_id = BACKBONE_SUFFIX_RE.sub("", stem)
    backbone_id = re.sub(r"_rf3$", "", backbone_id)
    backbone_id = BACKBONE_SUFFIX_RE.sub("", backbone_id)

    row: dict[str, Any] = {
        "id": stem,
        "backbone_id": backbone_id,
        "filename": f"{stem}.cif",
    }

    chain_ptm = data.get("chain_ptm") or []
    it, ib = _chain_ptm_indices_for_target_binder(data, chain_ptm, target_chain, binder_chain)
    if it is not None and ib is not None:
        row["ptm_target"] = chain_ptm[it] if it < len(chain_ptm) else None
        row["ptm_binder"] = chain_ptm[ib] if ib < len(chain_ptm) else None
    else:
        row["ptm_target"] = chain_ptm[0] if len(chain_ptm) > 0 else None
        row["ptm_binder"] = chain_ptm[1] if len(chain_ptm) > 1 else None

    for src_key, out_key in [
        ("chain_pair_pae_min", "pair_pae_min"),
        ("chain_pair_pde_min", "pair_pde_min"),
        ("chain_pair_pae", "pair_pae"),
        ("chain_pair_pde", "pair_pde"),
    ]:
        matrix = data.get(src_key)
        row[out_key] = _first_off_diagonal_non_null(matrix) if isinstance(matrix, list) else None

    row["plddt"] = data.get("overall_plddt")
    row["overall_pde"] = data.get("overall_pde")
    row["overall_pae"] = data.get("overall_pae")
    row["ptm"] = data.get("ptm")
    row["iptm"] = data.get("iptm")
    row["has_clash"] = data.get("has_clash")
    row["ranking_score"] = data.get("ranking_score")
    row["rf3_ipsae_binder_target"] = data.get("rf3_ipsae_binder_target", data.get("ipsae_binder_target"))
    row["rf3_ipsae_target_binder"] = data.get("rf3_ipsae_target_binder", data.get("ipsae_target_binder"))
    row["rf3_ipsae_min"] = data.get("rf3_ipsae_min", data.get("ipsae_min"))

    return row


def _write_rows(rows: list[dict[str, Any]], output: str) -> int:
    if not rows:
        log.warning("No rows extracted")
        if output != "-":
            Path(output).touch()
        return 0
    df = pd.DataFrame(rows)
    if output == "-":
        df.to_csv(sys.stdout, sep="\t", index=False, lineterminator="\n")
    else:
        df.to_csv(output, sep="\t", index=False, lineterminator="\n")
        log.info("Wrote %s", output)
    return 0


def cmd_rfdiffusion3(args: argparse.Namespace) -> int:
    rows: list[dict[str, Any]] = []
    for jp in args.json:
        try:
            rows.append(extract_rfdiffusion3(jp))
        except Exception as e:
            log.warning("Skipping %s: %s", jp, e)
    return _write_rows(rows, args.output)


def cmd_rf3(args: argparse.Namespace) -> int:
    rows: list[dict[str, Any]] = []
    for jp in args.json:
        try:
            rows.append(extract_rf3(jp, args.target_chain, args.binder_chain))
        except Exception as e:
            log.warning("Skipping %s: %s", jp, e)
    return _write_rows(rows, args.output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract RFD3/RF3 metrics from JSON to TSV")
    subparsers = parser.add_subparsers(dest="command", required=True)

    rfd3_parser = subparsers.add_parser("rfdiffusion3", help="Extract from RFDiffusion3 output JSON(s)")
    rfd3_parser.add_argument("json", type=Path, nargs="+", help="Path(s) to RFD3 output JSON(s)")
    rfd3_parser.add_argument("-o", "--output", default="-", help="Output TSV (default: stdout)")
    rfd3_parser.set_defaults(func=cmd_rfdiffusion3)

    rosettafold3_parser = subparsers.add_parser("rosettafold3", help="Extract from RosettaFold3 summary_confidences JSON(s)")
    rosettafold3_parser.add_argument("json", type=Path, nargs="+", help="Path(s) to *_summary_confidences.json")
    rosettafold3_parser.add_argument("--target-chain", default=None, help="Target chain ID in RF3 output")
    rosettafold3_parser.add_argument("--binder-chain", default=None, help="Binder chain ID in RF3 output")
    rosettafold3_parser.add_argument("-o", "--output", default="-", help="Output TSV (default: stdout)")
    rosettafold3_parser.set_defaults(func=cmd_rf3)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
