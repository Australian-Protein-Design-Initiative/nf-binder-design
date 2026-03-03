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
BACKBONE_SUFFIX_RE = re.compile(r"\.cif_b\d+_d\d+$")


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


def extract_rf3(json_path: Path) -> dict[str, Any]:
    with open(json_path) as f:
        data = json.load(f)

    name = json_path.name
    if name.endswith("_summary_confidences.json"):
        stem = name[: -len("_summary_confidences.json")]
    else:
        stem = json_path.stem

    backbone_id = BACKBONE_SUFFIX_RE.sub("", stem)

    row: dict[str, Any] = {
        "id": stem,
        "backbone_id": backbone_id,
        "filename": f"{stem}.cif",
    }

    # RF3/RFD3 output: chain A = target (index 0), chain B = binder (index 1)
    chain_ptm = data.get("chain_ptm") or []
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


def cmd_rfdiffusion3(args: argparse.Namespace) -> int:
    row = extract_rfdiffusion3(args.json)
    df = pd.DataFrame([row])
    out = args.output
    if out == "-":
        df.to_csv(sys.stdout, sep="\t", index=False, lineterminator="\n")
    else:
        df.to_csv(out, sep="\t", index=False, lineterminator="\n")
        log.info("Wrote %s", out)
    return 0


def cmd_rf3(args: argparse.Namespace) -> int:
    row = extract_rf3(args.json)
    df = pd.DataFrame([row])
    out = args.output
    if out == "-":
        df.to_csv(sys.stdout, sep="\t", index=False, lineterminator="\n")
    else:
        df.to_csv(out, sep="\t", index=False, lineterminator="\n")
        log.info("Wrote %s", out)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract RFD3/RF3 metrics from JSON to TSV")
    subparsers = parser.add_subparsers(dest="command", required=True)

    rfd3_parser = subparsers.add_parser("rfdiffusion3", help="Extract from RFDiffusion3 output JSON")
    rfd3_parser.add_argument("json", type=Path, help="Path to RFD3 output JSON")
    rfd3_parser.add_argument("-o", "--output", default="-", help="Output TSV (default: stdout)")
    rfd3_parser.set_defaults(func=cmd_rfdiffusion3)

    rf3_parser = subparsers.add_parser("rf3", help="Extract from RosettaFold3 summary_confidences JSON")
    rf3_parser.add_argument("json", type=Path, help="Path to *_summary_confidences.json")
    rf3_parser.add_argument("-o", "--output", default="-", help="Output TSV (default: stdout)")
    rf3_parser.set_defaults(func=cmd_rf3)

    rosettafold3_parser = subparsers.add_parser("rosettafold3", help="Extract from RosettaFold3 summary_confidences JSON")
    rosettafold3_parser.add_argument("json", type=Path, help="Path to *_summary_confidences.json")
    rosettafold3_parser.add_argument("-o", "--output", default="-", help="Output TSV (default: stdout)")
    rosettafold3_parser.set_defaults(func=cmd_rf3)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
