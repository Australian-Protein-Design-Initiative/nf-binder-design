#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# ///

"""After rf3 fold on a multi-example JSON, map each input name to outputs, run IPSAE + score extract, stage per_design/<meta_id>/."""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
log = logging.getLogger(__name__)

def _assign_summaries(output_root: Path, names_in_order: list[str]) -> dict[str, Path]:
    all_summaries = sorted(output_root.rglob("*_summary_confidences.json"))
    assigned: dict[str, Path] = {}
    used: set[Path] = set()
    for name in sorted(names_in_order, key=len, reverse=True):
        for p in all_summaries:
            if p in used:
                continue
            stem = p.name[: -len("_summary_confidences.json")]
            if stem == name or stem.startswith(name + "_"):
                assigned[name] = p
                used.add(p)
                break
    missing = [n for n in names_in_order if n not in assigned]
    if missing:
        raise SystemExit(f"No summary_confidences JSON matched RF3 input name(s): {missing!r} under {output_root}")
    return assigned


def _find_full_confidences(outdir: Path) -> Optional[Path]:
    for p in sorted(outdir.glob("*_confidences.json")):
        if "summary" not in p.name:
            return p
    return None


def _find_model_cif(outdir: Path) -> Optional[Path]:
    for p in sorted(outdir.glob("*_model.cif")):
        return p
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rf3-input-json", type=Path, required=True, help="rf3 fold inputs JSON (list of configs)")
    parser.add_argument("--rf3-output-dir", type=Path, required=True, help="rf3 fold out_dir (e.g. output/)")
    parser.add_argument("--target-chain", required=True)
    parser.add_argument("--binder-chain", required=True)
    parser.add_argument("--project-dir", type=Path, required=True, help="Pipeline root (ipsae, extract_rfd3_scores)")
    parser.add_argument("--batch-scores-out", type=Path, default=None,
                        help="If set, concatenate all per-design score TSVs into this file")
    parser.add_argument("--work-cwd", type=Path, default=Path("."), help="Where to write per-design dirs and scores")
    args = parser.parse_args()

    data = json.loads(args.rf3_input_json.read_text())
    if not isinstance(data, list):
        raise SystemExit("RF3 inputs JSON must be a list")
    names_in_order = [str(item.get("name", "")).strip() for item in data]
    if not all(names_in_order):
        raise SystemExit("Each RF3 input must have a non-empty name")
    name_to_summary = _assign_summaries(args.rf3_output_dir, names_in_order)

    per_design = args.work_cwd / "per_design"
    per_design.mkdir(parents=True, exist_ok=True)

    all_summaries: list[Path] = []

    for name in names_in_order:
        summary_path = name_to_summary[name]
        outdir = summary_path.parent
        if name.endswith("_rf3"):
            meta_id = name[: -len("_rf3")]
        else:
            meta_id = name

        full_conf = _find_full_confidences(outdir)
        model_cif = _find_model_cif(outdir)
        if not full_conf or not model_cif:
            raise SystemExit(f"Missing full confidences or model CIF in {outdir}")

        subprocess.run(
            [
                sys.executable,
                str(args.project_dir / "bin/ipsae.py"),
                "--format",
                "rf3",
                "--update-summary",
                summary_path.name,
                "--binder-chain",
                args.binder_chain,
                "--target-chain",
                args.target_chain,
                full_conf.name,
                model_cif.name,
                "10",
                "10",
            ],
            cwd=outdir,
            check=True,
        )

        all_summaries.append(summary_path)

        dest_dir = per_design / meta_id
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest_model = dest_dir / "model.cif"
        if dest_model.exists() or dest_model.is_symlink():
            dest_model.unlink()
        try:
            dest_model.symlink_to(model_cif.resolve())
        except OSError:
            shutil.copy2(model_cif, dest_model)

    scores_out = args.batch_scores_out or (args.work_cwd / "rf3_batch_scores.tsv")
    subprocess.run(
        [
            sys.executable,
            str(args.project_dir / "bin/rfd3/extract_rfd3_scores.py"),
            "rosettafold3",
            *[str(p) for p in all_summaries],
            "--target-chain",
            args.target_chain,
            "--binder-chain",
            args.binder_chain,
            "-o",
            str(scores_out),
        ],
        check=True,
    )

    # Write per-design scores.tsv (header + that design's row only)
    batch_lines = scores_out.read_text().splitlines()
    if batch_lines:
        header = batch_lines[0]
        cols = header.split("\t")
        id_col = cols.index("id") if "id" in cols else 0
        row_by_name: dict[str, str] = {}
        for line in batch_lines[1:]:
            fields = line.split("\t")
            row_by_name[fields[id_col]] = line

        for name in names_in_order:
            if name.endswith("_rf3"):
                meta_id = name[: -len("_rf3")]
            else:
                meta_id = name
            data_line = row_by_name.get(name, "")
            dest_scores = per_design / meta_id / "scores.tsv"
            dest_scores.write_text(header + "\n" + data_line + "\n")
            log.info("Staged %s -> per_design/%s/model.cif / scores.tsv", meta_id, meta_id)

    return 0


if __name__ == "__main__":
    sys.exit(main())
