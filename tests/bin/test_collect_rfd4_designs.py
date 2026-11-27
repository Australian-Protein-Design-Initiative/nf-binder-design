#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Tests for bin/rfd4_proteina/collect_rfd4_designs.py.

    uv run --with pytest pytest tests/bin/test_collect_rfd4_designs.py -q

The inline metrics join happens inside one run directory, because rfd4-generate's example
ids repeat across batches and only the collector knows which batch it is holding.
"""

import csv
import importlib.util
import sys
from pathlib import Path
from typing import List, Optional
from unittest import mock

import pytest

_MODPATH = Path(__file__).resolve().parents[2] / "bin" / "rfd4_proteina" / "collect_rfd4_designs.py"
_spec = importlib.util.spec_from_file_location("collect_rfd4_designs", _MODPATH)
collect = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(collect)

# A cut-down per_example_metrics.csv: Hydra config columns, then the metrics themselves
METRICS_CSV = """\
example_id,designed_structure_path,task,dataset_name,chain_of_interest,ckpt_path,seed,\
input_json_path,generation.args.nsteps,output_base_dir,metrics_tags,\
hotspot_obedience.hotspot_obedience,pre_oracle_pass,oracle_ppi.overall_iptm,\
oracle_ppi.refolded_structure_path,ppi.pass,oracle_apo.apo.overall_plddt,\
oracle_apo.apo.refolded_structure_path,refolded_apo_holo_rmsd.apo.rmsd
pdl1_0,rundir/inference_samples/pdl1_0/pdl1_0_generated.cif,pdl1,nvidia/RFD4-Proteina,B_1,\
model.safetensors,42.0,spec.json,60.0,rundir/inference_samples,"ppi_refolding,apo_refolding",\
0.6078,1.0,89.0,rundir/inference_samples/pdl1_0/pdl1_0_af3_ppi_refolded.cif,0.0,91.2,\
rundir/inference_samples/pdl1_0/pdl1_0_apo_af3_apo_refolded.cif,1.14
pdl1_1,rundir/inference_samples/pdl1_1/pdl1_1_generated.cif,pdl1,nvidia/RFD4-Proteina,B_1,\
model.safetensors,42.0,spec.json,60.0,rundir/inference_samples,"ppi_refolding,apo_refolding",\
0.3922,1.0,79.0,rundir/inference_samples/pdl1_1/pdl1_1_af3_ppi_refolded.cif,1.0,72.5,\
rundir/inference_samples/pdl1_1/pdl1_1_apo_af3_apo_refolded.cif,4.60
"""


def make_run_dir(tmp_path: Path, with_metrics: bool = True) -> Path:
    run_dir = tmp_path / "rundir"
    for example in ("pdl1_0", "pdl1_1"):
        sample_dir = run_dir / "inference_samples" / example
        sample_dir.mkdir(parents=True)
        (sample_dir / f"{example}_generated.cif").write_text("data_x\n")
        (sample_dir / f"{example}_af3_ppi_refolded.cif").write_text("data_ppi\n")
        (sample_dir / f"{example}_apo_af3_apo_refolded.cif").write_text("data_apo\n")
    if with_metrics:
        metrics_dir = run_dir / "inference_metrics"
        metrics_dir.mkdir()
        (metrics_dir / "per_example_metrics.csv").write_text(METRICS_CSV)
    return run_dir


def run_collect(
    run_dir: Path,
    tmp_path: Path,
    batch: str = "3",
    extra_args: Optional[List[str]] = None,
) -> List[dict]:
    tsv = tmp_path / "designs.tsv"
    argv = [
        "collect_rfd4_designs.py",
        str(run_dir),
        "--design-name",
        "pdl1",
        "--batch",
        batch,
        "--seed",
        "42",
        "-o",
        str(tmp_path / "output"),
        "--tsv",
        str(tsv),
    ] + (extra_args or [])
    with mock.patch.object(sys, "argv", argv):
        collect.main()
    with open(tsv, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_metrics_are_merged_onto_the_design_rows(tmp_path: Path) -> None:
    rows = run_collect(make_run_dir(tmp_path), tmp_path)

    assert [r["id"] for r in rows] == ["pdl1_b3_pdl1_0", "pdl1_b3_pdl1_1"]
    assert rows[0]["rfd4_oracle_ppi_overall_iptm"] == "89.0"
    assert rows[1]["rfd4_oracle_ppi_overall_iptm"] == "79.0"
    assert rows[0]["rfd4_hotspot_obedience_hotspot_obedience"] == "0.6078"
    assert rows[1]["rfd4_ppi_pass"] == "1.0"


def test_id_is_the_structure_filename_without_its_extension(tmp_path: Path) -> None:
    """The rfd4_proteina workflow joins the design table to the structures on this."""
    for row in run_collect(make_run_dir(tmp_path), tmp_path):
        assert row["filename"].split(".cif")[0] == row["id"]


def test_run_paths_point_at_the_published_batch(tmp_path: Path) -> None:
    rows = run_collect(make_run_dir(tmp_path), tmp_path, batch="7")
    assert (
        rows[0]["rfd4_oracle_ppi_refolded_structure_path"]
        == "batches/batch7/inference_samples/pdl1_0/pdl1_0_af3_ppi_refolded.cif"
    )


def test_apo_refold_columns_survive_the_merge(tmp_path: Path) -> None:
    """apo_refolding's save_id makes its keys doubly dotted; every dot has to flatten."""
    rows = run_collect(make_run_dir(tmp_path), tmp_path, batch="7")

    assert rows[0]["rfd4_oracle_apo_apo_overall_plddt"] == "91.2"
    # The apo/holo comparison, which only exists because ppi_refolding supplies the holo refold
    assert rows[0]["rfd4_refolded_apo_holo_rmsd_apo_rmsd"] == "1.14"
    assert rows[1]["rfd4_refolded_apo_holo_rmsd_apo_rmsd"] == "4.60"
    assert (
        rows[0]["rfd4_oracle_apo_apo_refolded_structure_path"]
        == "batches/batch7/inference_samples/pdl1_0/pdl1_0_apo_af3_apo_refolded.cif"
    )


def test_hydra_config_and_identity_columns_are_dropped(tmp_path: Path) -> None:
    rows = run_collect(make_run_dir(tmp_path), tmp_path)
    for dropped in ("ckpt_path", "task", "dataset_name", "metrics_tags"):
        assert f"rfd4_{dropped}" not in rows[0]
    assert not any(name.startswith("rfd4_generation") for name in rows[0])
    # The design table's own identity columns are untouched
    assert rows[0]["example_id"] == "pdl1_0"
    assert rows[0]["batch"] == "3"


def test_every_row_in_a_batch_gets_the_same_columns(tmp_path: Path) -> None:
    """collectFile(keepHeader: true) downstream assumes one shared header."""
    run_dir = make_run_dir(tmp_path)
    extra = run_dir / "inference_samples" / "pdl1_2"
    extra.mkdir()
    (extra / "pdl1_2_generated.cif").write_text("data_x\n")

    rows = run_collect(run_dir, tmp_path)
    assert len(rows) == 3
    assert set(rows[2]) == set(rows[0])
    # The design with no metrics row gets empty values, not missing columns
    assert rows[2]["rfd4_oracle_ppi_overall_iptm"] == ""


def test_designs_are_collected_without_metrics(tmp_path: Path) -> None:
    rows = run_collect(make_run_dir(tmp_path, with_metrics=False), tmp_path)
    assert [r["id"] for r in rows] == ["pdl1_b3_pdl1_0", "pdl1_b3_pdl1_1"]
    assert not any(name.startswith("rfd4_") for name in rows[0])


def test_a_metrics_csv_without_example_id_is_rejected(tmp_path: Path) -> None:
    run_dir = make_run_dir(tmp_path)
    (run_dir / "inference_metrics" / "per_example_metrics.csv").write_text("a,b\n1,2\n")
    with pytest.raises(SystemExit, match="no example_id column"):
        run_collect(run_dir, tmp_path)


def test_refolds_are_named_after_the_design_they_refold(tmp_path: Path) -> None:
    refolded = tmp_path / "refolded"
    rows = run_collect(
        make_run_dir(tmp_path), tmp_path, batch="0", extra_args=["--refolded-dir", str(refolded)]
    )

    # Each name is the design id plus the suffix naming the oracle run
    design_ids = [r["id"] for r in rows]
    assert sorted(p.name for p in refolded.iterdir()) == sorted(
        f"{design_id}{suffix}"
        for design_id in design_ids
        for suffix in ("_af3_ppi_refolded.cif", "_apo_af3_apo_refolded.cif")
    )
    assert design_ids == ["pdl1_b0_pdl1_0", "pdl1_b0_pdl1_1"]
    # The generated designs themselves are not swept up as refolds
    assert (refolded / "pdl1_b0_pdl1_0_af3_ppi_refolded.cif").read_text() == "data_ppi\n"


def test_run_paths_gain_a_gz_suffix_when_the_batch_is_gzipped(tmp_path: Path) -> None:
    rows = run_collect(
        make_run_dir(tmp_path), tmp_path, batch="7", extra_args=["--run-cifs-gzipped"]
    )
    assert (
        rows[0]["rfd4_oracle_ppi_refolded_structure_path"]
        == "batches/batch7/inference_samples/pdl1_0/pdl1_0_af3_ppi_refolded.cif.gz"
    )
