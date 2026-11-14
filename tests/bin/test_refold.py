#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "biopython"]
# ///

"""
Tests for the refold helpers: bin/refold/prepare_inputs.py, bin/refold/rmsd_rows.py,
bin/refold/summary_to_wide.py and bin/refold/sample_scores.py.

    uv run --with pytest --with biopython pytest tests/bin/test_refold.py -q
"""

import csv
import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("Bio")

_BIN = Path(__file__).resolve().parents[2] / "bin" / "refold"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare = _load("prepare_inputs")
to_wide = _load("summary_to_wide")
rmsd_rows = _load("rmsd_rows")
sample_scores = _load("sample_scores")


# Target chain A (2 residues), binder chain B (2 residues)
def _pdb(binder_seq: str) -> str:
    binder_3 = {"A": "ALA", "G": "GLY", "S": "SER"}
    lines = [
        "ATOM      1  CA  ALA A  18      32.394 -38.055 -22.849  1.00 37.28           C",
        "ATOM      2  CA  GLY A  19      31.602 -39.292 -22.650  1.00 35.59           C",
    ]
    for i, residue in enumerate(binder_seq):
        lines.append(
            f"ATOM  {3 + i:5d}  CA  {binder_3[residue]} B{1 + i:4d}"
            f"      12.000  -8.000  -{2 + i}.000  1.00 30.00           C"
        )
    return "\n".join(lines) + "\nEND\n"


@pytest.fixture()
def designs(tmp_path: Path) -> Path:
    d = tmp_path / "designs"
    d.mkdir()
    (d / "rfd4_x_0.pdb").write_text(_pdb("AG"))
    (d / "rfd4_x_1.pdb").write_text(_pdb("SA"))
    return d


def test_binder_record_per_design_named_after_the_file(designs: Path, tmp_path: Path):
    sequences = prepare.collect_sequences(sorted(designs.glob("*.pdb")), "B", "A")
    assert set(sequences) == {"rfd4_x_0", "rfd4_x_1"}
    assert sequences["rfd4_x_0"]["binder"] == "AG"
    assert sequences["rfd4_x_1"]["binder"] == "SA"


def test_target_is_shared_across_designs(designs: Path):
    sequences = prepare.collect_sequences(sorted(designs.glob("*.pdb")), "B", "A")
    assert {v["target"] for v in sequences.values()} == {"AG"}
    assert prepare.consensus_target(sequences) == "AG"


def test_consensus_target_takes_the_majority_and_warns(caplog):
    sequences = {
        "a": {"binder": "A", "target": "AG"},
        "b": {"binder": "A", "target": "AG"},
        "c": {"binder": "A", "target": "AGS"},
    }
    assert prepare.consensus_target(sequences) == "AG"
    assert "disagree on the target sequence" in caplog.text


def test_design_id_strips_compound_suffixes():
    assert prepare._design_id(Path("x/rfd4_a_0.cif.gz")) == "rfd4_a_0"
    assert prepare._design_id(Path("x/rfd4_a_0.cif")) == "rfd4_a_0"
    assert prepare._design_id(Path("x/rfd4_a_0.pdb")) == "rfd4_a_0"


def test_a_design_missing_the_binder_chain_is_skipped(tmp_path: Path, caplog):
    only_target = tmp_path / "no_binder.pdb"
    only_target.write_text(
        "ATOM      1  CA  ALA A  18      32.394 -38.055 -22.849  1.00 37.28           C\nEND\n"
    )
    assert prepare.collect_sequences([only_target], "B", "A") == {}
    assert "no chain B" in caplog.text


SUMMARY_ROWS = [
    {
        "target": "target", "binder": "rfd4_x_0", "tool": "boltz", "n": "2",
        "iptm_max": "0.78", "ipsae_max": "0.55",
        "consensus_z": "1.42", "consensus_z_metric": "ipsae", "n_tools": "2",
    },
    {
        "target": "target", "binder": "rfd4_x_0", "tool": "chai1", "n": "2",
        "iptm_max": "0.81", "ipsae_max": "0.62",
        "consensus_z": "1.42", "consensus_z_metric": "ipsae", "n_tools": "2",
    },
    {
        "target": "target", "binder": "rfd4_x_1", "tool": "boltz", "n": "2",
        "iptm_max": "0.40", "ipsae_max": "0.11",
        "consensus_z": "-0.9", "consensus_z_metric": "ipsae", "n_tools": "1",
    },
]


def test_pivot_gives_one_row_per_design_with_tool_prefixed_columns():
    by_id, columns = to_wide.pivot(SUMMARY_ROWS, ["iptm_max", "ipsae_max"])
    assert set(by_id) == {"rfd4_x_0", "rfd4_x_1"}
    assert by_id["rfd4_x_0"]["chai1_iptm_max"] == "0.81"
    assert by_id["rfd4_x_0"]["boltz_ipsae_max"] == "0.55"
    # A design scored by only one tool leaves the other tool's columns blank rather
    # than carrying the first tool's value across.
    assert "chai1_iptm_max" not in by_id["rfd4_x_1"]
    assert "chai1_iptm_max" in columns


def test_pair_level_columns_are_written_once_not_per_tool():
    by_id, columns = to_wide.pivot(SUMMARY_ROWS, ["iptm_max"])
    assert "consensus_z" in columns
    assert not any(c.startswith("boltz_consensus_z") for c in columns)
    assert by_id["rfd4_x_0"]["consensus_z"] == "1.42"


def test_several_targets_keep_their_columns_apart():
    rows = SUMMARY_ROWS + [
        {
            "target": "target2", "binder": "rfd4_x_0", "tool": "boltz", "n": "2",
            "iptm_max": "0.20", "ipsae_max": "0.05",
            "consensus_z": "0.1", "consensus_z_metric": "ipsae", "n_tools": "1",
        }
    ]
    pivoted, columns = to_wide.pivot(rows, ["iptm_max"])
    assert "target_boltz_iptm_max" in columns
    assert "target2_boltz_iptm_max" in columns
    row = pivoted["rfd4_x_0"]
    assert row["target_boltz_iptm_max"] == "0.78"
    assert row["target2_boltz_iptm_max"] == "0.20"


def test_empty_summary_still_writes_a_header(tmp_path: Path, monkeypatch):
    summary = tmp_path / "summary.tsv"
    summary.write_text("target\tbinder\ttool\n")
    out = tmp_path / "out.tsv"
    monkeypatch.setattr("sys.argv", ["x", str(summary), "-o", str(out)])
    assert to_wide.main() == 0
    assert out.read_text() == "id\n"


def test_end_to_end_writes_both_fastas(designs: Path, tmp_path: Path, monkeypatch):
    targets = tmp_path / "targets.fasta"
    binders = tmp_path / "binders.fasta"
    monkeypatch.setattr(
        "sys.argv",
        [
            "x", str(designs),
            "--binder-chain", "B", "--target-chain", "A",
            "--targets-out", str(targets), "--binders-out", str(binders),
        ],
    )
    assert prepare.main() == 0
    assert targets.read_text() == ">target\nAG\n"
    assert binders.read_text() == ">rfd4_x_0\nAG\n>rfd4_x_1\nSA\n"


def test_wide_output_is_valid_tsv(tmp_path: Path, monkeypatch):
    summary = tmp_path / "summary.tsv"
    with open(summary, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(SUMMARY_ROWS[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(SUMMARY_ROWS)
    out = tmp_path / "wide.tsv"
    monkeypatch.setattr("sys.argv", ["x", str(summary), "-o", str(out)])
    assert to_wide.main() == 0
    with open(out, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["id"] for r in rows] == ["rfd4_x_0", "rfd4_x_1"]


# --- rmsd_rows: rmsd4all tables -> long rows keyed by (id, tool, model) ---

def _rmsd_table(path: Path, rows):
    header = (
        "structure1\tstructure2\trmsd_pruned\tn_pairs_rmsd_pruned\trmsd_all\t"
        "n_pairs_rmsd_all\ttm_score_pruned\tn_pairs_tm_score_pruned\ttm_score_all\t"
        "n_pairs_tm_score_all\n"
    )
    path.write_text(header + "".join(rows))


def test_rmsd_rows_joins_the_tables_on_the_model(tmp_path: Path, monkeypatch):
    aligned = tmp_path / "aligned.tsv"
    complexes = tmp_path / "complex.tsv"
    _rmsd_table(aligned, ["d0.cif\td0_model_0.cif\t0.94\t50\t1.02\t60\t0.87\t60\t0.87\t60\n"])
    _rmsd_table(complexes, ["d0.cif\td0_model_0.cif\t1.30\t50\t1.41\t60\t0.80\t60\t0.80\t60\n"])
    out = tmp_path / "rows.tsv"
    monkeypatch.setattr("sys.argv", [
        "x", "--id", "d0", "--tool", "boltz", "--batch", "4",
        "--table", f"target_aligned_binder={aligned}",
        "--table", f"complex={complexes}",
        "-o", str(out),
    ])
    assert rmsd_rows.main() == 0
    row = list(csv.DictReader(open(out, newline=""), delimiter="\t"))[0]
    assert row["id"] == "d0"
    assert row["tool"] == "boltz"
    assert row["batch"] == "4"
    assert row["model"] == "0"
    assert row["target_aligned_binder_rmsd_all"] == "1.02"
    assert row["complex_rmsd_all"] == "1.41"


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("d0_model_3.cif", "3"),
        ("d0_seed_42_sample_2_model.cif", "2"),
        ("d0_sample-1_model.cif", "1"),
        ("relaxed_model_1_multimer.cif", "7"),  # no sample index: falls back to position
    ],
)
def test_model_index_reads_the_sample_number_each_engine_writes(filename, expected):
    assert rmsd_rows.model_index(filename, 7) == expected


# --- summary_to_wide: long tables merged in beside the confidence summary ---

def test_long_tables_become_tool_prefixed_columns():
    by_design, columns = to_wide.pivot(SUMMARY_ROWS, ["iptm_max"])
    to_wide.add_long_rows(by_design, columns, [
        {"id": "rfd4_x_0", "tool": "boltz", "model": "0", "complex_rmsd_all": "1.41"},
        {"id": "rfd4_x_0", "tool": "boltz", "model": "1", "complex_rmsd_all": "0.90"},
    ], "min")
    assert by_design["rfd4_x_0"]["boltz_complex_rmsd_all"] == "0.9"
    assert "boltz_complex_rmsd_all" in columns


def test_rmsd_stat_chooses_which_model_speaks_for_the_design():
    rows = [
        {"id": "d0", "tool": "boltz", "model": "0", "complex_rmsd_all": "2.0"},
        {"id": "d0", "tool": "boltz", "model": "1", "complex_rmsd_all": "1.0"},
    ]
    for stat, expected in [("min", "1"), ("max", "2"), ("first", "2"), ("mean", "1.5")]:
        by_design, columns = {}, []
        to_wide.add_long_rows(by_design, columns, rows, stat)
        assert by_design["d0"]["boltz_complex_rmsd_all"] == expected


def test_interface_scores_are_attributed_by_their_filename(tmp_path: Path):
    scores = tmp_path / "interface.tsv"
    scores.write_text(
        "filename\tinterface_dG\nrfd4_x_0__boltz__b3__seed_44_sample_0.pdb\t-44.7\n"
    )
    rows = to_wide.rows_from_bindcraft(str(scores))
    assert rows == [{
        "id": "rfd4_x_0", "tool": "boltz", "batch": "3",
        "model": "seed_44_sample_0", "interface_dG": "-44.7",
    }]


@pytest.mark.parametrize(
    "filename,expected",
    [
        # The current shape, and the two earlier ones, so an older score table merges.
        ("rfd4_x_0__boltz__b3__0.pdb", ("rfd4_x_0", "boltz", "3", "0")),
        ("rfd4_x_0__boltz__b3.pdb", ("rfd4_x_0", "boltz", "3", "")),
        ("rfd4_x_0__boltz.pdb", ("rfd4_x_0", "boltz", "1", "")),
        ("something_else.pdb", None),
    ],
)
def test_the_interface_filename_shapes_that_are_read(filename, expected):
    assert to_wide._parse_interface_filename(filename) == expected


def test_first_means_the_first_batch_not_the_first_row_written(tmp_path: Path):
    # collectFile gives no ordering guarantee, so the rows can arrive in any order.
    rows = [
        {"id": "d0", "tool": "boltz", "batch": "3", "model": "0", "interface_dG": "-30.0"},
        {"id": "d0", "tool": "boltz", "batch": "1", "model": "0", "interface_dG": "-44.0"},
    ]
    by_design, columns = {}, []
    to_wide.add_long_rows(by_design, columns, rows, "first")
    assert by_design["d0"]["boltz_interface_dG"] == "-44"


def test_the_batch_is_a_key_not_a_metric():
    rows = [
        {"id": "d0", "tool": "boltz", "batch": "2", "model": "0", "complex_rmsd_all": "1.0"},
    ]
    by_design, columns = {}, []
    to_wide.add_long_rows(by_design, columns, rows, "min")
    assert "boltz_batch" not in columns
    assert by_design["d0"]["boltz_complex_rmsd_all"] == "1"


def test_an_interface_score_that_cannot_name_its_design_is_an_error(tmp_path: Path):
    scores = tmp_path / "interface.tsv"
    scores.write_text("filename\tinterface_dG\nsomething_else.pdb\t-44.7\n")
    with pytest.raises(SystemExit):
        to_wide.rows_from_bindcraft(str(scores))


# --- sample_scores: one row per refolded sample, nothing reduced ---

def _sample_inputs(tmp_path: Path):
    scores = tmp_path / "fold_pulldown_scores.tsv"
    scores.write_text(
        "id\ttarget\tbinder\ttool\tmodel\tbatch\tiptm\tipsae\n"
        "target_and_d0\ttarget\td0\tboltz\t0\t1\t0.95\t0.89\n"
        "target_and_d0\ttarget\td0\tboltz\t0\t2\t0.91\t0.84\n"
    )
    rmsd = tmp_path / "refold_rmsd.tsv"
    rmsd.write_text(
        "id\ttool\tbatch\tmodel\tcomplex_rmsd_all\n"
        "d0\tboltz\t1\t0\t1.27\n"
        "d0\tboltz\t2\t0\t1.44\n"
    )
    interface = tmp_path / "refold_interface_scores.tsv"
    interface.write_text(
        "filename\tinterface_dG\n"
        "d0__boltz__b1.pdb\t-43.55\n"
        "d0__boltz__b2.pdb\t-41.02\n"
    )
    return scores, rmsd, interface


def _run_sample_scores(tmp_path: Path, monkeypatch, scores, rmsd, interface):
    out = tmp_path / "refold_sample_scores.tsv"
    monkeypatch.setattr("sys.argv", [
        "x", str(scores), "--long", str(rmsd), "--bindcraft", str(interface), "-o", str(out),
    ])
    assert sample_scores.main() == 0
    return list(csv.DictReader(open(out, newline=""), delimiter="\t"))


def test_every_sample_keeps_its_own_scores(tmp_path: Path, monkeypatch):
    rows = _run_sample_scores(tmp_path, monkeypatch, *_sample_inputs(tmp_path))
    assert [r["batch"] for r in rows] == ["1", "2"]
    assert [r["iptm"] for r in rows] == ["0.95", "0.91"]
    assert [r["complex_rmsd_all"] for r in rows] == ["1.27", "1.44"]
    assert [r["interface_dG"] for r in rows] == ["-43.55", "-41.02"]


def test_the_design_is_the_id_and_the_complex_is_the_pair_id(tmp_path: Path, monkeypatch):
    # The score tables call the complex `id` and the design `binder`; every other
    # refold table calls the design `id`, and this one follows those.
    rows = _run_sample_scores(tmp_path, monkeypatch, *_sample_inputs(tmp_path))
    assert rows[0]["id"] == "d0"
    assert rows[0]["pair_id"] == "target_and_d0"
    assert rows[0]["target"] == "target"


def test_a_sample_with_no_rmsd_row_keeps_its_confidence(tmp_path: Path, monkeypatch):
    scores, rmsd, interface = _sample_inputs(tmp_path)
    rmsd.write_text("id\ttool\tbatch\tmodel\tcomplex_rmsd_all\nd0\tboltz\t1\t0\t1.27\n")
    rows = _run_sample_scores(tmp_path, monkeypatch, scores, rmsd, interface)
    assert rows[1]["iptm"] == "0.91"
    assert rows[1]["complex_rmsd_all"] == ""


def test_a_model_index_that_disagrees_still_joins_within_the_batch(tmp_path: Path, monkeypatch):
    # The model index is parsed out of a filename and the engines do not all number
    # their samples the same way, so a batch holding one model joins on the batch.
    scores, rmsd, interface = _sample_inputs(tmp_path)
    rmsd.write_text(
        "id\ttool\tbatch\tmodel\tcomplex_rmsd_all\n"
        "d0\tboltz\t1\t3\t1.27\n"
        "d0\tboltz\t2\t3\t1.44\n"
    )
    rows = _run_sample_scores(tmp_path, monkeypatch, scores, rmsd, interface)
    assert [r["complex_rmsd_all"] for r in rows] == ["1.27", "1.44"]


def test_an_empty_score_table_still_writes_a_header(tmp_path: Path, monkeypatch):
    empty = tmp_path / "empty_scores.tsv"
    empty.write_text("id\ttarget\tbinder\ttool\tmodel\tbatch\n")
    out = tmp_path / "out.tsv"
    monkeypatch.setattr("sys.argv", ["x", str(empty), "-o", str(out)])
    assert sample_scores.main() == 0
    assert out.read_text().startswith("id\ttool\tbatch\tmodel")


# --- summary_to_wide: the structural columns all come from one sample ---

def _sample_table(tmp_path: Path, rows: str) -> Path:
    path = tmp_path / "refold_sample_scores.tsv"
    path.write_text("id\ttool\tbatch\tmodel\tipsae\tiptm\n" + rows)
    return path


def test_the_best_ranked_sample_is_the_one_selected(tmp_path: Path):
    samples = _sample_table(tmp_path, (
        "d0\tboltz\t1\t0\t0.70\t0.80\n"
        "d0\tboltz\t2\t0\t0.91\t0.75\n"
        "d0\tboltz\t3\t0\t0.60\t0.99\n"
    ))
    assert to_wide.selected_samples(str(samples), "ipsae") == {("d0", "boltz"): ("2", "0")}


def test_the_rank_metric_falls_back_when_the_engine_has_none_of_it(tmp_path: Path):
    samples = _sample_table(tmp_path, (
        "d0\tboltz\t1\t0\t\t0.80\n"
        "d0\tboltz\t2\t0\t\t0.99\n"
    ))
    assert to_wide.selected_samples(str(samples), "ipsae") == {("d0", "boltz"): ("2", "0")}


def test_a_tie_and_a_missing_metric_both_pick_the_first_sample(tmp_path: Path):
    samples = _sample_table(tmp_path, (
        "d0\tboltz\t2\t0\t0.80\t0.80\n"
        "d0\tboltz\t1\t0\t0.80\t0.80\n"
        "d1\tboltz\t2\t0\t\t\n"
        "d1\tboltz\t1\t0\t\t\n"
    ))
    chosen = to_wide.selected_samples(str(samples), "ipsae")
    assert chosen == {("d0", "boltz"): ("1", "0"), ("d1", "boltz"): ("1", "0")}


def test_rmsd_and_interface_come_from_the_same_sample(tmp_path: Path):
    # Batch 2 ranks best on ipsae but is not the closest to the design, so a
    # per-metric statistic would report batch 3's RMSD beside batch 2's interface.
    samples = _sample_table(tmp_path, (
        "d0\tboltz\t1\t0\t0.70\t0.80\n"
        "d0\tboltz\t2\t0\t0.91\t0.75\n"
        "d0\tboltz\t3\t0\t0.60\t0.99\n"
    ))
    selected = to_wide.selected_samples(str(samples), "ipsae")
    rmsd = [
        {"id": "d0", "tool": "boltz", "batch": "1", "model": "0", "complex_rmsd_all": "2.0"},
        {"id": "d0", "tool": "boltz", "batch": "2", "model": "0", "complex_rmsd_all": "1.5"},
        {"id": "d0", "tool": "boltz", "batch": "3", "model": "0", "complex_rmsd_all": "0.9"},
    ]
    interface = [
        {"id": "d0", "tool": "boltz", "batch": "1", "model": "0", "interface_dG": "-30.0"},
        {"id": "d0", "tool": "boltz", "batch": "2", "model": "0", "interface_dG": "-41.0"},
        {"id": "d0", "tool": "boltz", "batch": "3", "model": "0", "interface_dG": "-55.0"},
    ]
    by_design, columns = {}, []
    to_wide.add_long_rows(by_design, columns, rmsd, "selected", selected)
    to_wide.add_long_rows(by_design, columns, interface, "selected", selected)
    assert by_design["d0"]["boltz_complex_rmsd_all"] == "1.5"
    assert by_design["d0"]["boltz_interface_dG"] == "-41"


def test_samples_sharing_a_batch_are_told_apart_by_the_model(tmp_path: Path):
    # An engine handed the whole seed list returns every sample from batch 1, so the
    # batch cannot identify one: without the model label this picked seed 42's row.
    samples = _sample_table(tmp_path, (
        "d0\tesmfold2_fast\t1\tseed_42_sample_0\t0.70\t0.80\n"
        "d0\tesmfold2_fast\t1\tseed_43_sample_0\t0.91\t0.75\n"
        "d0\tesmfold2_fast\t1\tseed_44_sample_0\t0.60\t0.99\n"
    ))
    selected = to_wide.selected_samples(str(samples), "ipsae")
    assert selected == {("d0", "esmfold2_fast"): ("1", "seed_43_sample_0")}

    rows = [
        {"id": "d0", "tool": "esmfold2_fast", "batch": "1", "model": "seed_42_sample_0",
         "complex_rmsd_all": "2.0"},
        {"id": "d0", "tool": "esmfold2_fast", "batch": "1", "model": "seed_43_sample_0",
         "complex_rmsd_all": "1.5"},
        {"id": "d0", "tool": "esmfold2_fast", "batch": "1", "model": "seed_44_sample_0",
         "complex_rmsd_all": "0.9"},
    ]
    by_design, columns = {}, []
    to_wide.add_long_rows(by_design, columns, rows, "selected", selected)
    assert by_design["d0"]["esmfold2_fast_complex_rmsd_all"] == "1.5"


def test_a_table_with_no_model_column_still_matches_on_the_batch(tmp_path: Path):
    samples = _sample_table(tmp_path, (
        "d0\tboltz\t1\t0\t0.70\t0.80\n"
        "d0\tboltz\t2\t0\t0.91\t0.75\n"
    ))
    selected = to_wide.selected_samples(str(samples), "ipsae")
    rows = [
        {"id": "d0", "tool": "boltz", "batch": "1", "interface_dG": "-30.0"},
        {"id": "d0", "tool": "boltz", "batch": "2", "interface_dG": "-41.0"},
    ]
    by_design, columns = {}, []
    to_wide.add_long_rows(by_design, columns, rows, "selected", selected)
    assert by_design["d0"]["boltz_interface_dG"] == "-41"


def test_selected_without_a_selection_falls_back_to_the_first_sample():
    rows = [
        {"id": "d0", "tool": "boltz", "batch": "2", "model": "0", "complex_rmsd_all": "0.9"},
        {"id": "d0", "tool": "boltz", "batch": "1", "model": "0", "complex_rmsd_all": "2.0"},
    ]
    by_design, columns = {}, []
    to_wide.add_long_rows(by_design, columns, rows, "selected", {})
    assert by_design["d0"]["boltz_complex_rmsd_all"] == "2"
