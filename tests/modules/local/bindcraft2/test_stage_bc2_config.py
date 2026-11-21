#!/usr/bin/env python
"""Tests for bin/bindcraft2/stage_bc2_config.py."""

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "bin" / "bindcraft2"))

import stage_bc2_config as sbc  # noqa: E402


@pytest.fixture
def campaign_dir(tmp_path: Path) -> Path:
    structures = tmp_path / "structures"
    structures.mkdir()
    (structures / "targetA.pdb").write_text("ATOM\n")
    (structures / "VHH.cif").write_text("data_x\n")
    return tmp_path


def _write(path: Path, config: dict) -> Path:
    path.write_text(json.dumps(config))
    return path


def test_stage_rewrites_nested_paths_to_basenames(campaign_dir: Path) -> None:
    src = _write(
        campaign_dir / "campaign.json",
        {
            "targets": [{"name": "targetA", "target_path": "structures/targetA.pdb"}],
            "binder_scaffold": "structures/VHH.cif",
        },
    )
    out = campaign_dir / "staged.json"
    sbc.stage_config(str(src), str(out))

    staged = json.loads(out.read_text())
    assert staged["targets"][0]["target_path"] == "targetA.pdb"
    assert staged["binder_scaffold"] == "VHH.cif"


def test_stage_leaves_a_preset_only_campaign_alone(campaign_dir: Path) -> None:
    src = _write(campaign_dir / "campaign.json", {"target": "hPDL1", "binder_lengths": [70, 90]})
    out = campaign_dir / "staged.json"
    sbc.stage_config(str(src), str(out))

    assert json.loads(out.read_text()) == {"target": "hPDL1", "binder_lengths": [70, 90]}


def test_extract_input_paths_resolves_against_the_config_directory(campaign_dir: Path) -> None:
    nested = campaign_dir / "nested"
    nested.mkdir()
    src = _write(
        nested / "campaign.json",
        {"targets": [{"target_path": "../structures/targetA.pdb"}]},
    )

    paths = sbc.extract_input_paths(str(src))
    assert len(paths) == 1
    assert Path(paths[0]).resolve() == (campaign_dir / "structures" / "targetA.pdb").resolve()


@pytest.mark.parametrize(
    ("spec", "expected"),
    [("60-150", [60, 150]), ("80", [80]), ("60,80,100", [60, 80, 100])],
)
def test_parse_binder_length_range(spec: str, expected: list) -> None:
    assert sbc.parse_binder_length_range(spec) == expected


def test_generate_config_omits_per_batch_settings(campaign_dir: Path) -> None:
    config = sbc.generate_config(
        design_name="bc2_design",
        input_pdb=str(campaign_dir / "structures" / "targetA.pdb"),
        target_chains="A",
        hotspot_res="A10,A12",
        binder_length_range="60-150",
    )

    assert config["targets"][0]["target_path"] == "targetA.pdb"
    assert config["binder_lengths"] == [60, 150]
    # these are supplied per batch via --set, so they must not be baked in
    for key in ("project_folder", "campaign_seed", "max_trajectories", "number_of_final_designs"):
        assert key not in config


def test_stage_over_a_staged_symlink_does_not_touch_the_source(tmp_path: Path) -> None:
    """The output name may collide with the staged input's.

    Nextflow stages inputs as symlinks. Writing through one would overwrite the
    user's source config outside the work directory - which it did, before
    save_config switched to a temp file plus os.replace.
    """
    source = tmp_path / "source.json"
    source.write_text(json.dumps({"targets": [{"target_path": "input/T.pdb"}]}))

    work = tmp_path / "work"
    work.mkdir()
    staged = work / "campaign.json"
    staged.symlink_to(source)

    sbc.stage_config(str(staged), str(staged))

    assert json.loads(source.read_text())["targets"][0]["target_path"] == "input/T.pdb"
    assert not staged.is_symlink()
    assert json.loads(staged.read_text())["targets"][0]["target_path"] == "T.pdb"


def test_subsample_keeps_ceil_of_the_fraction_in_original_order() -> None:
    kept = sbc.subsample_hotspots("A54,A56,A66,A115", 0.5, seed=1)
    assert kept.split(",") == [h for h in ["A54", "A56", "A66", "A115"] if h in kept.split(",")]
    assert len(kept.split(",")) == 2


def test_subsample_treats_a_range_as_one_item() -> None:
    """A BC2 range like B12-16 is kept or dropped whole, never split."""
    for seed in range(20):
        kept = sbc.subsample_hotspots("A54,B12-16,A66", 0.4, seed=seed).split(",")
        assert len(kept) == 2
        assert set(kept) <= {"A54", "B12-16", "A66"}


def test_subsample_keeps_at_least_one() -> None:
    assert sbc.subsample_hotspots("A54,A56,A66", 0.01, seed=0) in ("A54", "A56", "A66")


def test_subsample_is_reproducible_but_varies_with_the_seed() -> None:
    hotspots = "A1,A2,A3,A4,A5,A6,A7,A8"
    assert sbc.subsample_hotspots(hotspots, 0.5, seed=7) == sbc.subsample_hotspots(
        hotspots, 0.5, seed=7
    )
    draws = {sbc.subsample_hotspots(hotspots, 0.5, seed=s) for s in range(10)}
    assert len(draws) > 1


def test_subsample_at_one_leaves_hotspots_untouched() -> None:
    assert sbc.subsample_hotspots("A54, A56", 1.0, seed=0) == "A54, A56"


def test_stage_subsamples_each_target_differently(tmp_path: Path) -> None:
    """One seed per batch must still give each target its own draw."""
    config = {
        "targets": [
            {"target_path": "a.pdb", "hotspots": "A1,A2,A3,A4"},
            {"target_path": "b.pdb", "hotspots": "A1,A2,A3,A4"},
        ]
    }
    source = _write(tmp_path / "campaign.json", config)
    out = tmp_path / "staged.json"

    # Only a seed that happens to draw differently proves the point, so try a few.
    differed = False
    for seed in range(10):
        sbc.stage_config(str(source), str(out), hotspot_subsample=0.5, subsample_seed=seed)
        targets = json.loads(out.read_text())["targets"]
        assert all(len(t["hotspots"].split(",")) == 2 for t in targets)
        differed = differed or targets[0]["hotspots"] != targets[1]["hotspots"]
    assert differed


def test_stage_without_subsampling_leaves_hotspots_alone(tmp_path: Path) -> None:
    source = _write(
        tmp_path / "campaign.json",
        {"targets": [{"target_path": "x/a.pdb", "hotspots": "A1,A2,A3,A4"}]},
    )
    out = tmp_path / "staged.json"
    sbc.stage_config(str(source), str(out))

    target = json.loads(out.read_text())["targets"][0]
    assert target["hotspots"] == "A1,A2,A3,A4"
    assert target["target_path"] == "a.pdb"
