#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest", "pyyaml"]
# ///

"""Checks on assets/rfd4_proteina/hydra_overrides, the overlay RFD4_GENERATE stages onto
rfd4-generate's Hydra search path (stdlib + PyYAML; no RFD4-Proteina container)."""

import re
from pathlib import Path

import pytest
import yaml

_REPO = Path(__file__).resolve().parents[2]
_OVERRIDES = _REPO / "assets" / "rfd4_proteina" / "hydra_overrides"
_CUSTOM_METRICS = _OVERRIDES / "metrics" / "custom_metrics" / "all_custom_metrics.yaml"

# The two metrics whose align_mask ignores C_CRDGRP-conditioned targets upstream
_FIXED_METRICS = (
    "refolded_target_aligned_binder_rmsd",
    "refolded_target_backbone_aligned_binder_rmsd",
)


def test_overlay_lands_where_hydra_composes_it() -> None:
    # The filename and group must match `optional custom_metrics: all_custom_metrics`,
    # the last entry of RFD4-Proteina's configs/metrics/all_metrics.yaml defaults list.
    assert _CUSTOM_METRICS.is_file()


def test_package_header_is_the_first_line() -> None:
    # Hydra stops scanning for the header at the first comment that is not a `@` directive,
    # so a `# @package` line further down is silently ignored and the overlay lands in
    # metrics.custom_metrics, where nothing reads it.
    first = _CUSTOM_METRICS.read_text().splitlines()[0]
    assert re.match(r"^#\s*@package\s+metrics\.refolding_metrics\s*$", first), first


@pytest.mark.parametrize("metric", _FIXED_METRICS)
def test_align_mask_accepts_either_coordinate_condition(metric) -> None:
    # A target conditioned with C_CRDGRP is flagged C_INTCRD, not C_CRD. Matching
    # oracle_ppi in configs/metrics/refolding_oracles/af3.yaml, both must be accepted,
    # or the mask selects no atoms and the RMSD is never written.
    masks = yaml.safe_load(_CUSTOM_METRICS.read_text())
    assert "(sa_mask('C_CRD') | sa_mask('C_INTCRD')" in masks[metric]["align_mask"]


@pytest.mark.parametrize("metric", _FIXED_METRICS)
def test_only_align_mask_is_restated(metric) -> None:
    # Everything else (_target_, subset_mask, the tag gates) must come from the base config
    # by deep merge; restating a key here would silently pin it to a stale value.
    masks = yaml.safe_load(_CUSTOM_METRICS.read_text())
    assert list(masks[metric]) == ["align_mask"]


def test_disabled_directory_holds_no_configs() -> None:
    # --rfd4_metrics_config_dir false stages this instead, so adding it to the search path
    # must be a no-op.
    none_dir = _REPO / "assets" / "rfd4_proteina" / "hydra_overrides_none"
    assert none_dir.is_dir()
    assert not list(none_dir.rglob("*.yaml"))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
