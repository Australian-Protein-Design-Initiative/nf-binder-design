#!/usr/bin/env python
# /// script
# requires-python = ">=3.9"
# dependencies = ["pytest"]
# ///

"""
Every process used by fold / fold_pulldown must have a withName selector in
nextflow.config, conf/platforms/m3.config and conf/platforms/m3_bdi.config
(AGENTS.md close-out checklist), and a container selector in
conf/platforms/monash_containers.config if its module declares a container.

Nextflow matches a selector as a regex against the WHOLE simple name (declared
name or include alias) or the whole fully-qualified name, so partial paths such
as 'BOLTZ_FOLD:X' match nothing once the subworkflow is nested under FOLD. The
checks below use the same whole-name rule.

    uv run --with pytest pytest tests/bin/test_withname_coverage.py -q
"""

import re
from pathlib import Path
from typing import Dict, List, Set

import pytest

ROOT = Path(__file__).resolve().parents[2]
PLATFORMS = ROOT / "conf" / "platforms"
CONFIGS = [ROOT / "nextflow.config", PLATFORMS / "m3.config", PLATFORMS / "m3_bdi.config"]
CONTAINER_CONFIG = PLATFORMS / "monash_containers.config"
# Declared containers with no Monash mirror image (see the header of monash_containers.config)
UNMIRRORED = {"MMSEQS_COLABFOLDSEARCH"}

_PROCESS_RE = re.compile(r"^\s*process\s+([A-Z0-9_]+)\s*\{", re.M)
_INCLUDE_RE = re.compile(r"include\s*\{([^}]*)\}\s*from\s*'([^']+)'")
_SELECTOR_RE = re.compile(r"withName:\s*(?:'([^']+)'|\"([^\"]+)\"|/([^/]+)/|([A-Za-z0-9_]+))")


def _selectors(path: Path) -> List[str]:
    text = re.sub(r"//[^\n]*", "", path.read_text())
    return [next(g for g in m.groups() if g) for m in _SELECTOR_RE.finditer(text)]


def _processes() -> Dict[str, Set[str]]:
    """Declared process name -> names it is included as, for fold/fold_pulldown."""
    declared: Dict[Path, List[str]] = {}
    for nf in (ROOT / "modules").rglob("*.nf"):
        declared[nf.resolve()] = _PROCESS_RE.findall(nf.read_text())

    entry = [ROOT / "workflows" / "fold.nf", ROOT / "workflows" / "fold_pulldown.nf"]
    seen: Set[Path] = set()
    names: Dict[str, Set[str]] = {}
    stack = [p.resolve() for p in entry]
    while stack:
        nf = stack.pop()
        if nf in seen or not nf.exists():
            continue
        seen.add(nf)
        for body, src in _INCLUDE_RE.findall(nf.read_text()):
            target = (nf.parent / src).resolve()
            target = target if target.suffix == ".nf" else target.with_suffix(".nf")
            for item in body.split(";"):
                parts = item.split(" as ")
                base = parts[0].strip()
                alias = parts[-1].strip()
                if base in declared.get(target, []):
                    names.setdefault(base, set()).update({base, alias})
            stack.append(target)
    return names


def _containerised() -> Set[str]:
    """Declared names of processes whose module has a `container` directive."""
    found: Set[str] = set()
    for nf in (ROOT / "modules").rglob("*.nf"):
        parts = re.split(r"^\s*process\s+([A-Z0-9_]+)\s*\{", nf.read_text(), flags=re.M)
        # parts = [preamble, name1, body1, name2, body2, ...]
        for name, body in zip(parts[1::2], parts[2::2]):
            if re.search(r"^\s*container\b", body, re.M):
                found.add(name)
    return found


PROCESSES = _processes()


def test_found_fold_processes():
    assert {"OPENFOLD3", "ANNOTATE_MSA", "FOLD_SCORE_AF2", "SPLIT_COMPLEX_FASTA"} <= set(PROCESSES)


def _missing_selectors(config: Path, processes: Dict[str, Set[str]]) -> List[str]:
    selectors = [re.compile(s) for s in _selectors(config)]
    return sorted(
        base for base, aliases in processes.items()
        # 'SUBWF:' stands in for any qualified path, for selectors like '.*:PROC'
        if not any(sel.fullmatch(name) or sel.fullmatch(f"SUBWF:{name}")
                   for sel in selectors for name in aliases)
    )


@pytest.mark.parametrize("config", CONFIGS, ids=lambda p: p.name)
def test_every_fold_process_has_a_selector(config):
    missing = _missing_selectors(config, PROCESSES)
    assert not missing, f"{config.name}: no withName selector for {missing}"


def test_every_containerised_fold_process_has_a_mirror_selector():
    containerised = _containerised() - UNMIRRORED
    needed = {base: aliases for base, aliases in PROCESSES.items() if base in containerised}
    assert needed, "no containerised fold processes found"
    missing = _missing_selectors(CONTAINER_CONFIG, needed)
    assert not missing, f"{CONTAINER_CONFIG.name}: no container selector for {missing}"


@pytest.mark.parametrize("config", CONFIGS, ids=lambda p: p.name)
def test_no_partial_path_selectors_for_nested_fold_subworkflows(config):
    # A 'SUBWF:PROC' selector only matches when SUBWF is the entry-level workflow;
    # fold's subworkflows sit under FOLD:FOLD_CORE:... or FOLD_PULLDOWN:...
    nested = {"BOLTZ_FOLD", "ALPHAFOLD2", "ROSETTAFOLD3_FOLD", "PROTENIX_FOLD",
              "ALPHAFOLD3_FOLD", "OPENFOLD3_FOLD", "FOLD_PREDICT", "FOLD_MSA"}
    bad = [s for s in _selectors(config) if ":" in s and s.split(":")[0] in nested]
    assert not bad, f"{config.name}: partial-path selectors match nothing: {bad}"
