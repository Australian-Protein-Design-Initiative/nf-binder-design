"""MkDocs hook that re-renders the metro-map diagrams before each build.

The ``.mmd`` files beside this one are the source of truth. The rendered SVG and
PNG stay committed under ``docs/images/`` so GitHub can show them when the
workflow pages are read in the repository rather than on the docs site, and this
hook keeps them in step with the sources on every ``mkdocs build``/``serve``.

Rendering is skipped when the outputs are already newer than their source, so a
normal build pays nothing. nf-metro is invoked as ``python -m nf_metro`` rather
than through the console script so it is found wherever mkdocs itself is
installed, without depending on PATH.
"""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

log = logging.getLogger("mkdocs.hooks.metro")

_HERE = Path(__file__).resolve().parent
_IMAGES = _HERE.parent / "docs" / "images"

# Matches the width the maps were tuned at; see README.md in this directory.
_RASTER_WIDTH = "2265"


def _run(args: list[str]) -> bool:
    proc = subprocess.run(
        [sys.executable, "-m", "nf_metro", *args],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        log.error("nf-metro %s failed:\n%s", " ".join(args), proc.stderr.strip())
        return False
    for line in proc.stderr.splitlines():
        if line.strip().startswith("- "):
            log.info("nf-metro: %s", line.strip()[2:])
    return True


def _render(source: Path, svg: Path, png: Path) -> None:
    if not _run(["render", str(source), "--mode", "light", "-o", str(svg)]):
        return
    # drawsvg omits the final newline; keep the committed file POSIX-clean.
    text = svg.read_text()
    if not text.endswith("\n"):
        svg.write_text(text + "\n")
    _run(
        [
            "render",
            str(source),
            "--mode",
            "light",
            "--raster-width",
            _RASTER_WIDTH,
            "-o",
            str(png),
        ]
    )


def on_pre_build(config) -> None:  # noqa: ARG001 - mkdocs passes the config
    try:
        import nf_metro  # noqa: F401
    except ImportError:
        log.warning(
            "nf-metro is not installed, so the metro-map diagrams were not "
            "re-rendered; the committed images under docs/images/ are used as "
            "they are. Install it with: pip install -r docs/requirements.txt"
        )
        return

    for source in sorted(_HERE.glob("*_metro_map.mmd")):
        svg = _IMAGES / f"{source.stem}.svg"
        png = _IMAGES / f"{source.stem}.png"
        fresh = (
            svg.exists()
            and png.exists()
            and svg.stat().st_mtime >= source.stat().st_mtime
            and png.stat().st_mtime >= source.stat().st_mtime
        )
        if fresh:
            log.debug("%s is up to date", svg.name)
            continue
        log.info("rendering %s", source.name)
        _render(source, svg, png)
