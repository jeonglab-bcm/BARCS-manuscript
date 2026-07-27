#!/usr/bin/env python3
"""Run official MAGeCK 0.5.9.5 with current NumPy/Python CNV compatibility.

MAGeCK's CNV reader expects ``numpy.genfromtxt`` to return byte strings and
uses the removed ``numpy.float`` alias.  This wrapper restores those two
legacy runtime behaviours without changing MAGeCK's statistical code.
"""

from __future__ import annotations

import os
import runpy
import shutil
import sys
from pathlib import Path

import numpy as np


_genfromtxt = np.genfromtxt


def _legacy_genfromtxt(*args, **kwargs):
    kwargs.setdefault("encoding", "bytes")
    return _genfromtxt(*args, **kwargs)


np.genfromtxt = _legacy_genfromtxt
np.float = float

def find_mageck_entrypoint() -> Path:
    """Locate the official MAGeCK console script.

    Mirrors the resolution order in ``R/mageck.R``: an explicit ``$MAGECK``
    override, then whatever is on ``PATH`` (which is what pixi provides), then
    the legacy ad-hoc virtualenv that this repository used to assume.

    The console script is executed with :mod:`runpy` rather than spawned as a
    subprocess, because the two NumPy behaviours patched above have to be in
    place inside MAGeCK's own process.
    """
    override = os.environ.get("MAGECK")
    if override:
        return Path(override)

    on_path = shutil.which("mageck")
    if on_path:
        return Path(on_path)

    legacy = Path(__file__).resolve().parents[1] / ".venv" / "bin" / "mageck"
    if legacy.exists():
        return legacy

    sys.exit(
        "mageck_compat.py: no MAGeCK installation found. "
        "Check the toolchain with `pixi run doctor`, or set $MAGECK."
    )


runpy.run_path(str(find_mageck_entrypoint()), run_name="__main__")
