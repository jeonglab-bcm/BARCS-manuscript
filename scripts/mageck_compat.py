#!/usr/bin/env python3
"""Run official MAGeCK 0.5.9.5 with current NumPy/Python CNV compatibility.

MAGeCK's CNV reader expects ``numpy.genfromtxt`` to return byte strings and
uses the removed ``numpy.float`` alias.  This wrapper restores those two
legacy runtime behaviours without changing MAGeCK's statistical code.
"""

from __future__ import annotations

import runpy
from pathlib import Path

import numpy as np


_genfromtxt = np.genfromtxt


def _legacy_genfromtxt(*args, **kwargs):
    kwargs.setdefault("encoding", "bytes")
    return _genfromtxt(*args, **kwargs)


np.genfromtxt = _legacy_genfromtxt
np.float = float

mageck_entrypoint = (
    Path(__file__).resolve().parents[1] / ".venv" / "bin" / "mageck"
)
runpy.run_path(str(mageck_entrypoint), run_name="__main__")
