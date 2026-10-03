"""The block-tridiagonal solve used by the faithful ports.

Byte-for-byte reproduction of the original programs needs their BAND/MATINV, operation by operation.
That code follows the listing in Newman's *Electrochemical Systems* (Appendix C) and is not distributed.
If a private transliteration is available (band_legacy.py in the folder named by the environment
variable ZNMNO2_ORACLE), it is used; otherwise bandsolver's legacy pivot is used, which reproduces the
original's pivoting but not its exact floating-point operation order: over full runs the output has the
same rows and exit, with values within 1e-9 relative of the original's (measured), not byte for byte. The
difference is only the compiler's fused multiply-adds: a bandsolver built with -ffp-contract=off
reproduces the original programs byte for byte.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import bandsolver
import numpy as np


def _private():
    root = os.environ.get("ZNMNO2_ORACLE")
    if not root:
        return None
    path = Path(root) / "band_legacy.py"
    if not path.is_file():
        return None
    spec = importlib.util.spec_from_file_location("znmno2_band_legacy", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_LEGACY = _private()
EXACT = _LEGACY is not None


def solve(A, B, D, G):
    if _LEGACY is not None:
        return _LEGACY.solve(A, B, D, G)
    try:
        return bandsolver.solve(A, B, D, G, pivot="legacy", singular="exact")
    except bandsolver.NonFiniteError:
        return np.full(G.shape, np.nan)       # the original propagates NaN; the run then ends on it
