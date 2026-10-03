"""The original programs' pH surrogate, pH(Zn, Mn): a bicubic B-spline fitted to PHREEQC runs.

Evaluated here in double precision for the corrected model's ph_mode = 'spline' (the faithful ports keep
the original single-precision evaluation). Same knots, coefficients and argument convention as the
original (its swapped arguments undo a transposed coefficient read). Inputs outside the knot range are
clamped to it. Needs the spline data file (data/ph_spline_phreeqc.txt).
"""
from __future__ import annotations


import numpy as np

K = 3


def _load():
    from .tables import ph_spline
    return tuple(np.asarray(a, float) for a in ph_spline())


_DATA = None


def _basis(t, i, x):
    """The K+1 nonzero cubic B-spline basis values at x in knot span i (t[i] <= x < t[i+1]), Cox-de Boor."""
    N = np.zeros(K + 1)
    N[0] = 1.0
    left, right = np.zeros(K + 1), np.zeros(K + 1)
    for j in range(1, K + 1):
        left[j] = x - t[i + 1 - j]
        right[j] = t[i + j] - x
        saved = 0.0
        for r in range(j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        N[j] = saved
    return N


def _span(t, x):
    n = len(t) - K - 1
    i = int(np.searchsorted(t, x, side="right")) - 1
    return min(max(i, K), n - 1)


def ph(zn_molar, mn_molar):
    """pH from total Zn and Mn [mol/L] (arrays)."""
    global _DATA
    if _DATA is None:
        _DATA = _load()
    tx, ty, C = _DATA
    zn = np.atleast_1d(np.asarray(zn_molar, float))
    mn = np.atleast_1d(np.asarray(mn_molar, float))
    out = np.empty(zn.shape)
    for q in range(zn.size):
        # the original evaluates its first knot axis at log10(Mn) and the second at log10(Zn)
        a = float(np.clip(np.log10(mn.flat[q]), tx[K], tx[-K - 1] - 1e-12))
        b = float(np.clip(np.log10(zn.flat[q]), ty[K], ty[-K - 1] - 1e-12))
        i, j = _span(tx, a), _span(ty, b)
        Bx, By = _basis(tx, i, a), _basis(ty, j, b)
        out.flat[q] = Bx @ C[i - K:i + 1, j - K:j + 1] @ By
    return out
