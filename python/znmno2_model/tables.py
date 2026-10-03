"""The spline tables (data/*.txt), shared with the Fortran and C++ programs. Values are float32 (as in the
original programs), stored with 9 significant digits, which round-trips them exactly."""
from __future__ import annotations

from importlib import resources

import numpy as np


def _numbers(name):
    text = (resources.files("znmno2_model") / "data" / name).read_text()
    vals = " ".join(line for line in text.splitlines() if not line.lstrip().startswith("#")).split()
    return vals


def ph_spline():
    """(zn_knots (153,), mn_knots (153,), coeffs (149, 149)) as float32."""
    v = _numbers("ph_spline_phreeqc.txt")
    nz, nm, cz, cm = (int(x) for x in v[:4])
    a = np.array(v[4:], dtype=np.float32)
    zk, mk = a[:nz], a[nz:nz + nm]
    c = a[nz + nm:].reshape(cz, cm)
    return zk, mk, c


def r3_ocp():
    """(theta (51,), coeffs (50, 4)) as float32."""
    v = _numbers("r3_ocp_spline.txt")
    n = int(v[0])
    a = np.array(v[1:], dtype=np.float32)
    return a[:n], a[n:].reshape(n - 1, 4)
