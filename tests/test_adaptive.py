"""Adaptive time stepping (dt_mode = 'adaptive', #7): backward Euler with the local error estimated from a
linear predictor and controlled against rtol/atol (scaled per variable); steps land on the output times. The
Fortran and C++ programs take the same steps (tests/test_ports.py, the adaptive cases)."""
import numpy as np
import pytest

from znmno2_model.params import Params
from znmno2_model.simulate import run

FAST = dict(n_probe=4, n_sep=6, n_cath=8)
S = "cc I=200 t=1800"


def voltages(r):
    return {round(t * 3600, 6): v for t, v in zip(r.column("t_h"), r.column("V"))}


def test_fixed_is_the_default():
    assert Params().dt_mode == "fixed"
    with pytest.raises(ValueError, match="dt_mode must be one of"):
        Params(dt_mode="variable")


def test_adaptive_converges_with_rtol_and_takes_fewer_steps():
    """Against fixed 1 s steps (rows every 60 s, on the same times as adaptive rows): the voltage error falls
    with the tolerance, at a fraction of the steps."""
    ref = run(Params(steps=S, dt=1.0, **FAST))
    vr = voltages(ref)

    def error(r):
        v = voltages(r)
        common = sorted(set(v) & set(vr))
        assert len(common) == 31
        return max(abs(v[t] - vr[t]) for t in common)

    loose = run(Params(steps=S, dt_mode="adaptive", rtol=1e-2, atol=1e-2, **FAST))
    tight = run(Params(steps=S, dt_mode="adaptive", rtol=1e-4, atol=1e-4, **FAST))
    e_loose, e_tight = error(loose), error(tight)
    assert e_tight < 3e-4 and e_tight < e_loose / 3                    # [V]
    assert loose.steps < tight.steps < ref.steps / 4
    assert np.all(np.abs(np.diff(tight.column("t_h") * 3600) - 60.0) < 1e-9)
