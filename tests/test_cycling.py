"""Cycling: constant-voltage steps, conservation over several cycles, and relaxation to equilibrium at rest."""
import numpy as np
import pytest

from znmno2_model.model import CATH, MN, SO, ZN, Model
from znmno2_model.params import Params
from znmno2_model.simulate import Stepper, run
from znmno2_model.driver import run_protocol
from znmno2_model.simulate import Result

FAST = dict(n_probe=4, n_sep=6, n_cath=8)


def test_cc_cv_charge():
    p = Params(steps="cc I=200 Vmin=1.2; cc I=-200 Vmax=1.75; cv V=1.75 Imin=20", end_on_cutoff=False,
               dt=30.0, write_interval=30.0, **FAST)
    r = run(p)
    assert r.exit_reason == "end_of_protocol"
    step, V, I = r.column("step"), r.column("V"), r.column("I_mAg")
    cv = step == 3
    assert cv.sum() >= 2
    assert np.allclose(V[cv], 1.75, atol=1e-6)
    assert np.all(I[cv] < 0) and abs(I[cv][-1]) <= 20.0 + 1e-6
    assert np.all(np.diff(np.abs(I[cv])) <= 1e-9)                    # the charging current decays


def test_conservation_over_cycles():
    """Three discharge/charge cycles: Mn, S, H conserved; Zn changes by exactly the charge passed / 2F."""
    p = Params(steps="cc I=200 Vmin=1.2; cc I=-200 Vmax=1.8", cycles=3, end_on_cutoff=False, dt=30.0, **FAST)
    st = Stepper(p)
    r = run_protocol(st, result=Result())
    assert r.exit_reason == "end_of_protocol"
    m = st.model
    a, b = m.inventory(st.x0), m.inventory(r.final_state)
    for k in ("Mn", "S"):
        assert abs(b[k] - a[k]) < 1e-11 * a[k]
    assert abs(b["H"] - a["H"]) < 1e-11 * a["S"]
    q = r.rows[-1][3] * p.mass * 3.6                                  # net charge passed [C] (mAh/g -> C)
    assert b["Zn"] - a["Zn"] == pytest.approx(q / (2 * p.F), rel=1e-5, abs=1e-6 * a["Zn"])


def test_rest_relaxes_to_equilibrium():
    """After a partial discharge, a long rest: the electrolyte becomes uniform and the net reaction rate in
    every cell vanishes (R1, irreversible, is switched off so that equilibrium exists)."""
    p = Params(steps="cc I=200 t=1800; rest t=345600", dt=600.0, R1_on=False, **FAST)
    st = Stepper(p)
    r = run_protocol(st, result=Result())
    x = r.final_state
    m = st.model
    for k in (ZN, MN, SO):
        assert np.ptp(x[:, k]) < 1e-6 * np.mean(x[:, k])
    i1, i2, i3, rz = m.reactions(x, m.free(x))
    scale = 1e-3 * p.mass / float(np.sum((m.mesh.area * m.mesh.dx)[m.mesh.region == CATH]))   # 1 mA/g per cm3
    assert np.max(np.abs(i1 + i2 + i3)) < 1e-6 * scale
    assert np.max(np.abs(i2)) < 1e-3 * scale and np.max(np.abs(i3)) < 1e-3 * scale


def test_physical_limit_ends_the_step_not_the_protocol():
    """A discharge that ends at a physical limit (host full, R1 and R2 off) is followed by the charge step."""
    p = Params(steps="cc I=400 Vmin=0.2; cc I=-400 Vmax=1.8", end_on_cutoff=False, R1_on=False, R2_on=False,
               zhs="off", dt=30.0, vf_host=0.005, **FAST)
    r = run(p)
    assert 2 in set(r.column("step").astype(int))
    assert r.exit_reason == "end_of_protocol"


def test_spline_nernst_ocp_keeps_the_spline_and_diverges_at_the_ends():
    """Default R3 OCP: the empirical spline in the middle (within 1 mV for 0.1 <= theta <= 0.9), with
    Nernstian ends: U3 -> +inf as theta -> 0 and -inf as theta -> 1, like -(RT/2F) ln(theta/(1-theta))."""
    from znmno2_model.model import logit
    p = Params(**FAST)
    assert p.r3_ocp == "spline_nernst"
    m = Model(p)
    ms = Model(Params(**{**FAST, "r3_ocp": "spline"}))
    th = np.linspace(0.1, 0.9, 81)
    c = np.full(th.shape, 0.3)
    assert np.max(np.abs(m.u3(logit(th), c) - ms.u3(logit(th), c))) < 1e-3
    s = np.array([-40.0, -60.0, 40.0, 60.0])
    u = m.u3(s, np.full(4, 0.3))
    assert u[1] - u[0] == pytest.approx(20.0 / (2 * m.f), rel=1e-3)      # slope RT/2F per unit of s
    assert u[2] - u[3] == pytest.approx(20.0 / (2 * m.f), rel=1e-3)


def test_charge_recovers_the_discharge_with_the_default_ocp():
    """Discharge to the cutoff, then charge: nearly all the capacity comes back (insertion and
    dissolution both reverse)."""
    p = Params(steps="cc I=100 Vmin=1.0; cc I=-100 Vmax=1.85", end_on_cutoff=False, dt=30.0,
               n_probe=6, n_sep=8, n_cath=10)
    r = run(p)
    q = r.column("mAhg")
    assert q.max() > 150 and q[-1] < 0.05 * q.max()


def test_cv_hold_proceeds_in_sub_steps():
    """A 1.9 V hold right after a discharge cannot be held for a whole 600 s step at first; it proceeds in
    sub-steps until the current falls to Imin."""
    p = Params(steps="cc I=200 Vmin=1.2; cv V=1.9 Imin=5", end_on_cutoff=False, dt=600.0, **FAST)
    r = run(p)
    assert r.exit_reason == "end_of_protocol"
    assert abs(r.column("I_mAg")[-1]) <= 5.0 and r.column("V")[-1] == pytest.approx(1.9, abs=1e-6)
