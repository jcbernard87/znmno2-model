"""zmc_transfer (#5): a nearly full insertion host turns into ZMO, the original programs' transfer (deviation M-9)
made continuous. The host's Mn is an unknown; the transferred Mn takes z_ZMO Zn along and the host keeps the rest."""
import numpy as np
import pytest

from znmno2_model.model import NH, TH, ZM, Model, sigmoid
from znmno2_model.params import Params
from znmno2_model.simulate import run

FAST = dict(n_probe=4, n_sep=6, n_cath=8)
DEEP = "cc I=400 Vmin=0.2; rest t=3600"           # R1 and R2 off: the host fills, then sits full


def deep(on, steps=DEEP, **kw):
    p = Params(steps=steps, end_on_cutoff=False, R1_on=False, R2_on=False, zhs="off", zmc_transfer=on, **FAST, **kw)
    return p, Model(p), run(p)


def test_off_keeps_the_host():
    p, m, r = deep(False)
    c = m.cath
    assert np.all(r.final_state[c, NH] == m.n_host)


def test_full_host_turns_into_zmo_conserving_everything():
    """The host's Mn goes to ZMO (Mn conserved to round-off); S and H are conserved, and Zn gains exactly the
    anode's charge / 2F, as without the transfer."""
    p, m, r = deep(True)
    x0, x = m.initial_state(), r.final_state
    c = m.cath
    assert np.all(x[c, NH] < 0.1 * m.n_host)                     # most of the host has turned
    dzm = float(np.sum((x[c, ZM] - x0[c, ZM]) * (m.mesh.area * m.mesh.dx)[c]))
    dnh = float(np.sum((x[c, NH] - x0[c, NH]) * (m.mesh.area * m.mesh.dx)[c]))
    assert dzm == pytest.approx(-dnh, rel=1e-10)
    inv0, inv = m.inventory(x0), m.inventory(x)
    for k in ("Mn", "S"):
        assert abs(inv[k] - inv0[k]) < 1e-12 * inv0[k]
    assert abs(inv["H"] - inv0["H"]) < 1e-12 * inv0["S"]
    q = r.column("mAhg")[-1] * 1e-3 * 3600.0 * p.mass                # charge passed [C]
    assert inv["Zn"] - inv0["Zn"] == pytest.approx(q / (2.0 * p.F), rel=1e-6)


def test_the_host_keeps_the_zinc_the_transfer_leaves():
    """Below z_ZMO the transferred Mn takes more Zn than the host holds per Mn, so theta falls back towards the
    onset theta_t, where the transfer switches off."""
    p, m, r = deep(True)
    th = sigmoid(r.final_state[m.cath, TH])
    assert np.all(th < m.theta_t + 3 * 0.005) and np.all(th > m.theta_t - 0.05)


def test_discharge_and_charge_with_and_without_the_transfer():
    """Deep discharge, rest, charge: with the transfer, part of the host's Mn ends in ZMO, so less capacity is
    left for insertion on the charge (and in the next discharge)."""
    steps = "cc I=400 Vmin=0.2; rest t=3600; cc I=-400 Vmax=1.9"
    _, m0, r0 = deep(False, steps)
    _, m1, r1 = deep(True, steps)
    assert r0.exit_reason == r1.exit_reason == "end_of_protocol"
    q0, q1 = (np.max(r.column("mAhg")) - r.column("mAhg")[-1] for r in (r0, r1))
    assert q1 < q0                                             # less comes back on the charge
    assert np.all(r1.final_state[m1.cath, ZM] > r0.final_state[m0.cath, ZM])
