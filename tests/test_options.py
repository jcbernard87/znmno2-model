"""The corrected model's options (plan: option matrix B): conservation, limits, and the Jacobian."""
import numpy as np
import pytest

from znmno2_model.model import CATH, H, MN, N, P2, SO, ZH, ZN, ZO, Model
from znmno2_model.params import Params

FAST = dict(n_probe=4, n_sep=6, n_cath=8)

CONFIGS = {
    "default": {},
    "no_probe": dict(L_probe=0.0),
    "ph_zhs_equilibrium": dict(ph_mode="zhs_equilibrium"),
    "ph_fixed": dict(ph_mode="fixed"),
    "ph_spline": dict(ph_mode="spline"),
    "no_H": dict(species="no_H"),
    "quasi": dict(transport="quasi"),
    "totals": dict(basis="totals"),
    "zhs_equilibrium": dict(zhs="equilibrium"),
    "zhs_lumped": dict(zhs="lumped"),
    "zhs_off": dict(zhs="off"),
    "nucleation": dict(zhs_nucleation=1.05),
    "zno": dict(ZnO_on=True, logK_ZnO=8.0),                 # low log K: supersaturated at pH ~4.3
    "znoh2": dict(ZnOH2_on=True, logK_ZnOH2=8.0),
    "r3_nernst": dict(r3_ocp="nernst"),
    "r3_spline": dict(r3_ocp="spline"),
    "r1_r2_off": dict(R1_on=False, R2_on=False),
    "acid": dict(c_H2SO4=0.1),
}


def _run(p, steps=5, dt=10.0, mAg=100.0):
    """Steps of dt through the driver's advance (sub-steps halved on Newton failure, as in a real run);
    q is the charge that passed the anode."""
    from znmno2_model.driver import advance

    class Counting:                                  # integrates the anode current over the sub-steps
        def __init__(self, m):
            self.m, self.q = m, 0.0

        def newton_step(self, x, h, I):
            new = self.m.newton_step(x, h, I)
            self.q += float(self.m.anode_current(new, self.m.free(new))) * self.m.mesh.area[0] * h
            return new

    m = Model(p)
    x = m.initial_state()
    I = mAg * 1e-3 * p.mass
    st = Counting(m)
    xs = [x]
    for _ in range(steps):
        x, h_done, _ = advance(st, x, dt, I)
        assert h_done == pytest.approx(dt)
        xs.append(x)
    return m, xs, I, st.q


@pytest.mark.parametrize("name", sorted(CONFIGS))
def test_option_runs_and_conserves(name):
    p = Params(**FAST, **CONFIGS[name])
    m, xs, I, q = _run(p)
    a, b = m.inventory(xs[0]), m.inventory(xs[-1])
    assert abs(b["Mn"] - a["Mn"]) < 1e-12 * a["Mn"]
    assert abs(b["S"] - a["S"]) < 1e-12 * a["S"]
    if p.species == "with_H":
        assert abs(b["H"] - a["H"]) < 1e-12 * a["S"]
    assert b["Zn"] - a["Zn"] == pytest.approx(q / (2 * p.F), rel=1e-8)
    x = xs[-1]
    if name in ("zno", "znoh2"):                                   # the case must exercise the precipitate
        assert x[m.cath, ZO if name == "zno" else 11].min() > 1e-9
    i1, i2, i3, _ = m.reactions(x, m.free(x))
    vol = (m.mesh.area * m.mesh.dx)[m.mesh.region == CATH]
    assert float(np.sum(vol * (i1 + i2 + i3))) == pytest.approx(-I, rel=1e-7)


def _fd_blocks(m, x, old, dt, I):
    """Central differences of the whole residual, by column colouring (3 colours)."""
    n = m.mesh.n
    A = np.zeros((n, N, N)); B = np.zeros((n, N, N)); D = np.zeros((n, N, N))
    h = 1e-6 * np.maximum(np.abs(x), np.array([1, 1, 1e-6, 1e-6, 1e-6, 1e-8, 1e-6, 1e-6, 1, 1e-6, 1e-6, 1e-6]))
    for k in range(N):
        for col in range(3):
            idx = np.arange(col, n, 3)
            xp, xm = x.copy(), x.copy()
            xp[idx, k] += h[idx, k]
            xm[idx, k] -= h[idx, k]
            d = m.residual(xp, old, dt, I) - m.residual(xm, old, dt, I)
            for j in idx:
                B[j, :, k] = d[j] / (2 * h[j, k])
                if j > 0:
                    D[j - 1, :, k] = d[j - 1] / (2 * h[j, k])
                if j < n - 1:
                    A[j + 1, :, k] = d[j + 1] / (2 * h[j, k])
    return A, B, D


@pytest.mark.parametrize("name", ["default", "quasi", "no_H", "zhs_equilibrium", "zno", "totals", "acid"])
def test_jacobian_matches_central_differences(name):
    p = Params(**FAST, **CONFIGS[name])
    m, xs, I, q = _run(p, steps=2)
    old, x = xs[-2], xs[-1].copy()
    rng = np.random.default_rng(0)
    x[:, P2] += 1e-3 * rng.standard_normal(m.mesh.n)               # off the solution, gradients everywhere
    x[:, ZN] *= 1 + 1e-2 * rng.random(m.mesh.n)
    c = m.mesh.region == CATH
    for k in (ZH, ZO, 11):                                         # away from the kink of max(n, 0) at n = 0
        x[c, k] = np.maximum(x[c, k], 1e-6)
    m.free(x)
    R, A, B, D = m.blocks(x, old, 10.0, I)
    Af, Bf, Df = _fd_blocks(m, x, old, 10.0, I)
    for a, f in ((A, Af), (B, Bf), (D, Df)):
        scale = np.max([np.max(np.abs(z) * m.typ, axis=2) for z in (Af, Bf, Df)], axis=0)[..., None] + 1e-300
        err = np.abs(a - f) * m.typ / scale
        assert err.max() < 1e-4, np.unravel_index(np.argmax(err), err.shape)


def test_fixed_ph_is_used():
    p = Params(**FAST, ph_mode="fixed", pH_fixed=4.8)
    m = Model(p)
    x = m.initial_state()
    cb = m._basis(x, m.free(x))
    assert np.allclose(-np.log10(cb[:, 0]), 4.8)


def test_zhs_equilibrium_ph_is_ksp_closed_form():
    p = Params(**FAST, ph_mode="zhs_equilibrium")
    m = Model(p)
    x = m.initial_state()
    free = m.free(x)
    cb = m._basis(x, free)
    w = np.exp((4 * np.log(free[:, 1]) + np.log(free[:, 3]) - p.logK_ZHS * np.log(10)) / 6) / cb[:, 0]
    assert np.allclose(w, 1.0)


def test_spline_ph_matches_original_surrogate():
    from znmno2_model import phspline
    p = Params(**FAST, ph_mode="spline")
    m = Model(p)
    x = m.initial_state()
    cb = m._basis(x, m.free(x))
    assert np.allclose(-np.log10(cb[:, 0]), phspline.ph(x[:, ZN] * 1e3, x[:, MN] * 1e3))


def test_no_H_equals_with_H_when_nothing_consumes_protons():
    """With R1, R2 and ZHS off and no acid, H_T stays 0: both species options give the same cell."""
    base = dict(R1_on=False, R2_on=False, zhs="off", **FAST)
    va = [Model(Params(**base)).voltage(x, 0) for x in _run(Params(**base))[1]]
    vb = [Model(Params(species="no_H", **base)).voltage(x, 0) for x in _run(Params(species="no_H", **base))[1]]
    assert np.allclose(va, vb, rtol=0, atol=1e-10)


def test_quasi_equals_ions_without_complexes(tmp_path):
    """With every complex suppressed and OH- negligible (acid), each total is one ion: the quasi-particle
    flux equals the ion flux."""
    from znmno2_model.eqchem import Equilibria
    names = [n for n in Equilibria().names[4:]]
    f = tmp_path / "k.txt"
    f.write_text("\n".join(f"{n} -60" for n in names if n != "OH-"))
    kw = dict(logk_file=str(f), c_H2SO4=0.05, D_OH=9.0e-5, **FAST)
    mi, mq = Model(Params(**kw)), Model(Params(transport="quasi", **kw))
    x = mi.initial_state()
    rng = np.random.default_rng(1)
    x[:, P2] += 1e-2 * rng.standard_normal(mi.mesh.n)
    for k in (ZN, MN, H):
        x[:, k] *= 1 + 0.1 * rng.random(mi.mesh.n)
    x[:, SO] = x[:, ZN] + x[:, MN] + x[:, H] / 2                     # electroneutral, H_T > 0
    Ri, Rq = mi.transport(x), mq.transport(x)
    for k in (ZN, MN, SO, H):
        assert np.allclose(Rq[:, k], Ri[:, k], rtol=1e-9, atol=1e-12 * np.max(np.abs(Ri[:, k])))


def test_r3_nernst_ocp():
    p = Params(**FAST, r3_ocp="nernst")
    m = Model(p)
    assert m.u3(np.array([0.0]), np.array([1.0]))[0] == pytest.approx(p.U3_nernst)
    assert m.u3(np.array([np.log(3.0)]), np.array([1.0]))[0] == pytest.approx(p.U3_nernst - np.log(3.0) / (2 * m.f))


def test_precipitate_forms_only_when_supersaturated():
    p = Params(**FAST, ZnO_on=True)
    m = Model(p)
    x = m.initial_state()
    free = m.free(x)
    assert np.all(m.precipitation(x, free, ZO) < 1e-30)            # pH 4.3: ZnO undersaturated, none present
    free2 = free.copy()
    free2[:, 0] = 10 ** -7.0                                       # pH 7: supersaturated
    assert np.all(m.precipitation(x, free2, ZO) > 0)


def test_equilibrium_zhs_is_complementary():
    """zhs = equilibrium: after steps that drive the pH to saturation, w <= 1 everywhere, and where ZHS
    exists it is at saturation."""
    p = Params(**FAST, zhs="equilibrium", c_ZnSO4=2.0, theta0=0.5, vf_host=0.0001)
    m, xs, I, q = _run(p, steps=30, dt=30.0, mAg=200.0)
    x = xs[-1]
    w = m.zhs_w(m.free(x))
    c = m.mesh.region == CATH
    assert np.all(w <= 1 + 1e-6)
    present = x[c, ZH] > 1e-9
    assert present.any()
    assert np.allclose(w[present], 1.0, atol=1e-6)


def test_invalid_option_rejected():
    with pytest.raises(ValueError):
        Params(ph_mode="phreeqc")
