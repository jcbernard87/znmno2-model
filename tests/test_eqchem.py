"""The model's speciation (eqchem) against the general solver (speciation.py), ideal activities."""
import numpy as np
import pytest

from znmno2_model.eqchem import Equilibria
from znmno2_model.speciation import Activity, Speciation, herrmann_2023


def _grid():
    zn = np.array([2.0, 2.0, 1.0, 0.5, 2.0, 1e-3, 2.0, 0.2, 3.0])
    mn = np.array([0.05, 0.5, 0.0001, 0.1, 0.05, 1e-3, 1e-6, 0.2, 0.1])
    acid = np.array([0.0, 0.0, 0.1, 0.5, 1e-4, 0.0, 0.0, 0.0, 0.0])
    s = zn + mn + acid
    h = 2.0 * acid
    h[4] = -1e-3                 # base added (H_T < 0): ZnOH+ and hydroxo complexes carry it
    h[7] = -2e-2
    return np.stack([h, zn, mn, s], axis=-1)


def test_matches_general_solver():
    T = _grid()
    eq = Equilibria()
    x = eq.solve(T)
    ref = Speciation(herrmann_2023(), Activity(kind="ideal"))
    q = T[:, 0] + 2 * T[:, 1] + 2 * T[:, 2] - 2 * T[:, 3]          # charge of the totals
    r = ref.solve(T[:, 1], T[:, 2], s=T[:, 3], h_excess=q, charge0=0.0)
    # speciation.py's root search picks the wrong branch with acid added (validated for salts only):
    # compare the salt and base rows here; the acid rows are checked by closure in test_balances_closed
    keep = T[:, 0] <= 0.0
    T, x, r = T[keep], x[keep], {k: (v[keep] if isinstance(v, np.ndarray) else v) for k, v in r.items()}
    assert np.all(r["converged"])
    np.testing.assert_allclose(-x[:, 0], r["pH"], atol=1e-9)
    m = eq.species(x)
    for name in ("Zn+2", "ZnSO4", "Zn(SO4)2-2", "SO4-2", "HSO4-", "Mn+2", "OH-", "ZnOH+"):
        a, b = m[:, eq.names.index(name)], r["molality"][:, r["names"].index(name)]
        np.testing.assert_allclose(a, b, rtol=1e-8, err_msg=name)


def test_balances_closed():
    T = _grid()
    eq = Equilibria()
    x = eq.solve(T)
    np.testing.assert_allclose(eq.totals(x), T, rtol=1e-12, atol=1e-15)
    m = eq.species(x)
    z = np.array([1, 2, 2, -2]) @ eq.nu.T                           # species charges
    np.testing.assert_allclose(m @ z, T[:, 0] + 2 * T[:, 1] + 2 * T[:, 2] - 2 * T[:, 3], atol=1e-12)
    assert -x[2, 0] == pytest.approx(1.15, abs=0.05)                  # 1 M ZnSO4 + 0.1 M H2SO4


def test_warm_start_from_neighbour():
    T = _grid()
    eq = Equilibria()
    x = eq.solve(T)
    T2 = T * (1 + 1e-6)
    x2 = eq.solve(T2, x0=x)
    np.testing.assert_allclose(eq.totals(x2), T2, rtol=1e-12, atol=1e-15)


def test_herrmann_figure_s1():
    """2 M ZnSO4 + 0.5 M MnSO4 at pH 4: 72.7 % ZnSO4, 15.2 % Zn(SO4)2, 12.0 % Zn2+ (Herrmann et al. 2023, Fig. S1)."""
    eq = Equilibria()
    # find H_T giving pH 4 by bisection on the solver itself
    lo, hi = -1e-2, 1e-1
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        x = eq.solve(np.array([[mid, 2.0, 0.5, 2.5 + mid / 2]]))
        lo, hi = (lo, mid) if -x[0, 0] < 4.0 else (mid, hi)
    m = eq.species(x)[0]
    f = {n: m[eq.names.index(n)] / 2.0 for n in ("ZnSO4", "Zn(SO4)2-2", "Zn+2")}
    assert f["ZnSO4"] == pytest.approx(0.727, abs=2e-3)
    assert f["Zn(SO4)2-2"] == pytest.approx(0.152, abs=2e-3)
    assert f["Zn+2"] == pytest.approx(0.120, abs=2e-3)


def test_sensitivity_matches_finite_differences():
    T = _grid()
    eq = Equilibria()
    x = eq.solve(T)
    S = eq.sensitivity(x, T)
    for k in range(4):
        h = 1e-6 * np.maximum(np.abs(T[:, k]), 1e-6)
        Tp, Tm = T.copy(), T.copy()
        Tp[:, k] += h; Tm[:, k] -= h
        d = (eq.solve(Tp, x0=x) - eq.solve(Tm, x0=x)) / (2 * h[:, None])
        np.testing.assert_allclose(S[:, :, k], d, rtol=1e-5, atol=1e-6 * np.max(np.abs(d)))


def test_logk_override_and_file(tmp_path):
    f = tmp_path / "k.txt"
    f.write_text("# optimized\nZnSO4 3.0\nHSO4-  2.5\n")
    from znmno2_model.eqchem import read_logk_file
    eq = Equilibria(read_logk_file(f))
    assert eq.logk[eq.names.index("ZnSO4")] == 3.0 and eq.logk[eq.names.index("HSO4-")] == 2.5
    with pytest.raises(ValueError):
        Equilibria({"NoSuch": 1.0})


def test_cold_solve_found_in_a_run():
    """A composition from a charge (HSO4- log K raised by 0.5) on which plain Newton diverged."""
    base = Equilibria()
    eq = Equilibria({"HSO4-": base.logk[base.names.index("HSO4-")] + 0.5})
    T = np.array([[1.17629069e-03, 2.02148198e+00, 4.87407860e-02, 2.07081091e+00]])
    x = eq.solve(T)
    np.testing.assert_allclose(eq.totals(x), T, rtol=1e-12, atol=1e-15)


def test_cold_solves_over_a_wide_range():
    """Random electroneutral compositions (Zn 1e-4..3 M, Mn 1e-6..1 M, H_T -0.05..0.5 M), log K shifted
    by up to +-1: every cold solve converges and closes its balances."""
    rng = np.random.default_rng(7)
    n = 400
    zn = 10 ** rng.uniform(-4, np.log10(3), n)
    mn = 10 ** rng.uniform(-6, 0, n)
    h = rng.uniform(-0.05, 0.5, n) * np.minimum(1.0, zn)
    T = np.stack([h, zn, mn, zn + mn + h / 2], axis=1)
    for seed in range(3):
        base = Equilibria()
        shift = np.random.default_rng(seed).uniform(-1, 1, len(base.names) - 4)
        eq = Equilibria({nm: base.logk[i + 4] + s for i, (nm, s) in enumerate(zip(base.names[4:], shift))})
        x = eq.solve(T)
        np.testing.assert_allclose(eq.totals(x), T, rtol=1e-11, atol=1e-14)


def test_spline_tables():
    from znmno2_model import tables
    zk, mk, c = tables.ph_spline()
    t, oc = tables.r3_ocp()
    assert zk.shape == mk.shape == (153,) and c.shape == (149, 149) and t.shape == (51,) and oc.shape == (50, 4)
    assert all(a.dtype == np.float32 for a in (zk, mk, c, t, oc))
    assert t[0] == 0 and t[-1] == 1 and oc[0, 3] == 1.0
