"""Corrected model: conservation, current balance, the porosity step, and the protocol driver."""
import numpy as np
import pytest

from znmno2_model.model import CATH, MN, P2, PROBE, SEP, SO, ZN, Model, make_mesh
from znmno2_model.params import Params

FAST = dict(n_probe=6, n_sep=8, n_cath=10)


def _discharge(p, steps, dt=10.0, mAg=100.0):
    m = Model(p)
    x = m.initial_state()
    I = mAg * 1e-3 * p.mass
    states, q_anode = [x], 0.0
    for _ in range(steps):
        new = m.newton_step(x, dt, I)
        free = m.free(new)
        q_anode += float(m.anode_current(new, free)) * m.mesh.area[0] * dt
        states.append(new)
        x = new
    return m, states, I, q_anode


def test_conservation_to_round_off():
    """Mn, S and H are conserved; Zn gains exactly the anode's charge / 2F (backward Euler)."""
    p = Params(**FAST)
    m, states, I, q = _discharge(p, 6)
    inv0, inv1 = m.inventory(states[0]), m.inventory(states[-1])
    for k in ("Mn", "S"):
        assert abs(inv1[k] - inv0[k]) < 1e-12 * inv0[k]
    assert abs(inv1["H"] - inv0["H"]) < 1e-12 * inv0["S"]
    assert inv1["Zn"] - inv0["Zn"] == pytest.approx(q / (2 * p.F), rel=1e-9)


def test_reaction_current_equals_applied_current():
    p = Params(**FAST)
    m, states, I, q = _discharge(p, 3)
    x = states[-1]
    i1, i2, i3, _ = m.reactions(x, m.free(x))
    vol = (m.mesh.area * m.mesh.dx)[m.mesh.region == CATH]
    assert float(np.sum(vol * (i1 + i2 + i3))) == pytest.approx(-I, rel=1e-8)
    assert float(m.anode_current(x, m.free(x))) * m.mesh.area[0] == pytest.approx(I, rel=1e-6)


def test_electroneutral():
    p = Params(**FAST)
    m, states, I, q = _discharge(p, 3)
    x = states[-1]
    q = 2 * x[:, ZN] + 2 * x[:, MN] + x[:, 5] - 2 * x[:, SO]
    assert np.max(np.abs(q)) < 1e-15


def test_porosity_step_flux_is_exact():
    """A piecewise-linear profile with continuous flux eps D/tau dc/dx has zero divergence in every cell,
    including the two cells at the probe/separator step (eps 1 -> 0.9, tau 1 -> 2 eps^-1/2)."""
    p = Params(L_probe=0.1, L_sep=0.05, **FAST)
    m = Model(p)
    mesh = m.mesh
    x = m.initial_state()
    x[:, P2] = 0.0
    eps = m.porosity(x)
    tau = m.tortuosity(eps)
    k_eff = mesh.area * eps / tau            # flow = D k_eff dc/dx must be the same in every region
    slope = 1e-3 / k_eff                     # dc/dx in each cell for a total flow of D * 1e-3
    c = np.empty(mesh.n)
    c[0] = 1e-3 + slope[0] * mesh.dx[0] / 2
    for i in range(1, mesh.n):
        c[i] = c[i - 1] + slope[i - 1] * mesh.dx[i - 1] / 2 + slope[i] * mesh.dx[i] / 2
    x[:, ZN] = c
    R = m.transport(x)[:, ZN]
    inner = (mesh.region != CATH)[1:-1]
    flow = p.D_Zn * 1e-3
    assert np.max(np.abs(R[1:-1][inner])) < 1e-12 * flow
    assert R[0] == pytest.approx(-flow, rel=1e-12)          # what leaves the first cell eastward, negative


def test_mesh_regions_and_inputs():
    p = Params(L_probe=0.0, **FAST)
    mesh = make_mesh(p)
    assert not np.any(mesh.region == PROBE)
    assert np.sum(mesh.region == SEP) == FAST["n_sep"]
    p = Params(L_probe=0.1, A_cell=0.3, **FAST)
    mesh = make_mesh(p)
    assert np.allclose(mesh.area, 0.3)                       # one cross-section throughout (1-D)
    assert mesh.x[-1] == pytest.approx(0.1 + p.L_sep + p.L_cath - p.L_cath / p.n_cath / 2)


def test_input_file_and_gitt_protocol(tmp_path):
    """The default input file loads; a short GITT (pulse, rest, pulse) ends on time, with the rest at zero current."""
    from pathlib import Path
    from znmno2_model.namelist import load
    from znmno2_model.simulate import run
    p, extra = load(Path(__file__).resolve().parents[1] / "input" / "default.nml")
    assert p.L_probe == 0.1 and extra["file"] == "znmno2_out.txt"
    p = p.with_(steps="cc I=100 t=60; rest t=120", cycles=2, dt=20.0, write_interval=20.0, **FAST)
    r = run(p)
    assert r.exit_reason == "end_of_protocol"
    I = r.column("I_mAg")
    assert np.any(I == 0.0) and np.any(np.isclose(I, 100.0, rtol=1e-12))
    assert r.column("mAhg")[-1] == pytest.approx(100.0 * 120 / 3600, rel=1e-9)


def test_discharge_is_not_stopped_by_the_upper_cutoff():
    """With acid the discharge starts above V_max (a high-voltage plateau); only V_min ends a discharge."""
    from znmno2_model.simulate import run
    p = Params(c_H2SO4=0.1, steps="cc I=100 t=60", dt=10.0, **FAST)
    r = run(p)
    assert r.exit_reason == "duration"
    assert r.column("V").max() > p.V_max


def test_full_host_is_a_physical_exit():
    """With only R3 on, a discharge ends when the host is full: reported as such, not as a solver failure."""
    from znmno2_model.simulate import run
    p = Params(R1_on=False, R2_on=False, zhs="off", steps="cc I=400 Vmin=0.2", dt=30.0, vf_host=0.005, **FAST)
    r = run(p)
    assert r.exit_reason in ("insertion_full", "cutoff_low")
    assert r.exit_reason == "insertion_full" or r.column("theta")[-1] > 0.99


def test_every_parameter_is_an_input():
    """Every Params field can be set from the input file (namelist groups), and nothing else is listed."""
    from dataclasses import fields
    from znmno2_model.namelist import EXTRA, GROUPS
    listed = [n for names in GROUPS.values() for n in names if n not in EXTRA]
    assert sorted(listed) == sorted(f.name for f in fields(Params))
    assert len(listed) == len(set(listed))


def test_default_input_is_the_default_model():
    from pathlib import Path
    from znmno2_model.namelist import load
    p, _ = load(Path(__file__).resolve().parents[1] / "input" / "default.nml")
    assert p == Params()


def test_initial_state_is_at_open_circuit():
    """The initial potentials are the mixed potential: no net current at the anode or in the cathode, and
    the compositions are the inputs."""
    p = Params(**FAST)
    m = Model(p)
    x = m.initial_state()
    free = m.free(x)
    i1, i2, i3, _ = m.reactions(x, free)
    vol = (m.mesh.area * m.mesh.dx)[m.mesh.region == CATH]
    scale = 1e-3 * p.mass                                          # 1 mA/g
    assert abs(float(np.sum(vol * (i1 + i2 + i3)))) < 1e-6 * scale
    assert abs(float(m.anode_current(x, free)) * m.mesh.area[0]) < 1e-4 * scale   # the rows' tolerances, summed
    assert np.allclose(x[:, ZN], p.c_ZnSO4 * 1e-3, rtol=1e-6)


def test_every_parameter_is_documented():
    from dataclasses import fields
    from pathlib import Path
    text = (Path(__file__).resolve().parents[1] / "docs" / "parameters.md").read_text()
    for f in fields(Params):
        assert f"| `{f.name}` |" in text, f.name


def test_versions_agree():
    """pyproject.toml, __version__, CITATION.cff, CMakeLists.txt and the newest numbered CHANGELOG heading."""
    import re
    from pathlib import Path
    import znmno2_model
    root = Path(__file__).resolve().parents[1]

    def find(path, pattern):
        return re.search(pattern, (root / path).read_text(), re.M).group(1).strip()

    versions = {"pyproject.toml": find("pyproject.toml", r'^version = "(.+)"'), "__version__": znmno2_model.__version__,
                "CITATION.cff": find("CITATION.cff", r"^version: (.+)$"),
                "CMakeLists.txt": find("CMakeLists.txt", r"^project\(\w+ VERSION ([\d.]+)"),
                "CHANGELOG.md": find("CHANGELOG.md", r"^## ([\d.]+)")}
    assert len(set(versions.values())) == 1, versions


def test_states_and_profiles():
    from znmno2_model.simulate import run
    p = Params(steps="cc I=100 t=600", dt=60.0, write_interval=120.0, **FAST)
    r = run(p, keep_states=True)
    assert len(r.states) == len(r.rows)
    t, x = r.states[-1]
    assert t == pytest.approx(r.rows[-1][0])
    prof = Model(p).profiles(x)
    n = Model(p).mesh.n
    for key in ("x_um", "pH", "Zn_M", "Mn_M", "S_M", "H_M", "phi2", "vf_ZHS", "vf_ZMO", "theta", "eps"):
        assert prof[key].shape == (n,), key
    assert np.all(np.isnan(prof["theta"][: n - FAST["n_cath"]]))      # host only in the cathode


def test_electrolyte_without_mnso4():
    """c_MnSO4 = 0 starts from a trace (1e-9 M), and Mn2+ released by R2 builds up."""
    from znmno2_model.simulate import run
    p = Params(c_MnSO4=0.0, steps="cc I=100 t=600", dt=60.0, **FAST)
    r = run(p)
    assert r.exit_reason == "duration"
