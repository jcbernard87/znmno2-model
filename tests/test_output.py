"""Configurable output (#15): chosen time-series columns, and profiles against x in a long table."""
import subprocess

import numpy as np
import pytest

from test_ports import ROOT, _input, _program, _run
from znmno2_model.model import Model
from znmno2_model.params import Params
from znmno2_model.simulate import COLUMNS, PROFILE_VARS, run

FAST = dict(n_probe=4, n_sep=6, n_cath=8, dt=30.0)
STEPS = "cc I=200 t=1800; rest t=600"
OUT = dict(profile_every=7, profile_times="100 450", columns="V, pH_cath, theta")


def _table(path):
    lines = path.read_text().splitlines()
    return lines[0], np.array([[float(v) for v in ln.split()] for ln in lines[1:]])


def test_profile_variables_are_model_profiles():
    m = Model(Params(**FAST))
    assert tuple(k for k in m.profiles(m.initial_state()) if k != "x_um") == PROFILE_VARS


def test_profiles_do_not_change_the_run(tmp_path):
    """A profile takes the row's speciation: the time series with profiles at every row is the one without."""
    a = run(Params(steps=STEPS, **FAST)).array
    b = run(Params(steps=STEPS, profile_file=str(tmp_path / "p.txt"), profile_every=1, **FAST)).array
    np.testing.assert_array_equal(a, b)


def test_profile_schedule_and_values(tmp_path):
    """Profiles with every 7th row, at the first row at or after 100 s and 450 s, and with the last row; the volume
    mean of a profile over the cathode is the time series' cathode mean."""
    pp = tmp_path / "prof.txt"
    p = Params(steps=STEPS, profile_file=str(pp), profile_every=7, profile_times="450, 100", **FAST)
    r = run(p)
    t_rows = r.array[:, 0] * 3600.0
    due = {i for i in range(len(t_rows)) if i % 7 == 0}
    for tp in (100.0, 450.0):
        due.add(int(np.argmax(t_rows >= tp)))
    due.add(len(t_rows) - 1)
    assert [round(pr[0], 6) for pr in r.profiles] == [round(t_rows[i], 6) for i in sorted(due)]
    head, a = _table(pp)
    assert head.split() == ["t_h", "x_um", *PROFILE_VARS]
    m = Model(p)
    n = m.mesh.n
    assert a.shape == (len(due) * n, 2 + len(PROFILE_VARS))
    c, vol = m.cath, m.mesh.area * m.mesh.dx
    for k, i in enumerate(sorted(due)):
        prof = a[k * n:(k + 1) * n]
        for var, col in (("Zn_M", "Zn_cath_M"), ("vf_MnO2", "vf_MnO2"), ("eps", "eps_cath")):
            v = prof[:, 2 + PROFILE_VARS.index(var)]
            assert np.sum(vol[c] * v[c]) / np.sum(vol[c]) == pytest.approx(r.array[i, COLUMNS.index(col)], rel=1e-6)
        assert np.all(np.isnan(prof[~c, 2 + PROFILE_VARS.index("theta")]))


def test_default_output_is_unchanged(tmp_path):
    r = run(Params(steps="cc I=100 t=120", **FAST))
    r.write(tmp_path / "o.txt")
    assert (tmp_path / "o.txt").read_text().splitlines()[0].split() == list(COLUMNS)
    assert r.profiles == []


@pytest.mark.parametrize("lang", ["fortran", "cpp"])
def test_native_columns_and_profiles(lang, tmp_path):
    """The programs: the time series and the profiles identical in Fortran and C++ (checked through each against
    Python to the printed precision, and against each other byte for byte)."""
    exe = _program(lang)
    inp, res, pf = tmp_path / "in.nml", tmp_path / "out.txt", tmp_path / "prof.txt"
    _input(inp, res, dict(FAST, steps=STEPS, profile_file=str(pf), **OUT))
    _run(exe, inp)
    pp = tmp_path / "py_prof.txt"
    r = run(Params(steps=STEPS, profile_file=str(pp), **FAST, **OUT))
    r.write(tmp_path / "py.txt")
    assert res.read_text().splitlines()[0] == (tmp_path / "py.txt").read_text().splitlines()[0]
    assert res.read_text().splitlines()[0].split() == ["t_h", "V", "pH_cath", "theta"]
    _, an = _table(res)
    _, ap = _table(tmp_path / "py.txt")
    np.testing.assert_allclose(an, ap, rtol=1e-6)
    hn, bn = _table(pf)
    hp, bp = _table(pp)
    assert hn == hp and bn.shape == bp.shape
    np.testing.assert_allclose(bn, bp, rtol=1e-6, atol=1e-12, equal_nan=True)


def test_fortran_and_cpp_write_the_same_bytes(tmp_path):
    outs = []
    for lang in ("fortran", "cpp"):
        exe = _program(lang)
        inp, res, pf = tmp_path / f"{lang}.nml", tmp_path / f"{lang}.txt", tmp_path / f"{lang}_p.txt"
        _input(inp, res, dict(FAST, steps=STEPS, profile_file=str(pf), **OUT))
        _run(exe, inp)
        outs.append((res.read_bytes(), pf.read_bytes()))
    assert outs[0][0] == outs[1][0]
    assert outs[0][1] == outs[1][1]


@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
@pytest.mark.parametrize("kw,msg", [(dict(columns="V, Vltage"), "unknown output column 'Vltage'"),
                                    (dict(profile_vars="pH, theta_surf"), "unknown profile variable 'theta_surf'")])
def test_unknown_names_are_rejected(tmp_path, which, kw, msg):
    if "profile_vars" in kw:
        kw = dict(kw, profile_file=str(tmp_path / "p.txt"))
    if which == "python":
        with pytest.raises(ValueError, match=msg):
            run(Params(steps="cc I=100 t=60", **FAST, **kw)).write(tmp_path / "o.txt")
        return
    exe = _program(which)
    inp = tmp_path / "in.nml"
    _input(inp, tmp_path / "o.txt", dict(FAST, steps="cc I=100 t=60", **kw))
    r = subprocess.run([str(exe), str(inp)], capture_output=True, text=True, cwd=ROOT, timeout=600)
    assert r.returncode != 0 and msg in r.stdout + r.stderr, r.stdout + r.stderr
