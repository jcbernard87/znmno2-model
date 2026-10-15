"""Restart files (#16): a run continued from a restart repeats the run it continues."""
import numpy as np
import pytest

from znmno2_model.params import Params
from znmno2_model.simulate import run

FAST = dict(n_probe=4, n_sep=6, n_cath=8, dt=30.0)
S1, S2 = "cc I=200 t=1800", "rest t=600; cc I=-200 Vmax=1.75; cv V=1.75 Imin=50"
S2_CV = "cv V=1.3 Imin=20"                          # a constant-voltage first step uses the stored current


def _split(tmp_path, s2, **kw):
    whole = run(Params(steps=f"{S1}; {s2}", end_on_cutoff=False, **FAST, **kw)).array
    f = tmp_path / "r.txt"
    run(Params(steps=S1, end_on_cutoff=False, restart_write=str(f), **FAST, **kw))
    part = run(Params(steps=s2, end_on_cutoff=False, restart_read=str(f), **FAST, **kw)).array
    step = 4                                         # the step-index column (COLUMNS)
    n1 = int(np.sum(whole[:, step] == 1))
    return whole[n1:], part[1:], step


@pytest.mark.parametrize("s2", [S2, S2_CV], ids=["rest_cc_cv", "cv_first"])
def test_restart_continues_the_run_exactly(tmp_path, s2):
    a, b, step = _split(tmp_path, s2)
    assert a.shape == b.shape
    cols = [c for c in range(a.shape[1]) if c != step]
    np.testing.assert_array_equal(a[:, cols], b[:, cols])


def test_restart_must_match_the_mesh(tmp_path):
    f = tmp_path / "r.txt"
    run(Params(steps="cc I=100 t=60", restart_write=str(f), **FAST))
    with pytest.raises(ValueError, match="does not match the input: array x"):
        run(Params(steps="rest t=60", restart_read=str(f), **{**FAST, "n_cath": 9}))


@pytest.mark.parametrize("lang", ["fortran", "cpp"])
def test_native_restart_continues_the_run_exactly(lang, tmp_path):
    """The Fortran and C++ programs: S1 then S2 from a restart equals S1; S2 in one run, byte for byte (step index
    aside); Python continues from the program's restart file to the printed precision."""
    from test_ports import _input, _program, _run
    exe = _program(lang)
    over = dict(FAST, end_on_cutoff=False)

    def go(name, steps, **out):
        inp, res = tmp_path / f"{name}.nml", tmp_path / f"{name}.txt"
        _input(inp, res, {**over, "steps": steps, **out})
        _run(exe, inp)
        return [ln.split() for ln in res.read_text().splitlines()[1:]]

    rf = tmp_path / f"{lang}_r.txt"
    whole = go("w", f"{S1}; {S2}")
    go("a", S1, restart_write=str(rf))
    part = go("b", S2, restart_read=str(rf))
    n1 = sum(1 for r in whole if r[4] == "1")
    drop = lambda r: r[:4] + r[5:]
    assert [drop(r) for r in whole[n1:]] == [drop(r) for r in part[1:]]
    py = run(Params(steps=S2, restart_read=str(rf), **over)).array
    b = np.array([[float(x) for x in r] for r in part])
    np.testing.assert_allclose(py[:, :4], b[:, :4], rtol=1e-6, atol=1e-6)
