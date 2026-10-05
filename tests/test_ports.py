"""The Fortran and C++ programs against each other and against the Python package.

The programs are found at build/fortran/znmno2_f and build/cpp/znmno2_cpp (cmake -S . -B build && cmake
--build build), or at $ZNMNO2_FORTRAN and $ZNMNO2_CPP; tests of a missing program are skipped.

Corrected model: the Fortran and C++ output tables are byte-identical (macOS, Windows; on Linux they differ
only in round-off of near-zero values), and agree with Python to the printed precision. Faithful ports: Fortran
and C++ are byte-identical to each other (and, with ZNMNO2_ORACLE set, to the original programs' outputs); they
agree with the Python faithful port over its first 1,100 steps.
"""
import csv
import dataclasses
import os
import platform
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

from znmno2_model.namelist import EXTRA, GROUPS
from znmno2_model.params import Params

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "python" / "znmno2_model" / "data"
PROGRAMS = {"fortran": Path(os.environ.get("ZNMNO2_FORTRAN", ROOT / "build" / "fortran" / "znmno2_f")),
            "cpp": Path(os.environ.get("ZNMNO2_CPP", ROOT / "build" / "cpp" / "znmno2_cpp"))}
ORACLE = os.environ.get("ZNMNO2_ORACLE")
# platforms where the Fortran and C++ tables were found byte-identical (macOS: gfortran/clang; Windows: MinGW)
EXACT = platform.system() in ("Darwin", "Windows")


def _program(lang):
    exe = PROGRAMS[lang]
    if not exe.exists():
        pytest.skip(f"{exe} not built")
    return exe


def _lit(v):
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, str):
        return f"'{v}'"
    return repr(v)


def _input(path, out, over, mode="corrected", faithful=None):
    lines = [f"&run mode = '{mode}', data_dir = '{DATA}' /"]
    if faithful:
        lines.append("&faithful " + ", ".join(f"{k} = {v}" for k, v in faithful.items()) + " /")
    for g, names in GROUPS.items():
        items = [f"{n} = {_lit(over[n])}" for n in names if n in over]
        if g == "output":
            items = [f"file = '{out}'"] + items
        if items:
            lines.append(f"&{g}\n  " + ",\n  ".join(items) + "\n/")
    path.write_text("\n".join(lines) + "\n")


def _run(exe, inp):
    r = subprocess.run([str(exe), str(inp)], capture_output=True, text=True, cwd=ROOT, timeout=1800)
    assert r.returncode == 0, r.stdout + r.stderr
    return r.stdout


# ---------------------------------------------------------------------- parameters
def _defaults_from(path, pattern):
    return {m.group(1): m.group(2) for m in re.finditer(pattern, path.read_text(encoding="utf-8"), re.M)}


def _same(py, text):
    if isinstance(py, bool):
        return text.strip(".").lower() == str(py).lower()
    if isinstance(py, str):
        return text.strip("'\"") == py
    return float(text.replace("_dp", "")) == float(py)


def test_fortran_and_cpp_parameters_match_python():
    """Every parameter, with the same default, in fortran/params.f90 and in cpp/znmno2.cpp's Params."""
    d = Params()
    names = [f.name for f in dataclasses.fields(Params)]
    f90 = _defaults_from(ROOT / "fortran" / "params.f90",
                         r"^\s+(?:real\(dp\)|integer|logical|character\(len=512\)) :: (\w+) = ('[^']*'|\S+)")
    cpp = _defaults_from(ROOT / "cpp" / "znmno2.cpp", r"^    (?:double|int|bool|std::string) (\w+) = ([^;]+);")
    for n in names:
        assert n in f90 and _same(getattr(d, n), f90[n]), f"fortran default of {n}"
        assert n in cpp and _same(getattr(d, n), cpp[n]), f"C++ default of {n}"
    groups = dict(re.findall(r'\{"(\w+)", "(\w+)"\}', (ROOT / "cpp" / "znmno2.cpp").read_text(encoding="utf-8")))
    for g, ns in GROUPS.items():
        for n in ns:
            if n not in EXTRA:
                assert groups.get(n.lower()) == g, f"C++ group of {n}"


# ---------------------------------------------------------------------- corrected model
FAST = dict(n_probe=4, n_sep=6, n_cath=8, steps="cc I=200 t=1800; rest t=600", dt=30.0)
CASES = {
    "default": {}, "no_probe": dict(L_probe=0.0), "ph_zhs_equilibrium": dict(ph_mode="zhs_equilibrium"),
    "ph_fixed": dict(ph_mode="fixed"), "ph_spline": dict(ph_mode="spline"), "no_H": dict(species="no_H"),
    "quasi": dict(transport="quasi"), "totals": dict(basis="totals"), "zhs_equilibrium": dict(zhs="equilibrium"),
    "zhs_lumped": dict(zhs="lumped"), "zhs_off": dict(zhs="off"), "nucleation": dict(zhs_nucleation=1.05),
    "zno": dict(ZnO_on=True, logK_ZnO=8.0), "znoh2": dict(ZnOH2_on=True, logK_ZnOH2=8.0),
    "r3_nernst": dict(r3_ocp="nernst"), "r3_spline": dict(r3_ocp="spline"),
    "r1_r2_off": dict(R1_on=False, R2_on=False), "acid": dict(c_H2SO4=0.1),
    "cc_cv": dict(steps="cc I=200 Vmin=1.2; cc I=-200 Vmax=1.75; cv V=1.75 Imin=50", end_on_cutoff=False),
    "cv_sub_steps": dict(steps="cc I=200 Vmin=1.2; cv V=1.9 Imin=5", end_on_cutoff=False, dt=600.0),
    # ends at a full host: the exit row is located at theta = 1 - 1e-3, where the voltage is determined (#14)
    "host_full": dict(R1_on=False, R2_on=False, steps="cc I=400 Vmin=0.2"),
}
_PY = {}


def _python_table(name, tmp_path_factory):
    if name not in _PY:
        from znmno2_model.simulate import run
        r = run(Params(**{**FAST, **CASES[name]}))
        out = tmp_path_factory.mktemp("py") / f"{name}.txt"
        r.write(out)
        _PY[name] = (out, r.exit_reason)
    return _PY[name]


def _port_table(lang, name, tmp_path):
    out, inp = tmp_path / f"{name}_{lang}.txt", tmp_path / f"{name}_{lang}.nml"
    _input(inp, out, {**FAST, **CASES[name]})
    msg = _run(_program(lang), inp)
    return out, re.search(r"exit '(\w+)'", msg).group(1)


@pytest.mark.parametrize("name", sorted(CASES))
@pytest.mark.parametrize("lang", ["fortran", "cpp"])
def test_corrected_agrees_with_python(lang, name, tmp_path, tmp_path_factory):
    out, reason = _port_table(lang, name, tmp_path)
    ref, ref_reason = _python_table(name, tmp_path_factory)
    assert reason == ref_reason
    a_txt, b_txt = out.read_text().splitlines(), ref.read_text().splitlines()
    assert a_txt[0] == b_txt[0] and len(a_txt) == len(b_txt)
    a, b = np.loadtxt(out, skiprows=1, ndmin=2), np.loadtxt(ref, skiprows=1, ndmin=2)
    # the printed precision (8 digits) relative to each column's scale; near-zero columns at an absolute floor
    tol = 1e-6 * np.abs(b).max(axis=0) + 1e-15
    assert np.all(np.abs(a - b) <= tol), np.argwhere(np.abs(a - b) > tol)[:5]


@pytest.mark.parametrize("name", sorted(CASES))
def test_corrected_fortran_equals_cpp(name, tmp_path):
    f, _ = _port_table("fortran", name, tmp_path)
    c, _ = _port_table("cpp", name, tmp_path)
    fl, cl = f.read_text().splitlines(), c.read_text().splitlines()
    assert len(fl) == len(cl) and fl[0] == cl[0]
    if not EXACT:
        # Linux: the two compilers leave different round-off in quantities that cancel to zero (e.g. a ZHS
        # volume fraction of 1e-26 against 2e-17); everything else is identical
        a, b = np.loadtxt(f, skiprows=1, ndmin=2), np.loadtxt(c, skiprows=1, ndmin=2)
        assert np.all(np.abs(a - b) <= 1e-12 * np.abs(b).max(axis=0) + 1e-15)
        return
    for i, (x, y) in enumerate(zip(fl, cl)):
        cols = [k for k, (u, v) in enumerate(zip(x.split(), y.split())) if u != v]
        assert x == y, f"row {i}, columns {cols}:\nfortran {x}\nC++     {y}"


def test_corrected_with_a_phreeqc_database_and_log_k_overrides(tmp_path):
    """equilibria_db (tabs between fields, a redox species and a foreign element to leave out) and logk_file."""
    db = tmp_path / "mini.dat"
    db.write_text("SOLUTION_SPECIES\nH+ = H+\n\tlog_k\t0\ne- = e-\n\tlog_k\t0\nH2O = H2O\nZn+2 = Zn+2\n"
                  "Mn+2 = Mn+2\nSO4-2 = SO4-2\nH2O = OH- + H+\n\t-analytic\t293.29227\t0.1360833\t-10576.913"
                  "\t-123.73158\t0\t-6.996455e-5\n2 H2O = O2 + 4 H+ + 4 e-\n\tlog_k\t-86.08\n"
                  "SO4-2 + H+ = HSO4-\n\tlog_k\t1.988\nZn+2 + H2O = ZnOH+ + H+\n\tlog_k\t-8.96\n"
                  "Zn+2 + SO4-2 = ZnSO4\n\tlog_k\t2.37\nMn+2 + SO4-2 = MnSO4\n\tlog_k\t2.25\n"
                  "Mn+2 + Cl- = MnCl+\n\tlog_k\t0.61\nPHASES\n")
    lk = tmp_path / "logk.txt"
    lk.write_text("# overrides\nZnSO4 2.1\n")
    over = {**FAST, "equilibria_db": str(db), "logk_file": str(lk), "transport": "quasi"}
    from znmno2_model.simulate import run
    r = run(Params(**over))
    ref = tmp_path / "py.txt"
    r.write(ref)
    outs = []
    for lang in ("fortran", "cpp"):
        if not PROGRAMS[lang].exists():
            continue
        out, inp = tmp_path / f"{lang}.txt", tmp_path / f"{lang}.nml"
        _input(inp, out, over)
        _run(PROGRAMS[lang], inp)
        a, b = np.loadtxt(out, skiprows=1), np.loadtxt(ref, skiprows=1)
        assert a.shape == b.shape and np.all(np.abs(a - b) <= 1e-6 * np.abs(b).max(axis=0) + 1e-15)
        outs.append(out.read_bytes())
    if not outs:
        pytest.skip("no program built")
    assert all(o == outs[0] for o in outs)


# ---------------------------------------------------------------------- faithful ports
CHARGE = dict(rxnk_2="-8.5", rxnk_3="-7.5", frac_zmcx="0.4", frac_zmcmax="0.03", xmax_t="0.022",
              applied_current="0.000121", porosity="0.8", volfrac_mno2="0.07", stated_mass_loading="0.00121",
              zhs_ksp="30")
PHCELL = dict(rxnk_5="-9.0", fraction_kmno2="0.7")


def _faithful(lang, line, values, tmp_path):
    out, inp = tmp_path / f"{line}_{lang}.txt", tmp_path / f"{line}_{lang}.nml"
    _input(inp, out, {}, mode=f"faithful_{line}", faithful=values)
    _run(_program(lang), inp)
    return out.read_text().splitlines()


@pytest.mark.parametrize("line", ["charge", "phcell"])
def test_faithful_fortran_equals_cpp(line, tmp_path):
    values = CHARGE if line == "charge" else PHCELL
    assert _faithful("fortran", line, values, tmp_path) == _faithful("cpp", line, values, tmp_path)


@pytest.mark.parametrize("line", ["charge", "phcell"])
@pytest.mark.parametrize("lang", ["fortran", "cpp"])
def test_faithful_agrees_with_python(lang, line, tmp_path):
    """The first 1,100 steps (ramp and the first 100 s) against the Python faithful port, which uses the
    installed bandsolver (with fused multiply-adds it differs from the original in the last digits)."""
    from znmno2_model.faithful.charge import ChargeModel
    from znmno2_model.faithful.constants import Template
    from znmno2_model.faithful.phcell import PhCellModel, PhTemplate
    rows = _faithful(lang, line, CHARGE if line == "charge" else PHCELL, tmp_path)
    if line == "charge":
        m = ChargeModel(Template(**{k.upper() if k != "xmax_t" else "XMAX": v for k, v in CHARGE.items()}))
    else:
        m = PhCellModel(PhTemplate(**{k.upper(): v for k, v in PHCELL.items()}))
    assert m.run(max_steps=1100) == "max_steps"
    assert rows[0] == m.rows[0]
    num = lambda r: np.array([float(v) for i, v in enumerate(r.split()) if i != 4])
    for got, want in zip(rows[1:len(m.rows)], m.rows[1:]):
        a, b = num(got), num(want)
        assert np.allclose(a, b, rtol=1e-6, atol=1e-12), (got, want)


@pytest.mark.skipif(not ORACLE, reason="ZNMNO2_ORACLE not set (private reference outputs)")
@pytest.mark.parametrize("lang", ["fortran", "cpp"])
@pytest.mark.parametrize("line", ["charge", "phcell"])
def test_faithful_reproduces_the_original_programs(lang, line, tmp_path):
    """Every reference run, start to finish, byte for byte."""
    path = next((Path(ORACLE) / line).glob("simulation_params_*.txt"))
    for row in csv.DictReader(path.open()):
        k = row["Sim_label"]
        ref = Path(ORACLE) / "local" / line / k / "Time_Voltage.txt"
        if not ref.exists():
            continue
        vals = {("xmax_t" if c.strip("_") == "XMAX" else c.strip("_").lower()): v
                for c, v in row.items() if c.startswith("_")}
        rows = _faithful(lang, line, vals, tmp_path)
        assert rows == ref.read_text().splitlines(), f"{line} run {k}"


# ---------------------------------------------------------------------- bad input
BAD = [
    ("&options\n  ph_mode = 'speciatio'\n/\n", "ph_mode must be one of"),
    ("&cell\n  n_sep = 6, foo = 1.0\n/\n", ("foo", "&cell")),                        # unknown name
    ("&options\n  n_sep = 6\n/\n", ("n_sep", "&options")),                          # name in the wrong group
    ("&cellz\n  n_sep = 6\n/\n", "unknown namelist group &cellz"),
    ("&cell\n  n_sep = 6\n/\n&cell\n  n_cath = 8\n/\n", "namelist group &cell appears twice"),
    ("&protocol\n  steps = 'cx I=1'\n/\n", "unknown step type 'cx'"),
]


@pytest.mark.parametrize("text,msg", BAD)
@pytest.mark.parametrize("lang", ["fortran", "cpp"])
def test_programs_reject_bad_input(lang, text, msg, tmp_path):
    """Each program stops with a non-zero exit and a message naming the problem (#10)."""
    nml = tmp_path / "bad.nml"
    nml.write_text(f"&run data_dir = '{DATA}' /\n" + text + "&output file = 'o.txt' /\n")
    r = subprocess.run([str(_program(lang)), str(nml)], cwd=tmp_path, capture_output=True, text=True)
    msgs = msg if isinstance(msg, tuple) else (msg,)
    assert r.returncode != 0 and all(m in r.stdout + r.stderr for m in msgs), r.stdout + r.stderr


def test_python_rejects_a_repeated_group(tmp_path):
    from znmno2_model.namelist import load
    nml = tmp_path / "bad.nml"
    nml.write_text("&cell\n  n_sep = 6\n/\n&cell\n  n_cath = 8\n/\n")
    with pytest.raises(ValueError, match="namelist group &cell appears twice"):
        load(nml)
