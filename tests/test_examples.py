"""The example inputs in input/examples are valid: the corrected-mode ones are the default parameters with the
changes their comments state, and the Fortran program runs them (GITT, 25 s, is only loaded)."""
import os
import re
import subprocess
from pathlib import Path

import pytest

from znmno2_model.namelist import load
from znmno2_model.params import Params

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "input" / "examples"
FORTRAN = Path(os.environ.get("ZNMNO2_FORTRAN", ROOT / "build" / "fortran" / "znmno2_f"))

CORRECTED = {
    "gitt": dict(steps="cc I=50 t=720; rest t=1800", cycles=200),
    "cccv": dict(steps="cc I=100 Vmin=1.0; rest t=600; cc I=-100 Vmax=1.8; cv V=1.8 Imin=5", end_on_cutoff=False),
    "acid_probe": dict(c_H2SO4=0.1, L_probe=0.2),
}
EXPECT = {"cccv": "exit 'end_of_protocol'", "acid_probe": "exit 'cutoff_low'",
          "faithful_charge": "NaN in cprev", "faithful_phcell": "EXIT BECAUSE LOWER VOLTAGE CUTOFF"}


@pytest.mark.parametrize("name", sorted(CORRECTED))
def test_corrected_examples_are_the_defaults_plus_their_changes(name):
    p, _ = load(EX / f"{name}.nml")
    assert p == Params(**CORRECTED[name])


@pytest.mark.parametrize("name", sorted(EXPECT))
def test_the_fortran_program_runs_the_examples(name, tmp_path):
    if not FORTRAN.exists():
        pytest.skip(f"{FORTRAN} not built")
    text = (EX / f"{name}.nml").read_text(encoding="utf-8")
    data, out = ROOT / "python" / "znmno2_model" / "data", tmp_path / "out.txt"
    # function replacements: a Windows path's backslashes are not escapes
    text = re.sub(r"data_dir = '[^']*'", lambda m: f"data_dir = '{data}'", text)
    text = re.sub(r"^(\s*)file = '[^']*'", lambda m: f"{m.group(1)}file = '{out}'", text, flags=re.M)
    if "&run" not in text:
        text = f"&run data_dir = '{data}' /\n" + text
    (tmp_path / "in.nml").write_text(text)
    r = subprocess.run([str(FORTRAN), str(tmp_path / "in.nml")], capture_output=True, text=True, cwd=tmp_path)
    assert r.returncode == 0 and EXPECT[name] in r.stdout, r.stdout + r.stderr
