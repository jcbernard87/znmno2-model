"""Faithful mode against the original programs' output (byte for byte).

The reference outputs and parameter rows are not distributed. Set ZNMNO2_ORACLE to the private folder
holding the parameter files (charge/simulation_params_*.txt, phcell/simulation_params_*.txt) and the
outputs rebuilt from the original sources on this machine (local/<line>/<k>/Time_Voltage.txt).
Without it these tests are skipped. Each test runs the first 1,100 steps: the current ramp and the
first 100 s of discharge.
"""
import csv
import os
from pathlib import Path

import pytest

ORACLE = os.environ.get("ZNMNO2_ORACLE")
pytestmark = pytest.mark.skipif(not ORACLE, reason="ZNMNO2_ORACLE not set (private reference outputs)")

STEPS = 1100


def _params(line, label):
    path = next((Path(ORACLE) / line).glob("simulation_params_*.txt"))
    row = next(r for r in csv.DictReader(path.open()) if r["Sim_label"] == str(label))
    return {k.strip("_"): v for k, v in row.items() if k.startswith("_")}


def _compare(model, line, label):
    assert model.run(max_steps=STEPS) == "max_steps"
    ref = (Path(ORACLE) / "local" / line / str(label) / "Time_Voltage.txt").read_text().splitlines()
    assert len(model.rows) > 2
    for i, row in enumerate(model.rows):
        assert row == ref[i], f"row {i} differs:\n{row}\n{ref[i]}"


@pytest.mark.parametrize("label", [1, 4])
def test_charge_first_rows(label):
    from znmno2_model.faithful.charge import ChargeModel
    from znmno2_model.faithful.constants import Template
    _compare(ChargeModel(Template(**_params("charge", label))), "charge", label)


@pytest.mark.parametrize("label", [1, 4])
def test_phcell_first_rows(label):
    from znmno2_model.faithful.phcell import PhCellModel, PhTemplate
    _compare(PhCellModel(PhTemplate(**_params("phcell", label))), "phcell", label)


@pytest.mark.parametrize("line,label", [("charge", 1), ("phcell", 1)])
def test_public_solver_matches_to_round_off(line, label, monkeypatch):
    """Without the private BAND/MATINV transliteration, the faithful ports (bandsolver's legacy pivot) agree
    with the original output to round-off, not byte for byte."""
    import numpy as np
    from znmno2_model.faithful import solver
    monkeypatch.setattr(solver, "_LEGACY", None)
    if line == "charge":
        from znmno2_model.faithful.charge import ChargeModel
        from znmno2_model.faithful.constants import Template
        m = ChargeModel(Template(**_params("charge", label)))
    else:
        from znmno2_model.faithful.phcell import PhCellModel, PhTemplate
        m = PhCellModel(PhTemplate(**_params("phcell", label)))
    assert m.run(max_steps=STEPS) == "max_steps"
    ref = (Path(ORACLE) / "local" / line / str(label) / "Time_Voltage.txt").read_text().splitlines()
    num = lambda r: np.array([float(v) for i, v in enumerate(r.split()) if i != 4])
    for a, b in zip(m.rows[1:], ref[1:len(m.rows)]):
        x, y = num(a), num(b)
        assert np.allclose(x, y, rtol=1e-9, atol=1e-300), np.max(np.abs(x - y) / np.maximum(np.abs(y), 1e-300))


def test_private_solver_is_in_use():
    from znmno2_model.faithful import solver
    assert solver.EXACT, "ZNMNO2_ORACLE has no band_legacy.py: byte-for-byte tests would be meaningless"
