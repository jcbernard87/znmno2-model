"""Faithful ports with the public solver (bandsolver's legacy pivot): they run, with illustrative inputs."""

import numpy as np

from znmno2_model.faithful import solver
from znmno2_model.faithful.charge import ChargeModel
from znmno2_model.faithful.constants import Template
from znmno2_model.faithful.phcell import PhCellModel, PhTemplate


def _numbers(rows):
    return np.array([[float(v) for i, v in enumerate(r.split()) if i != 4] for r in rows[1:]])


def test_charge_line_runs(monkeypatch):
    monkeypatch.setattr(solver, "_LEGACY", None)
    t = Template(RXNK_2="-8.5", RXNK_3="-7.5", FRAC_ZMCX="0.4", FRAC_ZMCMAX="0.03", XMAX="0.022",
                 APPLIED_CURRENT="0.000121", POROSITY="0.8", VOLFRAC_MNO2="0.07", STATED_MASS_LOADING="0.00121",
                 ZHS_KSP="30")
    m = ChargeModel(t)
    assert m.run(max_steps=1100) == "max_steps"
    a = _numbers(m.rows)
    assert len(a) > 5 and np.all(np.isfinite(a[:, :4]))


def test_phcell_line_runs(monkeypatch):
    monkeypatch.setattr(solver, "_LEGACY", None)
    m = PhCellModel(PhTemplate(RXNK_5="-9.0", FRACTION_KMNO2="0.7"))
    assert m.run(max_steps=1100) == "max_steps"
    a = _numbers(m.rows)
    assert len(a) > 5 and np.all(np.isfinite(a[:, :4]))
