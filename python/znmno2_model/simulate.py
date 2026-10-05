"""Running the corrected model: the protocol driver and the output table."""
from __future__ import annotations

import math

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .driver import EVENT_DV, run_protocol
from .model import CATH, MN, MO, PROBE, SO, TH, ZH, ZM, ZN, Model, sigmoid
from .params import Params

LIMIT = 1.0e-3        # a physical limit is reached within this fraction (see Stepper.limit_reason)
THETA_TOL = 1.0e-6    # a full (or empty) host is located to within this theta, as a cutoff to within EVENT_DV

COLUMNS = ("t_h", "V", "I_mAg", "mAhg", "step", "pH_cath", "pH_probe", "pH_anode", "Zn_cath_M", "Mn_cath_M",
           "S_cath_M", "vf_MnO2", "vf_ZMO", "vf_ZHS", "theta", "i_R1_mAg", "i_R2_mAg", "i_R3_mAg", "eps_cath")


@dataclass
class Result:
    rows: list = field(default_factory=list)
    exit_reason: str = ""
    steps: int = 0
    final_state: Any = None
    states: list = field(default_factory=list)      # (t_h, state) at every row, with run(..., keep_states=True)

    @property
    def array(self) -> np.ndarray:
        return np.array(self.rows, dtype=float)

    def column(self, name: str) -> np.ndarray:
        return self.array[:, COLUMNS.index(name)]

    def write(self, path) -> None:
        with open(path, "w") as fh:
            fh.write(" ".join(f"{c:>15}" for c in COLUMNS) + "\n")
            for r in self.rows:
                fh.write(" ".join(f"{int(v):>15d}" if i == 4 else f"{v:15.7E}" for i, v in enumerate(r)) + "\n")


class Stepper:
    def __init__(self, p: Params, keep_states: bool = False):
        self.p = p
        self.keep_states, self.states = keep_states, []
        self.model = Model(p)
        m = self.model.mesh
        self.vol = m.area * m.dx
        self.c = m.region == CATH
        probe = np.flatnonzero(m.region == PROBE)
        self.i_probe = int(probe[len(probe) // 2]) if len(probe) else 0

    def initial_state(self):
        self.x0 = self.model.initial_state()
        return self.x0

    def newton_step(self, x, h, I):
        return self.model.newton_step(x, h, I)

    def voltage(self, x, I):
        return self.model.voltage(x, I)

    def limit_reason(self, x, I=0.0):
        """The physical limit a state has reached when a step at current I cannot be solved, or None.
        A full host limits a discharge (I > 0), an empty one a charge (I < 0)."""
        m, p = self.model, self.p
        c = self.c
        th = sigmoid(x[c, TH])
        x0 = self.x0
        if np.min(x[:, ZN]) < LIMIT * x0[:, ZN].min():
            return "zinc_depleted"
        if p.c_MnSO4 > 0 and np.min(x[:, MN]) < LIMIT * x0[:, MN].min():
            return "manganese_depleted"
        if p.R3_on and I > 0 and th.max() > 1.0 - LIMIT:
            return "insertion_full"
        if p.R3_on and I < 0 and th.min() < LIMIT:
            return "insertion_empty"
        if np.max(x[c, MO] + x[c, ZM]) < LIMIT * max(np.max(x0[c, MO] + x0[c, ZM]), 1e-300):
            return "dissolvable_mno2_exhausted"
        if np.min(m.porosity(x)[c]) < LIMIT:
            return "pores_clogged"
        return None

    def event_margin(self, x, I):
        """Distance to a full host on discharge (theta = 1 - LIMIT) or an empty one on charge (theta = LIMIT),
        scaled so that the driver's EVENT_DV is THETA_TOL in theta; negative once crossed, inf if none applies.
        When a run ends at a host limit, its last step is redone and stopped there (driver.run_protocol): closer
        to the end the rate no longer depends on 1 - theta (or theta) and the voltage is not determined."""
        if not self.p.R3_on or I == 0:
            return math.inf
        th = sigmoid(x[self.c, TH])
        d = (1.0 - LIMIT) - th.max() if I > 0 else th.min() - LIMIT
        return d * EVENT_DV / THETA_TOL

    @staticmethod
    def finite(x):
        return bool(np.all(np.isfinite(x)))

    def row(self, t, x, mAhg, I, k):
        m, p = self.model, self.p
        if self.keep_states:
            self.states.append((t / 3600.0, x.copy()))
        free = m.free(x)
        c, vol = self.c, self.vol
        vc = vol[c]
        eps = m.porosity(x)

        def cmean(q):                       # volume average over the cathode
            return float(np.sum(vc * q) / np.sum(vc))

        i1, i2, i3, _ = m.reactions(x, free)
        to_mAg = 1.0e3 / p.mass
        return (t / 3600.0, m.voltage(x, I), I * 1.0e3 / p.mass, mAhg, k,
                -np.log10(cmean(free[c, 0])), -np.log10(free[self.i_probe, 0]), -np.log10(free[0, 0]),
                cmean(x[c, ZN]) * 1e3, cmean(x[c, MN]) * 1e3, cmean(x[c, SO]) * 1e3,
                cmean(x[c, MO] * m.V_MO), cmean(x[c, ZM] * m.V_ZM), cmean(x[c, ZH] * m.V_ZH),
                cmean(sigmoid(x[c, TH])),
                float(np.sum(vc * i1)) * to_mAg, float(np.sum(vc * i2)) * to_mAg, float(np.sum(vc * i3)) * to_mAg,
                cmean(eps[c]))


def run(p: Params, *, max_steps=None, keep_states=False) -> Result:
    """Run the protocol. keep_states=True also keeps the full state at every output row (Result.states);
    Model(p).profiles(state) turns one into profiles across the cell."""
    st = Stepper(p, keep_states)
    r = run_protocol(st, max_steps=max_steps, result=Result())
    r.states = st.states
    return r
