"""Mesh and time-step convergence of the corrected model (docs/validation.md).

Runs the default discharge at mesh factors 1, 2, 4 (dt = 2.5 s) and at dt = 10, 5, 2.5, 1.25 s (mesh 1),
and prints the voltage at fixed capacities with the observed order p = log2(e1/e2).

usage: python scripts/convergence.py
"""
import numpy as np

from znmno2_model.params import Params
from znmno2_model.simulate import run

Q = (25.0, 50.0, 100.0, 140.0, 170.0)


def voltages(mesh=1, dt=10.0):
    p = Params(n_probe=20 * mesh, n_sep=30 * mesh, n_cath=40 * mesh, dt=dt, write_interval=dt)
    a = run(p).array
    return np.array([np.interp(q, a[:, 3], a[:, 1]) for q in Q])


def table(title, labels, V):
    print(title)
    for q, col in zip(Q, V.T):
        e1, e2 = abs(col[0] - col[1]), abs(col[1] - col[2])
        print(f"  {q:5.0f} mAh/g: " + "  ".join(f"{lab}: {v:.9f}" for lab, v in zip(labels, col))
              + f"   order {np.log2(e1 / e2):.2f}")


if __name__ == "__main__":
    table("space (dt = 2.5 s)", ["x1", "x2", "x4"], np.array([voltages(m, 2.5) for m in (1, 2, 4)]))
    table("time (mesh x1)", ["10 s", "5 s", "2.5 s"], np.array([voltages(1, dt) for dt in (10.0, 5.0, 2.5)]))
