"""Speciation for the corrected model: free-ion concentrations from the four transported totals.

Ideal activities (concentrations in mol/L) and the 21 complexation equilibria of Herrmann et al.,
Adv. Energy Mater. (2023), Table S1 (speciation.HERRMANN_2023). The unknowns are the log10 free
concentrations of the masters H+, Zn2+, Mn2+ and SO4 2-; the equations are the balances on the totals
of the masters:

    H_T  = sum_s nu_sH  m_s   (master protons: H+, HSO4-, 2 H2SO4, minus OH- and hydroxo complexes)
    Zn_T = sum_s nu_sZn m_s,  Mn_T = ...,  S_T = ...

H_T can be negative (more hydroxide than protons). The charge of a solution of totals is
H_T + 2 Zn_T + 2 Mn_T - 2 S_T, which the model's electroneutrality keeps at zero.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .speciation import MASTERS, Activity, Database, Speciation, herrmann_2023

LN10 = np.log(10.0)


def read_logk_file(path) -> dict:
    """log K overrides: one 'species_name  log10_K' pair per line ('#' starts a comment)."""
    out = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#")[0].strip()
        if line:
            name, val = line.split()
            out[name] = float(val)
    return out


class Equilibria:
    def __init__(self, logk_overrides: dict | None = None, database: str = ""):
        """database: '' for Herrmann et al. (2023) Table S1, or the path of a PHREEQC database whose
        SOLUTION_SPECIES made of H+, Zn+2, Mn+2, SO4-2 and H2O are used (e.g. phreeqc.dat)."""
        db = herrmann_2023() if not database else Database.from_phreeqc(database)
        sp = Speciation(db, Activity(kind="ideal"))
        self.names = [s.name for s in sp.species]
        self.nu = sp.nu.copy()                      # (n_species, 4): masters H+, Zn+2, Mn+2, SO4-2
        self.logk = sp.logk.copy()
        for name, val in (logk_overrides or {}).items():
            if name not in self.names:
                raise ValueError(f"log K override for unknown species {name!r}")
            self.logk[self.names.index(name)] = val
        self.z = sp.z.copy()
        self.abs_nu_h = np.abs(self.nu[:, 0])
        assert tuple(self.names[:4]) == MASTERS

    def species(self, x):
        """Concentrations of every species (..., n_species) [mol/L] from x = log10 master concentrations."""
        return 10.0 ** np.clip(self.logk + x @ self.nu.T, -300.0, 300.0)

    def totals(self, x):
        """Totals (..., 4) [H_T, Zn_T, Mn_T, S_T] of a speciation x."""
        return self.species(x) @ self.nu

    def _residual(self, x, T):
        m = self.species(x)
        bal = m @ self.nu
        d = m @ self.abs_nu_h
        F = np.empty_like(x)
        F[..., 0] = (bal[..., 0] - T[..., 0]) / d
        F[..., 1:] = np.log10(bal[..., 1:]) - np.log10(T[..., 1:])
        return F, m, bal, d

    def _jacobian(self, m, bal, d, T):
        dm = m[..., :, None] * self.nu * LN10                         # d m_s / d x_j
        J = np.empty(m.shape[:-1] + (4, 4))
        db = np.einsum("sk,...sj->...kj", self.nu, dm)                # d bal_k / d x_j
        dd = np.einsum("s,...sj->...j", self.abs_nu_h, dm)
        J[..., 0, :] = (db[..., 0, :] * d[..., None] - (bal[..., 0] - T[..., 0])[..., None] * dd) / (d ** 2)[..., None]
        J[..., 1:, :] = db[..., 1:, :] / (bal[..., 1:, None] * LN10)
        return J

    def initial_guess(self, T):
        T = np.asarray(T, float)
        x = np.empty(T.shape)
        x[..., 0] = np.log10(np.maximum(T[..., 0], 0.0) + 1e-5)
        x[..., 1:] = np.log10(np.maximum(T[..., 1:], 1e-300)) - 0.5
        return x

    def _inner(self, xh, T, x, tol=1e-13):
        """Zn, Mn and S balances at fixed log10 c_H (xh): Newton on log10 c_Zn, c_Mn, c_SO4 (monotone in
        each; well behaved from any start). Returns x with x[..., 0] = xh."""
        x = np.array(x, float)
        x[..., 0] = xh
        for _ in range(200):
            F, m, bal, d = self._residual(x, T)
            f = F[..., 1:]
            if np.all(np.abs(f) < tol):
                break
            J = self._jacobian(m, bal, d, T)[..., 1:, 1:]
            dx = np.linalg.solve(J, -f[..., None])[..., 0]
            dx = np.where(np.isfinite(dx), dx, 0.0)
            big = np.max(np.abs(dx), axis=-1, keepdims=True)
            x[..., 1:] += dx * np.minimum(1.0, 1.0 / np.maximum(big, 1e-300))
        return x

    def _bracketed(self, T, x):
        """Robust fallback: bisection on log10 c_H in [-16, 2] (the proton balance increases with c_H),
        with the other balances solved at each trial pH; then a Newton polish."""
        lo = np.full(T.shape[:-1], -16.0)
        hi = np.full(T.shape[:-1], 2.0)
        x = self.initial_guess(T)
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            x = self._inner(mid, T, x)
            g = self.species(x) @ self.nu[:, 0] - T[..., 0]
            lo, hi = np.where(g < 0, mid, lo), np.where(g < 0, hi, mid)
            if np.max(hi - lo) < 1e-12:
                break
        return self._inner(0.5 * (lo + hi), T, x)

    def solve(self, T, x0=None, tol=1e-13, max_iter=100):
        """Free-ion log10 concentrations (..., 4) for totals T (..., 4) = [H_T, Zn_T, Mn_T, S_T] in mol/L.

        Newton with the step scaled to at most one decade and backtracking on |F|. Zn_T, Mn_T and S_T
        must be positive. Raises RuntimeError if some solution does not converge.
        """
        T = np.array(T, float)
        T[..., 1:] = np.maximum(T[..., 1:], 1e-30)                    # an absent element: a trace amount
        x = self.initial_guess(T) if x0 is None else np.array(x0, float)
        F, m, bal, d = self._residual(x, T)
        for _ in range(max_iter):
            err = np.max(np.abs(F), axis=-1)
            if np.all(err < tol):
                return x
            J = self._jacobian(m, bal, d, T)
            dx = np.linalg.solve(J, -F[..., None])[..., 0]
            dx = np.where(np.isfinite(dx), dx, 0.0)
            big = np.max(np.abs(dx), axis=-1, keepdims=True)
            dx = dx * np.minimum(1.0, 1.0 / np.maximum(big, 1e-300))
            dx = np.where((err >= tol)[..., None], dx, 0.0)
            f0 = np.linalg.norm(F, axis=-1)
            lam = np.ones(err.shape)
            for _ in range(40):
                F1, m1, bal1, d1 = self._residual(x + lam[..., None] * dx, T)
                worse = ~(np.linalg.norm(F1, axis=-1) < f0) & (err >= tol)
                if not worse.any():
                    break
                lam = np.where(worse, 0.5 * lam, lam)
            x = x + lam[..., None] * dx
            F, m, bal, d = F1, m1, bal1, d1
        bad = ~(np.max(np.abs(F), axis=-1) < 1e-10)
        if not bad.any():
            return x
        # Newton lost its way (the scaled proton balance saturates far from the solution): bracket the pH
        xb = self._bracketed(T[bad], x[bad])
        F2 = self._residual(xb, T[bad])[0]
        if not np.all(np.max(np.abs(F2), axis=-1) < 1e-10):
            raise RuntimeError("speciation did not converge")
        x[bad] = xb
        return x

    def sensitivity(self, x, T):
        """d x / d T (..., 4, 4) at a converged speciation x of totals T (implicit function theorem)."""
        T = np.array(T, float)
        T[..., 1:] = np.maximum(T[..., 1:], 1e-30)
        F, m, bal, d = self._residual(x, T)
        J = self._jacobian(m, bal, d, T)
        dFdT = np.zeros(J.shape)
        dFdT[..., 0, 0] = -1.0 / d
        for k in range(1, 4):
            dFdT[..., k, k] = -1.0 / (T[..., k] * LN10)
        return -np.linalg.solve(J, dFdT)
