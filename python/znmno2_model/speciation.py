"""Aqueous speciation of the Zn/Mn sulfate electrolyte: pH from the equilibrium equations, solved directly.

This replaces the original program's 2-D pH spline, which was a fit to PHREEQC runs (pH as a function of
Zn and Mn only). Here the same equilibria are solved at every node, for any composition:

  unknowns (per node):  log10 activities of the master species H+, Zn+2, Mn+2, SO4-2
  equations:            mass balances of Zn, Mn and S, and electroneutrality
  species:              log10 a_s = log10 K_s + sum_i nu_si log10 a_i  (nu over master species and H2O)
  activity:             PHREEQC's model (phreeqc.dat): WATEQ Debye-Hueckel for species with -gamma a b,
                        Davies otherwise, log10 gamma = 0.1 I for uncharged species

The species set and constants are read from a PHREEQC database file (the subset of elements used here).
Redox species (those whose reaction involves e-) are left out: in these solutions they are below 1e-15 mol/kgw.

Concentrations are molal (mol/kg water), as PHREEQC uses them.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

MASTERS = ("H+", "Zn+2", "Mn+2", "SO4-2")      # order of the unknowns
CHARGE_OF = {"H+": 1, "Zn+2": 2, "Mn+2": 2, "SO4-2": -2, "H2O": 0}
LN10 = math.log(10.0)


def _charge(name: str) -> int:
    m = re.search(r"([+-])(\d*)$", name)
    if not m:
        return 0
    return (1 if m.group(1) == "+" else -1) * (int(m.group(2)) if m.group(2) else 1)


def _parse_side(side: str) -> list[tuple[float, str]]:
    terms = []
    for tok in side.split(" + "):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^(\d*\.?\d*)\s*(.+)$", tok)
        coef = float(m.group(1)) if m.group(1) else 1.0
        terms.append((coef, m.group(2).strip()))
    return terms


@dataclass
class Species:
    name: str
    nu: dict                    # master species (incl. H2O) -> stoichiometric coefficient
    logk: float
    z: int
    gamma: tuple | None = None  # (a, b) of -gamma, or None (Davies / neutral)


@dataclass
class Database:
    species: list = field(default_factory=list)
    master_gamma: dict = field(default_factory=dict)   # -gamma (a, b) of the master species

    @classmethod
    def from_phreeqc(cls, path, temperature_c: float = 25.0, overrides: dict | None = None) -> "Database":
        """Read the SOLUTION_SPECIES of a PHREEQC database for species made of H+, H2O, Zn+2, Mn+2, SO4-2.

        `overrides` maps a reaction string (as written in the database) to a log_k value.
        """
        T = temperature_c + 273.15
        text = Path(path).read_text(errors="replace")
        block = re.search(r"^SOLUTION_SPECIES\s*$(.*?)^PHASES\s*$", text, re.M | re.S).group(1)
        allowed = set(MASTERS) | {"H2O"}
        found: dict[str, Species] = {}
        cur = None

        master_gamma: dict[str, tuple] = {}

        def finish(c):
            if c is None:
                return
            if c["ok"]:
                found[c["name"]] = Species(c["name"], c["nu"], c["logk"], _charge(c["name"]), c["gamma"])
            elif c.get("identity") and c["gamma"] is not None:
                master_gamma[c["name"]] = c["gamma"]

        for raw in block.splitlines():
            line = raw.split("#")[0].rstrip()
            if not line.strip():
                continue
            if "=" in line and not line.lstrip().startswith("-") and not re.match(r"^\s*log_k", line):
                finish(cur)
                lhs, rhs = line.split("=")
                L, R = _parse_side(lhs), _parse_side(rhs)
                name = R[0][1]
                nu: dict[str, float] = {}
                ok = R[0][0] == 1.0
                for c, sp in L:
                    nu[sp] = nu.get(sp, 0.0) + c
                for c, sp in R[1:]:
                    nu[sp] = nu.get(sp, 0.0) - c
                ok = ok and all(sp in allowed for sp in nu) and name not in allowed and "e-" not in line
                identity = name in allowed and len(L) == 1 and L[0][1] == name
                if name in allowed:                      # identity / master definitions
                    ok = False
                key = re.sub(r"\s+", " ", line.strip())
                cur = dict(name=name, nu=nu, logk=0.0, gamma=None, ok=ok, key=key, identity=identity)
                if overrides and key in overrides:
                    cur["logk"] = overrides[key]
                    cur["fixed"] = True
                continue
            if cur is None:
                continue
            opt = line.strip().split()
            o = opt[0].lstrip("-").lower()
            if o in ("log_k", "logk") and not cur.get("fixed"):
                cur["logk"] = float(opt[1])
            elif o.startswith("analytic") and not cur.get("fixed"):
                a = [float(x) for x in opt[1:]] + [0.0] * 6
                cur["logk"] = a[0] + a[1] * T + a[2] / T + a[3] * math.log10(T) + a[4] / T ** 2 + a[5] * T ** 2
            elif o == "gamma":
                cur["gamma"] = (float(opt[1]), float(opt[2]))
        finish(cur)
        return cls(list(found.values()), master_gamma)


# Herrmann, Euchner, Gross, Horstmann, "The Cycling Mechanism of Manganese-Oxide Cathodes in Zinc Batteries:
# A Theory-Based Approach", Adv. Energy Mater. (2023), Supporting Information Table S1 (log10 beta, 25 C).
# Reactions are written as in the table; "2 Mn+2 + 3 H2O = Mn2(OH)3+3 + H+" in the table is unbalanced and is
# taken as 2 Mn+2 + H2O = Mn2OH+3 + H+ (log10 beta -10.56).
HERRMANN_2023 = [
    ("OH-", {"H2O": 1, "H+": -1}, -14.0),          # H+ + OH- = H2O, log beta = -14.0 as tabulated (i.e. Kw)
    ("H2SO4", {"H+": 2, "SO4-2": 1}, 0.0),
    ("HSO4-", {"H+": 1, "SO4-2": 1}, 1.98),
    ("ZnOH+", {"Zn+2": 1, "H2O": 1, "H+": -1}, -7.5),
    ("Zn(OH)2", {"Zn+2": 1, "H2O": 2, "H+": -2}, -16.4),
    ("Zn(OH)3-", {"Zn+2": 1, "H2O": 3, "H+": -3}, -28.2),
    ("Zn(OH)4-2", {"Zn+2": 1, "H2O": 4, "H+": -4}, -41.3),
    ("Zn2OH+3", {"Zn+2": 2, "H2O": 1, "H+": -1}, -9.0),
    ("Zn2(OH)6-2", {"Zn+2": 2, "H2O": 6, "H+": -6}, -54.3),
    ("Zn4(OH)4+4", {"Zn+2": 4, "H2O": 4, "H+": -4}, -27.0),
    ("ZnSO4", {"Zn+2": 1, "SO4-2": 1}, 2.37),
    ("Zn(SO4)2-2", {"Zn+2": 1, "SO4-2": 2}, 3.28),
    ("Zn(SO4)3-4", {"Zn+2": 1, "SO4-2": 3}, 1.7),
    ("Zn(SO4)4-6", {"Zn+2": 1, "SO4-2": 4}, 1.7),
    ("MnOH+", {"Mn+2": 1, "H2O": 1, "H+": -1}, -10.59),
    ("Mn(OH)2", {"Mn+2": 1, "H2O": 2, "H+": -2}, -18.54),
    ("Mn(OH)3-", {"Mn+2": 1, "H2O": 3, "H+": -3}, -34.8),
    ("Mn(OH)4-2", {"Mn+2": 1, "H2O": 4, "H+": -4}, -48.3),
    ("Mn2(OH)3+", {"Mn+2": 2, "H2O": 3, "H+": -3}, -23.9),
    ("Mn2OH+3", {"Mn+2": 2, "H2O": 1, "H+": -1}, -10.56),
    ("MnSO4", {"Mn+2": 1, "SO4-2": 1}, 2.25),
]


def herrmann_2023() -> "Database":
    """The 21 complexation equilibria of Herrmann et al. (2023), Table S1."""
    return Database([Species(n, dict(nu), float(k), _charge(n)) for n, nu, k in HERRMANN_2023])


@dataclass
class Activity:
    """Activity coefficients.

    kind = "phreeqc": PHREEQC's model at 25 C (A, B of the Debye-Hueckel equation; WATEQ for species with
    -gamma a b, Davies otherwise, 0.1 I for neutral species), with a_w = 1 - 0.017 sum m.
    kind = "ideal": all activity coefficients 1 and a_w = 1 (concentrations), as in Herrmann et al.,
    Adv. Energy Mater. (2023), Supporting Information, Table S1 and Eqs. S1-S4.
    """
    A: float = 0.5100247894     # PHREEQC at 25 C (computed from water properties; checked against PHREEQC output)
    B: float = 0.3284906340
    neutral_b: float = 0.1
    kind: str = "phreeqc"

    def log10_gamma(self, sp: Species, I):
        if self.kind == "ideal":
            return np.zeros_like(np.asarray(I, float))
        sI = np.sqrt(I)
        if sp.z == 0:
            return self.neutral_b * I
        if sp.gamma is not None:
            a, b = sp.gamma
            return -self.A * sp.z ** 2 * sI / (1.0 + self.B * a * sI) + b * I
        return -self.A * sp.z ** 2 * (sI / (1.0 + sI) - 0.3 * I)


class Speciation:
    """Vectorized solver: pH and species distribution for arrays of total Zn, Mn, S (mol/kgw)."""

    def __init__(self, db: Database, activity: Activity | None = None):
        self.act = activity or Activity()
        masters = [Species(m, {m: 1.0}, 0.0, CHARGE_OF[m], db.master_gamma.get(m)) for m in MASTERS]
        self.species = masters + [s for s in db.species if s.name not in MASTERS]
        ns, nm = len(self.species), len(MASTERS)
        self.nu = np.zeros((ns, nm))
        self.nu_w = np.zeros(ns)
        self.logk = np.array([s.logk for s in self.species])
        self.z = np.array([s.z for s in self.species], dtype=float)
        for r, s in enumerate(self.species):
            for m, c in s.nu.items():
                if m == "H2O":
                    self.nu_w[r] = c
                else:
                    self.nu[r, MASTERS.index(m)] = c
        # element content of each species: Zn, Mn, S
        self.e_zn = self.nu[:, 1]
        self.e_mn = self.nu[:, 2]
        self.e_s = self.nu[:, 3]

    def _log_gamma(self, I):
        return np.stack([self.act.log10_gamma(s, I) for s in self.species], axis=-1)

    def initial_water_charge(self) -> float:
        """Charge imbalance (eq/kgw) of PHREEQC's starting solution: pure water at pH 7 exactly.

        add_solution_simple() starts from SOLUTION 1 (pure water, pH 7, not charge-balanced) and adds the
        salts as a REACTION; PHREEQC then conserves that small imbalance. Solving for it here reproduces
        PHREEQC's pH at high dilution (about 0.0015 pH units at 1e-8 mol/kgw).
        """
        iOH = next(r for r, sp in enumerate(self.species) if sp.name == "OH-")
        I = 1e-7
        for _ in range(50):
            gH = 10.0 ** self.act.log10_gamma(self.species[0], I)
            gOH = 10.0 ** self.act.log10_gamma(self.species[iOH], I)
            mH = 1e-7 / gH
            mOH = 10.0 ** (self.logk[iOH] + 7.0) / gOH
            I = 0.5 * (mH + mOH)
        return mH - mOH

    def solve(self, zn, mn, s=None, h_excess=0.0, x0=None, tol=1e-12, max_iter=200, charge0=None):
        """Totals (arrays, mol/kgw) -> dict(pH, I, aw, molality (..., n_species), x = log10 master activities).

        s defaults to zn + mn (the sulfate salts). h_excess (mol/kgw) is acid added as H+ with its sulfate
        counted in s. charge0 is the conserved charge of the starting water (PHREEQC's convention); by
        default the pure-water-at-pH-7 value, so results match PHREEQC's add_solution_simple.
        """
        zn = np.atleast_1d(np.asarray(zn, float))
        mn = np.broadcast_to(np.asarray(mn, float), zn.shape)
        s = zn + mn if s is None else np.broadcast_to(np.asarray(s, float), zn.shape)
        q0 = self.initial_water_charge() if charge0 is None else charge0
        q = q0 + np.broadcast_to(np.asarray(h_excess, float), zn.shape)
        n = zn.shape
        if x0 is None:
            x = np.stack([np.full(n, -6.0), np.log10(np.maximum(zn, 1e-300)) - 0.3,
                          np.log10(np.maximum(mn, 1e-300)) - 0.3,
                          np.log10(np.maximum(s, 1e-300)) - 0.3], axis=-1)
        else:
            x = np.array(x0, float)
        I = np.minimum(2.0 * (zn + mn + s), 0.1) + 1e-7
        aw = np.ones(n)
        E = np.stack([self.e_zn, self.e_mn, self.e_s, self.z], axis=0)          # (4, ns)
        tot = np.stack([zn, mn, s], axis=-1)
        absz = np.abs(self.z)

        def residual(x, lg, aw):
            m = 10.0 ** np.clip(self.logk + x @ self.nu.T + self.nu_w * np.log10(aw)[..., None] - lg, -300.0, 300.0)
            bal = m @ E[:3].T
            qs = m @ self.z - q
            zsum = m @ absz
            F = np.concatenate([np.log10(bal) - np.log10(tot), (qs / zsum)[..., None]], axis=-1)
            return F, m, bal, qs, zsum

        def inner(x, lg, aw, active):
            """Newton with backtracking at fixed gamma and a_w, on the nodes in `active`."""
            for _ in range(60):
                F, m, bal, qs, zsum = residual(x, lg, aw)
                conv = np.max(np.abs(F), axis=-1) < 1e-13
                act = active & ~conv
                if not act.any():
                    break
                dm = m[..., :, None] * self.nu * LN10
                Jb = np.einsum("es,...sk->...ek", E[:3], dm) / (bal[..., :, None] * LN10)
                dq = np.einsum("s,...sk->...k", self.z, dm)
                dz = np.einsum("s,...sk->...k", absz, dm)
                Jq = (dq * zsum[..., None] - qs[..., None] * dz) / (zsum ** 2)[..., None]
                J = np.concatenate([Jb, Jq[..., None, :]], axis=-2)
                J = np.where(np.isfinite(J), J, 0.0) + 1e-14 * np.eye(4)
                dx = np.linalg.solve(J, -np.where(np.isfinite(F), F, 0.0)[..., None])[..., 0]
                dx = np.where(act[..., None] & np.isfinite(dx), dx, 0.0)
                big = np.max(np.abs(dx), axis=-1, keepdims=True)
                dx = dx * np.minimum(1.0, 1.0 / np.maximum(big, 1e-300))   # scale, don't clip: keep the direction
                f0 = np.linalg.norm(F, axis=-1)
                lam = np.ones(n)
                for _ in range(30):
                    F1 = residual(x + lam[..., None] * dx, lg, aw)[0]
                    worse = act & ~(np.linalg.norm(F1, axis=-1) < f0)
                    if not worse.any():
                        break
                    lam = np.where(worse, 0.5 * lam, lam)
                x = x + lam[..., None] * dx
            return x

        def g_of(lnI, x):
            """ln I_new - ln I for given ln I (inner solve warm-started from x)."""
            I_ = np.exp(lnI)
            lg = self._log_gamma(I_)
            x = inner(x, lg, aw_of(I_, x, lg), np.ones(n, bool))
            m = residual(x, lg, aw_of(I_, x, lg))[1]
            return np.log(np.maximum(0.5 * (m @ (self.z ** 2)), 1e-300)) - lnI, x

        def aw_of(I_, x, lg):
            return aw                                            # a_w is iterated in the outer loop below

        lo = np.log(np.maximum(1e-9, 1e-3 * (zn + mn)))
        hi = np.log(0.5 * (36.0 * (zn + mn + s)) + 1.0)
        it = 0
        for outer in range(8):                                   # a_w fixed point around the I root
            # scan for the first sign change of g from low I (the branch PHREEQC reaches from pure water)
            grid = np.linspace(0.0, 1.0, 25)
            g_prev, x_scan = g_of(lo, x)
            a_, b_ = lo.copy(), hi.copy()
            found = np.zeros(n, bool)
            x_lo = x_scan.copy()
            for t in grid[1:]:
                l_t = lo + t * (hi - lo)
                g_t, x_t = g_of(l_t, x_scan)
                cross = ~found & (g_prev > 0) & (g_t <= 0)
                a_ = np.where(cross, lo + (t - grid[1]) * (hi - lo) if False else a_, a_)
                a_ = np.where(cross, l_t - (hi - lo) / 24.0, a_)
                b_ = np.where(cross, l_t, b_)
                x_lo = np.where(cross[..., None], x_scan, x_lo)
                found |= cross
                g_prev, x_scan = g_t, x_t
            # bisection on [a_, b_]
            xb = x_lo
            for it in range(60):
                mid = 0.5 * (a_ + b_)
                g_m, xb = g_of(mid, xb)
                a_ = np.where(g_m > 0, mid, a_)
                b_ = np.where(g_m > 0, b_, mid)
                if np.max(b_ - a_) < 1e-13:
                    break
            I = np.exp(0.5 * (a_ + b_))
            lg = self._log_gamma(I)
            x = inner(xb, lg, aw, np.ones(n, bool))
            m = residual(x, lg, aw)[1]
            aw_new = np.clip(1.0 - 0.017 * m.sum(axis=-1), 0.05, 1.0) if self.act.kind != "ideal" else np.ones(n)
            if np.max(np.abs(aw_new - aw)) < tol:
                aw = aw_new
                break
            aw = aw_new
        done = found | (np.abs(np.log(np.maximum(0.5 * (m @ (self.z ** 2)), 1e-300)) - np.log(I)) < 1e-9)
        lg = self._log_gamma(I)
        log_a = self.logk + x @ self.nu.T + self.nu_w * np.log10(aw)[..., None]
        m = 10.0 ** (log_a - lg)
        return dict(pH=-x[..., 0], I=I, aw=aw, molality=m, x=x, iterations=it + 1, converged=done,
                    names=[sp.name for sp in self.species])

    def solve_warm(self, zn, mn, x0, I0, aw0=None, s=None, h_excess=0.0, max_iter=30, tol=1e-12, charge0=None):
        """Warm-started solve (the previous time step's x, I and a_w): Newton on the master activities and a
        secant iteration on ln I. Returns the same dict as solve(); falls back to solve() where it fails."""
        zn = np.atleast_1d(np.asarray(zn, float))
        mn = np.broadcast_to(np.asarray(mn, float), zn.shape)
        s = zn + mn if s is None else np.broadcast_to(np.asarray(s, float), zn.shape)
        q0 = self.initial_water_charge() if charge0 is None else charge0
        q = q0 + np.broadcast_to(np.asarray(h_excess, float), zn.shape)
        E = np.stack([self.e_zn, self.e_mn, self.e_s, self.z], axis=0)
        tot = np.stack([zn, mn, s], axis=-1)
        absz = np.abs(self.z)
        x = np.array(x0, float)
        aw = np.ones(zn.shape) if aw0 is None else np.array(aw0, float)

        def inner(x, lg, aw):
            for _ in range(20):
                m = 10.0 ** np.clip(self.logk + x @ self.nu.T + self.nu_w * np.log10(aw)[..., None] - lg, -300, 300)
                bal = m @ E[:3].T
                qs = m @ self.z - q
                zsum = m @ absz
                F = np.concatenate([np.log10(bal) - np.log10(tot), (qs / zsum)[..., None]], axis=-1)
                if np.max(np.abs(F)) < 1e-13:
                    return x, m
                dm = m[..., :, None] * self.nu * LN10
                Jb = np.einsum("es,...sk->...ek", E[:3], dm) / (bal[..., :, None] * LN10)
                dq = np.einsum("s,...sk->...k", self.z, dm)
                dz = np.einsum("s,...sk->...k", absz, dm)
                Jq = (dq * zsum[..., None] - qs[..., None] * dz) / (zsum ** 2)[..., None]
                J = np.concatenate([Jb, Jq[..., None, :]], axis=-2)
                dx = np.linalg.solve(J, -F[..., None])[..., 0]
                big = np.max(np.abs(dx), axis=-1, keepdims=True)
                x = x + dx * np.minimum(1.0, 1.0 / np.maximum(big, 1e-300))
            return x, m

        l1 = np.log(np.asarray(I0, float))
        x, m = inner(x, self._log_gamma(np.exp(l1)), aw)
        g1 = np.log(0.5 * (m @ self.z ** 2)) - l1
        l2 = l1 + g1
        it = 0
        for it in range(max_iter):
            x, m = inner(x, self._log_gamma(np.exp(l2)), aw)
            g2 = np.log(0.5 * (m @ self.z ** 2)) - l2
            aw = np.clip(1.0 - 0.017 * m.sum(axis=-1), 0.05, 1.0)
            if np.max(np.abs(g2)) < tol:
                break
            d = np.where(np.abs(g2 - g1) > 1e-300, (l2 - l1) / (g2 - g1), 0.0)
            l_new = np.clip(l2 - g2 * d, l2 - 1.0, l2 + 1.0)
            l1, g1, l2 = l2, g2, l_new
        I = np.exp(l2)
        ok = np.abs(g2) < 1e-9
        out = dict(pH=-x[..., 0], I=I, aw=aw, molality=m, x=x, iterations=it + 1, converged=ok,
                   names=[sp.name for sp in self.species])
        if not ok.all():                                            # cold-start the failures
            bad = ~ok
            r = self.solve(zn[bad], mn[bad], s=s[bad], charge0=charge0,
                           h_excess=np.broadcast_to(np.asarray(h_excess, float), zn.shape)[bad])
            for key in ("pH", "I", "aw", "x", "converged"):
                out[key] = np.array(out[key]); out[key][bad] = r[key]
            out["molality"] = np.array(out["molality"]); out["molality"][bad] = r["molality"]
        return out

