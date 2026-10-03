"""Corrected model: finite volumes, implicit chemistry, Newton (plan: Phase 2 spec and option matrix).

Unknowns per cell (column):
    0 phi1    solid potential [V] (cathode; 0 elsewhere)
    1 phi2    electrolyte potential [V]
    2 Zn_T, 3 Mn_T, 4 S_T, 5 H_T   totals in the electrolyte [mol/cm3]; H_T may be negative
    6 n_MnO2  pristine MnO2 (R1) [mol/cm3 of electrode]
    7 n_ZMO   Zn_z MnO2 (R2) [mol Mn/cm3]
    8 s       insertion log-odds ln(theta/(1-theta)) (R3)
    9 n_ZHS   zinc hydroxide sulfate [mol/cm3]
   10 n_ZnO, 11 n_Zn(OH)2  optional precipitates [mol/cm3]

Cells run from the Zn anode (x = 0) through the probe region, the separator and the cathode to the
current collector. The anode is a boundary with Butler-Volmer kinetics; the anode metal is the
potential reference (0 V), so the cell voltage is phi1 at the collector face.

The residual is linear(x) + sources(x) + transport(x). The linear part (storage, electroneutrality,
trivial rows) and the transport have analytic Jacobians; only the local, nonlinear sources (reactions,
speciation, anode) are differentiated by finite differences, with steps relative to each value.
Differencing the linear terms loses every digit when a species is depleted: a step small enough for its
logarithm is below the round-off of the sums it appears in.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import bandsolver

from .driver import SolverFailure
from .eqchem import LN10, Equilibria, read_logk_file
from .params import Params

N = 12
P1, P2, ZN, MN, SO, H, MO, ZM, TH, ZH, ZO, ZX = range(12)
SPECIES = (ZN, MN, SO, H)
CHARGE = {ZN: 2.0, MN: 2.0, SO: -2.0, H: 1.0}
SOLIDS = (MO, ZM, ZH, ZO, ZX)             # amounts that set the porosity
PRECIP = (ZH, ZO, ZX)                     # precipitates: Zn, SO4 consumed and H+ released per unit
PROBE, SEP, CATH = 0, 1, 2
TOT_COLS = (H, ZN, MN, SO)                # column of each speciation total (masters H+, Zn2+, Mn2+, SO4 2-)
W_SMOOTH = 0.01      # width of the smoothed onset of seed precipitation (in w = c_H,eq / c_H)
S_MAX = 200.0         # |log-odds| bound of the insertion fraction (theta, 1 - theta >= 1e-87, representable)
S_STEP = 20.0         # largest step in s per Newton iteration where theta <= 1/2
TH_FLOOR = 1e-8       # theta (or 1 - theta) falls by at most this factor in one Newton iteration
RES_TOL = 1e-6        # residual at convergence, relative to the 1 mA/g equivalent of each row (res_scale)
TRACE = 1.0e-9        # [mol/L] initial concentration of a salt given as 0
FB_REF = 1.0e-6       # amount scale [mol/cm3] of the equilibrium-ZHS complementarity condition


@dataclass(frozen=True)
class Mesh:
    x: np.ndarray        # cell centres [cm]
    dx: np.ndarray       # cell widths
    area: np.ndarray     # cross-section [cm2]
    region: np.ndarray   # PROBE, SEP, CATH
    c0: int              # first cathode cell

    @property
    def n(self) -> int:
        return len(self.x)

    @property
    def cath(self) -> np.ndarray:
        return self.region == CATH


def make_mesh(p: Params) -> Mesh:
    parts = []
    if p.L_probe > 0:
        parts.append((PROBE, p.L_probe, p.n_probe))
    parts.append((SEP, p.L_sep, p.n_sep))
    parts.append((CATH, p.L_cath, p.n_cath))
    x, dx, region = [], [], []
    start = 0.0
    for reg, L, n in parts:
        h = L / n
        x += list(start + h * (np.arange(n) + 0.5))
        dx += [h] * n
        region += [reg] * n
        start += L
    region = np.array(region)
    return Mesh(np.array(x), np.array(dx), np.full(len(x), p.A_cell), region, int(np.argmax(region == CATH)))


def bernoulli(a):
    """B(a) = a / (exp(a) - 1), with B(0) = 1."""
    a = np.asarray(a, float)
    small = np.abs(a) < 1e-6
    safe = np.where(small, 1.0, a)
    return np.where(small, 1.0 - a / 2.0 + a * a / 12.0, safe / np.expm1(safe))


def dbernoulli(a):
    """dB/da = (exp(a) - 1 - a exp(a)) / (exp(a) - 1)**2, with dB/da(0) = -1/2."""
    a = np.asarray(a, float)
    small = np.abs(a) < 1e-4
    safe = np.where(small, 1.0, a)
    em = np.expm1(safe)
    big = (em - safe * (em + 1.0)) / em ** 2
    return np.where(small, -0.5 + a / 6.0, big)


def _load_ocp():
    from .tables import r3_ocp
    t, c = r3_ocp()
    return np.asarray(t, float), np.asarray(c, float)


def logit(t):
    return np.log(t / (1.0 - t))


def theta_pair(s):
    """(theta, 1 - theta) = (1/(1+exp(-s)), 1/(1+exp(s))), each to full relative precision."""
    s = np.asarray(s, float)
    e = np.exp(-np.abs(s))
    big, small = 1.0 / (1.0 + e), e / (1.0 + e)
    return np.where(s >= 0, big, small), np.where(s >= 0, small, big)


def sigmoid(s):
    return theta_pair(s)[0]


def pos(z, width=1.0e-2):
    """A C1 max(z, 0) that is exactly zero for z <= 0: z^2 / (z + width) for z > 0.

    Used for one-directional terms (growth on foreign surfaces, irreversible dissolution), so that every
    rate vanishes at equilibrium (detailed balance: no current circulates in a cell at rest)."""
    z = np.asarray(z, float)
    zp = np.maximum(z, 0.0)
    return zp * zp / (zp + width)


class Model:
    def __init__(self, p: Params):
        self.p = p
        self.mesh = make_mesh(p)
        over = read_logk_file(p.logk_file) if p.logk_file else None
        self.eq = Equilibria(over, database=p.equilibria_db)
        self.f = p.F / (p.R * p.T)
        m = self.mesh
        self.cath = m.cath
        # molar volumes [cm3/mol] and the fixed solid volume of the cathode (inert + host)
        self.V = np.zeros(N)
        self.V[MO], self.V[ZM], self.V[ZH] = p.M_MnO2 / p.rho_MnO2, p.M_ZMO / p.rho_ZMO, p.M_ZHS / p.rho_ZHS
        self.V[ZO], self.V[ZX] = p.M_ZnO / p.rho_ZnO, p.M_ZnOH2 / p.rho_ZnOH2
        self.V_MO, self.V_ZM, self.V_ZH = self.V[MO], self.V[ZM], self.V[ZH]
        self.eps_fixed = 1.0 - p.eps_cath - p.vf_MnO2 - p.vf_ZMO - p.vf_ZHS   # inert + host
        self.n_host = p.vf_host * p.rho_host / p.M_host                       # mol Mn of the host per cm3
        self.a_host = 3.0 * p.vf_host / p.r_host
        self.D = {ZN: p.D_Zn, MN: p.D_Mn, SO: p.D_SO4, H: p.D_H}
        # precipitates: (on, log K, Zn, SO4, H+ per unit, rate constant)
        self.precip = {ZH: (p.zhs == "kinetic", p.logK_ZHS, 4.0, 1.0, 6.0, p.k_ZHS),
                       ZO: (p.ZnO_on, p.logK_ZnO, 1.0, 0.0, 2.0, p.k_ZnO),
                       ZX: (p.ZnOH2_on, p.logK_ZnOH2, 1.0, 0.0, 2.0, p.k_ZnOH2)}
        self.transported = (ZN, MN, SO) if p.species == "no_H" else SPECIES
        if p.r3_ocp != "nernst":
            self.ocp_t, self.ocp_c = _load_ocp()
        if p.transport == "quasi":
            dmap = {"H+": p.D_H, "Zn+2": p.D_Zn, "Mn+2": p.D_Mn, "SO4-2": p.D_SO4, "OH-": p.D_OH, "HSO4-": p.D_HSO4}
            self.D_sp = np.array([dmap.get(nm, p.D_complex) for nm in self.eq.names])
        self.H0 = None
        self._spec_x = None                                                     # warm start for the speciation
        self._cache = (None, None)
        self.typ = np.array([1.0, 1.0, 1e-3, 1e-3, 1e-3, 1e-6, 1e-3, 1e-3, 1.0, 1e-3, 1e-3, 1e-3])
        # physical residual scales: the 1 mA/g equivalent of each row (current [A], species [mol/s],
        # solids and insertion [mol/cm3/s of cathode]); electroneutrality in mol/L
        i_ref = 1.0e-3 * p.mass
        v_cath = float(np.sum((m.area * m.dx)[self.cath]))
        self.res_scale = np.full(N, i_ref / (p.F * v_cath))
        self.res_scale[P1] = i_ref
        self.res_scale[P2] = 1.0
        for k in SPECIES:
            self.res_scale[k] = i_ref / p.F
        if p.zhs == "equilibrium":
            self.res_scale[ZH] = 1.0                    # the complementarity row is dimensionless

    # ------------------------------------------------------------------ state
    def initial_state(self) -> np.ndarray:
        p, m = self.p, self.mesh
        x = np.zeros((m.n, N))
        # an absent salt starts as a trace (TRACE mol/L): the equilibrium potentials need ln c
        c_zn, c_mn = max(p.c_ZnSO4, TRACE), max(p.c_MnSO4, TRACE)
        x[:, ZN] = c_zn * 1e-3
        x[:, MN] = c_mn * 1e-3
        x[:, SO] = (c_zn + c_mn + p.c_H2SO4) * 1e-3
        x[:, H] = 2.0 * p.c_H2SO4 * 1e-3
        self.H0 = x[:, H].copy()
        c = self.cath
        x[c, MO] = p.vf_MnO2 / self.V[MO]
        x[c, ZM] = p.vf_ZMO / self.V[ZM]
        x[c, ZH] = p.vf_ZHS / self.V[ZH]
        x[:, TH] = logit(p.theta0)
        free = self.free(x)
        u_zn = 0.5 / self.f * np.log(free[0, 1])
        x[:, P2] = -u_zn                                             # anode at equilibrium
        cb = self._basis(x, free)
        x[c, P1] = x[c, P2] + self.u3(x[c, TH], cb[c, 1])            # insertion at equilibrium (a guess)
        # settle the potentials at open circuit: a step so short (1e-4 s) that the compositions barely change
        # the reactions' mixed potential (net current zero in every cell) and the anode equilibrium
        return self.newton_step(x, 1.0e-4, 0.0)

    def _lx(self, x, warm=True):
        """log10 free concentrations [mol/L] (n, 4) of H+, Zn2+, Mn2+, SO4 2-, from the totals."""
        T = np.stack([x[:, k] for k in TOT_COLS], axis=-1) * 1e3
        key = T.tobytes()
        if self._cache[0] == key:
            return self._cache[1], T
        x0 = self._spec_x if (self._spec_x is not None and self._spec_x.shape == T.shape) else None
        try:
            lx = self.eq.solve(T, x0=x0)
        except RuntimeError:
            try:
                lx = self.eq.solve(T)
            except RuntimeError as e:                   # a Newton trial state the chemistry cannot solve
                raise SolverFailure(str(e)) from e
        if warm:
            self._spec_x = lx
            self._cache = (key, lx)
        return lx, T

    def free(self, x, warm=True) -> np.ndarray:
        """Free concentrations [mol/L] (n, 4): H+, Zn2+, Mn2+, SO4 2-."""
        return 10.0 ** self._lx(x, warm)[0]

    def _basis(self, x, free):
        """Concentrations [mol/L] (n, 4) of H+, Zn, Mn, SO4 in the Nernst and rate terms (options basis, ph_mode)."""
        p = self.p
        c = free.copy()
        if p.basis == "totals":
            c[:, 1], c[:, 2], c[:, 3] = x[:, ZN] * 1e3, x[:, MN] * 1e3, x[:, SO] * 1e3
        if p.ph_mode == "zhs_equilibrium":
            c[:, 0] = np.exp((4.0 * np.log(c[:, 1]) + np.log(c[:, 3]) - p.logK_ZHS * LN10) / 6.0)
        elif p.ph_mode == "fixed":
            c[:, 0] = 10.0 ** (-p.pH_fixed)
        elif p.ph_mode == "spline":
            from . import phspline
            c[:, 0] = 10.0 ** (-phspline.ph(x[:, ZN] * 1e3, x[:, MN] * 1e3))
        return c

    # ------------------------------------------------------------------ constitutive
    def porosity(self, x):
        p, m = self.p, self.mesh
        eps = np.where(m.region == PROBE, p.eps_probe, p.eps_sep).astype(float)
        c = self.cath
        eps[c] = 1.0 - self.eps_fixed - sum(self.V[k] * x[c, k] for k in SOLIDS)
        return eps

    def tortuosity(self, eps):
        p, m = self.p, self.mesh
        if np.any(eps <= 0.0):                          # a trial state with the pores overfilled by solids
            raise SolverFailure("non-positive porosity")
        tau = np.full(m.n, p.tau_probe)
        s = m.region == SEP
        tau[s] = p.tau_factor_sep * eps[s] ** p.bruggeman_sep
        c = self.cath
        tau[c] = p.tau_factor_cath * eps[c] ** p.bruggeman_cath
        return tau

    def u3(self, s, c_zn):
        """R3 equilibrium potential from the insertion log-odds s and Zn2+ [mol/L]."""
        p = self.p
        if p.r3_ocp == "nernst":
            return p.U3_nernst + 0.5 / self.f * (np.log(c_zn) - s)
        theta, om = theta_pair(s)
        end = 0.0
        if p.r3_ocp == "spline_nernst":
            # Nernstian ends: -(RT/2F) ln(theta/(1-theta)) beyond about r3_end_width of either end, added
            # smoothly (a soft min of ln(theta/eps) and 0, in log space, width 1/2) so the middle of the
            # empirical curve is unchanged. A finite OCP at the ends lets the host fill or empty exactly in
            # finite time and then never react again (the charge failure); this keeps both directions open.
            eps, w = p.r3_end_width, 0.5
            lt = np.minimum(s, 0.0) - np.log1p(np.exp(-np.abs(s)))         # ln(theta), accurate for s << 0
            lo = np.minimum(-s, 0.0) - np.log1p(np.exp(-np.abs(s)))        # ln(1 - theta)
            soft = lambda a: -w * np.logaddexp(0.0, -a / w)              # smooth min(a, 0)
            end = -(0.5 / self.f) * (soft(lt - math.log(eps)) - soft(lo - math.log(eps)))
        t, cf = self.ocp_t, self.ocp_c
        k = np.clip(np.searchsorted(t, theta, side="left") - 1, 0, len(t) - 2)
        d = theta - t[k]
        vs = ((cf[k, 0] * d + cf[k, 1]) * d + cf[k, 2]) * d + cf[k, 3]
        return (p.V_at_zmin - p.V_at_zmax) * vs + p.V_at_zmax + 0.5 / self.f * np.log(c_zn / p.c_ref3) + end

    def area(self, x, k):
        p = self.p
        r = {MO: p.r_MnO2, ZM: p.r_ZMO, ZH: p.r_ZHS, ZO: p.r_ZHS, ZX: p.r_ZHS}[k]
        return 3.0 * self.V[k] * np.maximum(x[self.cath, k], 0.0) / r

    def precipitation(self, x, free, k):
        """Kinetic precipitation rate of solid k [mol/cm3/s] in the cathode (positive = forms)."""
        p = self.p
        on, logk, nzn, ns, nh, kr = self.precip[k]
        c = self.cath
        if not on:
            return np.zeros(int(c.sum()))
        cH, cZn, cS = free[c, 0], free[c, 1], free[c, 3]
        lq = nzn * np.log(cZn) + (ns * np.log(cS) if ns else 0.0) - logk * LN10
        w = np.exp(lq / nh) / cH                                   # = c_H,eq / c_H
        w_nuc = p.zhs_nucleation ** (nzn / nh)                     # Zn supersaturation -> w
        return kr * (self.area(x, k) * (w - 1.0) + p.a_seed_ZHS * pos(w - w_nuc, W_SMOOTH))

    def zhs_w(self, free):
        p = self.p
        c = self.cath
        return np.exp((4.0 * np.log(free[c, 1]) + np.log(free[c, 3]) - p.logK_ZHS * LN10) / 6.0) / free[c, 0]

    def reactions(self, x, free):
        """Reaction currents (A/cm3, positive anodic) and the kinetic ZHS rate (mol/cm3/s), per cathode cell."""
        p, f = self.p, self.f
        c = self.cath
        cb = self._basis(x, free)
        cH, cZn, cMn = cb[c, 0], cb[c, 1], cb[c, 2]
        p1, p2 = x[c, P1], x[c, P2]
        a_MO, a_ZM, a_ZH = self.area(x, MO), self.area(x, ZM), self.area(x, ZH)
        # R1: irreversible dissolution, only below its equilibrium potential
        u1 = p.U1 - 0.5 / f * (np.log(cMn) - 4.0 * np.log(cH))
        e1 = p1 - p2 - u1
        bv1 = np.exp(p.alpha1 * 2.0 * f * e1) - np.exp(-(1.0 - p.alpha1) * 2.0 * f * e1)
        i1 = -a_MO * p.F * p.k1 * pos(-bv1) * p.R1_on
        # R2: oxidation on ZMO + ZHS + seed, reduction on ZMO
        n2 = 2.0 - 2.0 * p.z_ZMO
        u2 = p.U2 - (1.0 / (n2 * f)) * (p.z_ZMO * np.log(cZn) + np.log(cMn) - 4.0 * np.log(cH))
        e2 = p1 - p2 - u2
        bv2 = np.exp(p.alpha2 * n2 * f * e2) - np.exp(-(1.0 - p.alpha2) * n2 * f * e2)
        if p.zhs == "lumped":
            # lumped ZHS (2024 paper): deposition consumes ZHS, so it proceeds on the ZHS present (the
            # original's charge-area switch, M-4, made smooth); dissolution on ZMO
            i2 = p.F * p.k2 * (a_ZH * pos(bv2) - a_ZM * pos(-bv2)) * p.R2_on
        else:
            # both directions on ZMO; deposition also grows on ZHS and the seed area
            i2 = p.F * p.k2 * (a_ZM * bv2 + (a_ZH + p.a_seed_R2) * pos(bv2)) * p.R2_on
        # R3: insertion
        th, om = theta_pair(x[c, TH])
        e3 = p1 - p2 - self.u3(x[c, TH], cZn)
        i0 = p.F * p.k3 * np.sqrt(cZn * th * om)
        i3 = self.a_host * i0 * (np.exp(p.alpha3 * 2.0 * f * e3) - np.exp(-(1.0 - p.alpha3) * 2.0 * f * e3)) * p.R3_on
        return i1, i2, i3, self.precipitation(x, free, ZH)

    def anode_current(self, x, free):
        """Zn anode current density [A/cm2] (positive = Zn dissolves) on the first cell's face."""
        p, f = self.p, self.f
        czn = free[0, 1] if p.basis == "free" else x[0, ZN] * 1e3
        eta = 0.0 - x[0, P2] - 0.5 / f * np.log(czn)
        return p.F * p.k_an * np.sqrt(czn) * (np.exp(p.alpha_an * 2.0 * f * eta)
                                              - np.exp(-(1.0 - p.alpha_an) * 2.0 * f * eta))

    @staticmethod
    def d_theta(s, s_old):
        """theta(s) - theta(s_old), from the complements when theta is near 1 (no cancellation)."""
        t, o = theta_pair(s)
        t0, o0 = theta_pair(s_old)
        return np.where((s > 0) & (s_old > 0), o0 - o, t - t0)

    def n_ins(self, s):
        p = self.p
        return self.n_host * (p.zmin + sigmoid(s) * (p.zmax - p.zmin))

    # ------------------------------------------------------------------ residual
    def linear(self, x, old, dt, I):
        """Linear terms and their Jacobian (n, N), (n, N, N): storage, electroneutrality, trivial rows,
        solid storage, the collector current, and (zhs = equilibrium) the ZHS sinks."""
        p, m = self.p, self.mesh
        n = m.n
        R = np.zeros((n, N))
        J = np.zeros((n, N, N))
        vol = m.area * m.dx
        eps, eps_old = self.porosity(x), self.porosity(old)
        c, nc = self.cath, ~self.cath
        for k in self.transported:
            R[:, k] = vol * (eps * x[:, k] - eps_old * old[:, k]) / dt
            J[:, k, k] = vol * eps / dt
            for s_ in SOLIDS:
                J[c, k, s_] = vol[c] * x[c, k] * (-self.V[s_]) / dt
        if p.species == "no_H":
            R[:, H] = x[:, H] - self.H0
            J[:, H, H] = 1.0
        R[:, P2] = (2.0 * x[:, ZN] + 2.0 * x[:, MN] + x[:, H] - 2.0 * x[:, SO]) * 1e3
        J[:, P2, ZN], J[:, P2, MN], J[:, P2, H], J[:, P2, SO] = 2e3, 2e3, 1e3, -2e3
        for k in (P1,) + SOLIDS:
            R[nc, k] = x[nc, k]
            J[nc, k, k] = 1.0
        R[nc, TH] = x[nc, TH] - old[nc, TH]
        J[nc, TH, TH] = 1.0
        for k in SOLIDS:
            if k == ZH and p.zhs == "equilibrium":
                continue                                  # its row is the complementarity condition (sources)
            R[c, k] = (x[c, k] - old[c, k]) / dt
            J[c, k, k] = 1.0 / dt
        if p.zhs == "equilibrium":                        # ZHS formed this step is taken from the solution
            vc = vol[c]
            for k, nu in ((ZN, 4.0), (SO, 1.0), (H, -6.0)):
                if k in self.transported:
                    R[c, k] += nu * vc * (x[c, ZH] - old[c, ZH]) / dt
                    J[c, k, ZH] += nu * vc / dt
        cap = self.n_host * (p.zmax - p.zmin)
        th, om = theta_pair(x[c, TH])
        R[c, TH] = cap * self.d_theta(x[c, TH], old[c, TH]) / dt
        J[c, TH, TH] = cap * th * om / dt
        R[-1, P1] += I                                   # the solid current I leaves the last cell
        return R, J

    def sources(self, x, warm=True):
        """Local nonlinear terms (n, N): reactions and precipitation in the cathode, the anode flux into cell 0."""
        p, m = self.p, self.mesh
        F = p.F
        R = np.zeros((m.n, N))
        free = self.free(x, warm=warm)
        c = self.cath
        vc = (m.area * m.dx)[c]
        i1, i2, i3, r = self.reactions(x, free)
        n2 = 2.0 - 2.0 * p.z_ZMO
        xi1, xi2, xi3 = -i1 / (2.0 * F), -i2 / (n2 * F), -i3 / (2.0 * F)
        if p.zhs == "lumped":                            # the 4 H+ per Mn come from ZHS forming (2024 paper)
            r = (2.0 / 3.0) * (xi1 + xi2)
        elif p.zhs in ("off", "equilibrium"):
            r = np.zeros_like(xi1)
        rates = {ZH: r, ZO: self.precipitation(x, free, ZO), ZX: self.precipitation(x, free, ZX)}
        dzn = p.z_ZMO * xi2 - xi3
        dso = 0.0 * xi2
        dh = -4.0 * xi1 - 4.0 * xi2
        for k, rk in rates.items():
            _, _, nzn, ns, nh, _ = self.precip[k]
            dzn, dso, dh = dzn - nzn * rk, dso - ns * rk, dh + nh * rk
        R[c, ZN] = -vc * dzn
        R[c, MN] = -vc * (xi1 + xi2)
        R[c, SO] = -vc * dso
        if p.species == "with_H":
            R[c, H] = -vc * dh
        R[c, P1] = vc * (i1 + i2 + i3)
        R[c, MO] = xi1
        R[c, ZM] = xi2
        R[c, TH] = -xi3
        for k, rk in rates.items():
            R[c, k] = -rk
        if p.zhs == "equilibrium":                       # n_ZHS >= 0, w <= 1, n_ZHS (1 - w) = 0 (Fischer-Burmeister)
            a, b = x[c, ZH] / FB_REF, 1.0 - self.zhs_w(free)
            R[c, ZH] = a + b - np.sqrt(a * a + b * b + 1e-20)
        R[0, ZN] -= self.anode_current(x, free) * m.area[0] / (2.0 * F)
        return R

    def local(self, x, old, dt, I, warm=True):
        return self.linear(x, old, dt, I)[0] + self.sources(x, warm=warm)

    def _face_geometry(self, x):
        p, m = self.p, self.mesh
        eps = self.porosity(x)
        tau = self.tortuosity(eps)
        half = 0.5 * m.dx * tau / (eps * m.area)
        G = 1.0 / (half[:-1] + half[1:])                                  # geometric face conductance [cm]
        c = self.cath
        bexp = np.where(c, p.bruggeman_cath, 0.0)
        dhalf = np.where(c, (bexp - 1.0) * half / eps, 0.0)              # d half / d eps
        return eps, G, dhalf

    def _add_face(self, A, B, D, row, dL, dR):
        """Add the derivatives of an east-face outflow of `row` (n-1 faces) to the blocks."""
        n = self.mesh.n
        L, Rr = slice(0, n - 1), slice(1, n)
        B[L, row, :] += dL; D[L, row, :] += dR
        A[Rr, row, :] -= dL; B[Rr, row, :] -= dR

    def transport(self, x, jacobian=False):
        """Face fluxes between neighbouring cells (n, N): outflow east minus inflow west.
        With jacobian=True also returns the blocks (A, B, D) of their derivatives."""
        p, m = self.p, self.mesh
        n = m.n
        R = np.zeros((n, N))
        A = np.zeros((n, N, N)); B = np.zeros((n, N, N)); D = np.zeros((n, N, N))
        eps, G, dhalf = self._face_geometry(x)
        dphi = x[1:, P2] - x[:-1, P2]

        def solid_terms(dL, dR, flow, Gf):
            for s_ in SOLIDS:
                dL[:, s_] += -flow * Gf * dhalf[:-1] * (-self.V[s_])
                dR[:, s_] += -flow * Gf * dhalf[1:] * (-self.V[s_])

        if p.transport == "ions":
            for k in self.transported:
                zf = CHARGE[k] * self.f
                a = zf * dphi
                Bp, Bm = bernoulli(a), bernoulli(-a)
                flow = self.D[k] * G * (Bp * x[:-1, k] - Bm * x[1:, k])
                R[:-1, k] += flow
                R[1:, k] -= flow
                if jacobian:
                    dL = np.zeros((n - 1, N)); dR = np.zeros((n - 1, N))
                    dL[:, k] = self.D[k] * G * Bp
                    dR[:, k] = -self.D[k] * G * Bm
                    dfa = self.D[k] * G * (dbernoulli(a) * x[:-1, k] + dbernoulli(-a) * x[1:, k]) * zf
                    dR[:, P2], dL[:, P2] = dfa, -dfa
                    solid_terms(dL, dR, flow, G)
                    self._add_face(A, B, D, k, dL, dR)
        else:                                            # quasi-particles: each species' flux, summed per total
            lx, T = self._lx(x)
            cs = self.eq.species(lx) * 1e-3                                # species [mol/cm3] (n, ns)
            nu = self.eq.nu                                                # (ns, 4)
            z = self.eq.z
            a = (z[None, :] * self.f) * dphi[:, None]                      # (n-1, ns)
            Bp, Bm = bernoulli(a), bernoulli(-a)
            fl = self.D_sp * G[:, None] * (Bp * cs[:-1] - Bm * cs[1:])     # species flows (n-1, ns)
            if jacobian:
                S = self.eq.sensitivity(lx, T)                             # d log10 free / d T [L/mol]
                dc = cs[:, :, None] * LN10 * np.einsum("sm,nmk->nsk", nu, S) * 1e3   # d c_s / d total (n, ns, 4)
                dfa = self.D_sp * G[:, None] * (dbernoulli(a) * cs[:-1] + dbernoulli(-a) * cs[1:]) * z * self.f
            for q, k in enumerate(TOT_COLS):
                if k not in self.transported:
                    continue
                flow = fl @ nu[:, q]
                R[:-1, k] += flow
                R[1:, k] -= flow
                if jacobian:
                    dL = np.zeros((n - 1, N)); dR = np.zeros((n - 1, N))
                    wL = self.D_sp * G[:, None] * Bp * nu[:, q]            # d flow / d c_s(L)
                    wR = -self.D_sp * G[:, None] * Bm * nu[:, q]
                    for qq, kk in enumerate(TOT_COLS):
                        dL[:, kk] = np.sum(wL * dc[:-1, :, qq], axis=1)
                        dR[:, kk] = np.sum(wR * dc[1:, :, qq], axis=1)
                    d = dfa @ nu[:, q]
                    dR[:, P2], dL[:, P2] = d, -d
                    solid_terms(dL, dR, flow, G)
                    self._add_face(A, B, D, k, dL, dR)
        # solid current (cathode faces only)
        c = self.cath
        sig = np.where(c, p.sigma * (1.0 - eps), 1.0)
        hs = 0.5 * m.dx / (sig * m.area)
        both = c[:-1] & c[1:]
        Gs = np.where(both, 1.0 / (hs[:-1] + hs[1:]), 0.0)
        flow = Gs * (x[:-1, P1] - x[1:, P1])
        R[:-1, P1] += flow
        R[1:, P1] -= flow
        if jacobian:
            dhs = np.where(c, hs / np.where(c, 1.0 - eps, 1.0), 0.0)      # d hs / d eps
            dL = np.zeros((n - 1, N)); dR = np.zeros((n - 1, N))
            dL[:, P1], dR[:, P1] = Gs, -Gs
            for s_ in SOLIDS:
                dL[:, s_] = -flow * Gs * dhs[:-1] * (-self.V[s_])
                dR[:, s_] = -flow * Gs * dhs[1:] * (-self.V[s_])
            self._add_face(A, B, D, P1, dL, dR)
            return R, A, B, D
        return R

    def residual(self, x, old, dt, I):
        return self.local(x, old, dt, I) + self.transport(x)

    # ------------------------------------------------------------------ Jacobian and Newton
    def _steps(self, x):
        """Finite-difference steps for the sources: relative to each value (a depleted species, e.g. Mn2+
        near the collector on charge, needs a step far below its typical size); the potentials and s use
        their typical size as the floor. H_T can be zero; its floor is the free-proton level (~1e-8 mol/cm3,
        pH 5), so the step stays well above the speciation solver's tolerance."""
        floor = np.array([1.0, 1.0, 1e-15, 1e-15, 1e-15, 1e-8, 1e-15, 1e-15, 1.0, 1e-15, 1e-15, 1e-15])
        return 1e-7 * np.maximum(np.abs(x), floor)

    def blocks(self, x, old, dt, I):
        RL, B = self.linear(x, old, dt, I)
        S0 = self.sources(x)
        RT, A, BT, D = self.transport(x, jacobian=True)
        B = B + BT
        h = self._steps(x)
        spec_base, cache = self._spec_x, self._cache
        for k in range(N):
            xp = x.copy()
            xp[:, k] += h[:, k]
            self._spec_x = spec_base
            B[:, :, k] += (self.sources(xp, warm=False) - S0) / h[:, k][:, None]
        self._spec_x, self._cache = spec_base, cache
        return RL + S0 + RT, A, B, D

    def _damping(self, x, dx):
        """Largest lam <= 1 keeping the Zn, Mn and S totals positive (fraction to the boundary) and limiting
        potential changes to 0.2 V. The solids are projected per cell instead (newton_step)."""
        lam = 1.0
        for k in (ZN, MN, SO):
            xv, dv = x[:, k], dx[:, k]
            cross = (xv + dv < 0) & (xv > 0)                 # only steps that would cross zero are limited
            if cross.any():
                lam = min(lam, float(np.min(0.9 * xv[cross] / -dv[cross])))
        big = np.max(np.abs(dx[:, (P1, P2)]))
        if big > 0.2:
            lam = min(lam, 0.2 / big)
        return lam

    def newton_step(self, old, dt, I, start=None):
        p = self.p
        if self.H0 is None:
            self.H0 = old[:, H].copy()
        x = (old if start is None else start).copy()
        self._cache = (None, None)
        self._spec_x = None
        self.free(old)                                   # warm start from the step's start (not a failed attempt)
        for it in range(p.newton_max_iter):
            R, A, B, D = self.blocks(x, old, dt, I)
            if not np.all(np.isfinite(R)):
                raise SolverFailure("non-finite residual")
            # round-off floor of each residual: ~1000 eps times the size of its terms (|J| |x| over the band);
            # it matters for tiny steps (storage ~ 1/dt) and large flows
            ax = np.abs(x)
            mag = np.einsum("nij,nj->ni", np.abs(B), ax)
            mag[1:] += np.einsum("nij,nj->ni", np.abs(A[1:]), ax[:-1])
            mag[:-1] += np.einsum("nij,nj->ni", np.abs(D[:-1]), ax[1:])
            res_tol = np.maximum(RES_TOL * self.res_scale, 1e3 * np.finfo(float).eps * mag)
            # the insertion unknown is stored as s = logit(theta). Where theta > 1/2 it is solved for in
            # theta: near a full host d theta / d s vanishes, and with a Nernstian OCP the rate does not
            # depend on s at all, so the s column would be singular. Where theta <= 1/2 it is solved for in
            # s: with sqrt(theta) kinetics the derivative in theta diverges as theta -> 0, in s it does not.
            th, om = theta_pair(x[:, TH])
            in_theta = x[:, TH] > 0.0
            B[:, :, TH] = B[:, :, TH] / np.where(in_theta, th * om, 1.0)[:, None]
            # scale the columns by the typical size of each unknown and every row to unit maximum
            A, B, D = A * self.typ, B * self.typ, D * self.typ
            rs = np.maximum(np.max(np.abs(np.concatenate([A, B, D], axis=2)), axis=2), 1e-300)
            A, B, D, G = A / rs[..., None], B / rs[..., None], D / rs[..., None], -R / rs
            try:
                dx = bandsolver.solve(A, B, D, G, pivot="partial") * self.typ
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
                raise SolverFailure(str(e)) from e
            dth = dx[:, TH].copy()                       # a step in theta (see above)
            dx[:, TH] = 0.0
            lam = self._damping(x, dx)
            new = x + lam * dx
            # solids are local: a step that would make one negative takes it to a tenth of its value instead
            for k in SOLIDS:
                neg = new[:, k] < 0
                new[neg, k] = 0.1 * x[neg, k]
            # theta > 1/2: theta and 1 - theta updated separately (each keeps its full precision), neither
            # below TH_FLOOR of its current value in one iteration; s recomputed from them.
            # theta <= 1/2: a step in s, at most S_STEP.
            d = lam * dth
            th_new, om_new = th + d, om - d
            low, high = th_new < TH_FLOOR * th, om_new < TH_FLOOR * om
            th_new = np.where(low, TH_FLOOR * th, th_new)
            om_new = np.where(low, 1.0 - th_new, om_new)
            om_new = np.where(high, TH_FLOOR * om, om_new)
            th_new = np.where(high, 1.0 - om_new, th_new)
            s_theta = np.log(np.maximum(th_new, 1e-300)) - np.log(np.maximum(om_new, 1e-300))
            s_s = x[:, TH] + np.clip(d, -S_STEP, S_STEP)
            new[:, TH] = np.clip(np.where(in_theta, s_theta, s_s), -S_MAX, S_MAX)
            x = new
            sc = np.abs(dx) / self.typ
            sc[:, TH] = np.where(in_theta, np.abs(dth), np.abs(dth) * th * om)   # counted by its change in theta
            upd = float(np.max(sc))
            if not math.isfinite(upd):
                break
            # converged: a full step, a small update, and a small residual (the residual test catches a
            # variable held at a bound while its equation is unsatisfied)
            if lam == 1.0 and upd < p.newton_tol and np.all(np.abs(R) <= res_tol):
                self.free(x)                       # refresh the warm start at the converged state
                return x
        raise SolverFailure("Newton did not converge")

    # ------------------------------------------------------------------ outputs
    def voltage(self, x, I):
        p, m = self.p, self.mesh
        eps = self.porosity(x)[-1]
        return float(x[-1, P1] - I * 0.5 * m.dx[-1] / (p.sigma * (1.0 - eps) * m.area[-1]))

    def profiles(self, x):
        """Profiles across the cell for plotting: position [um] and per-cell quantities (NaN where a
        quantity does not exist, e.g. solids outside the cathode)."""
        c = self.cath
        free = self.free(x)
        nanc = lambda v: np.where(c, v, np.nan)
        th, _ = theta_pair(x[:, TH])
        return dict(x_um=self.mesh.x * 1e4, pH=-np.log10(free[:, 0]), Zn_M=x[:, ZN] * 1e3, Mn_M=x[:, MN] * 1e3,
                    S_M=x[:, SO] * 1e3, H_M=x[:, H] * 1e3, Zn_free_M=free[:, 1], phi2=x[:, P2],
                    vf_ZHS=nanc(x[:, ZH] * self.V[ZH]), vf_ZMO=nanc(x[:, ZM] * self.V[ZM]),
                    vf_MnO2=nanc(x[:, MO] * self.V[MO]), theta=nanc(th), eps=self.porosity(x))

    def inventory(self, x):
        """Moles of Zn, Mn, S and H in the cell (electrolyte + solids)."""
        p, m = self.p, self.mesh
        vol = m.area * m.dx
        eps = self.porosity(x)
        c = self.cath
        el = {k: float(np.sum(vol * eps * x[:, k])) for k in SPECIES}
        vc = vol[c]
        s = lambda a: float(np.sum(vc * a))
        zn = el[ZN] + s(p.z_ZMO * x[c, ZM] + self.n_ins(x[c, TH]) + 4.0 * x[c, ZH] + x[c, ZO] + x[c, ZX])
        mn = el[MN] + s(x[c, MO] + x[c, ZM])
        so = el[SO] + s(x[c, ZH])
        h = el[H] - s(6.0 * x[c, ZH] + 2.0 * (x[c, ZO] + x[c, ZX]) + 4.0 * (x[c, MO] + x[c, ZM]))
        return dict(Zn=zn, Mn=mn, S=so, H=h)
