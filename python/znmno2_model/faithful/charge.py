"""Faithful port of the original charge-line program (ZnMn02_v3, N = 5).

Reproduces the original's output byte for byte, defects included (see docs/deviations.md when written;
leads M-1 .. M-17 in the private plan). Written as a close transliteration of the Fortran so that every
floating-point operation happens in the same order:

- scalar `math` functions (the same libm as gfortran), never numpy's vectorised transcendental functions;
- single-precision literals and implicitly single-precision variables rounded with f32();
- the pH spline evaluated in float32 (numpy.float32 scalars);
- Newman's BAND with the legacy pivot, one linearized solve per time step (faithful/solver.py).

Indices follow the Fortran: node j = 1..NJ is array index j-1; unknown k = 1..5 is index k-1
(1 phi1, 2 phi2, 3 Zn2+, 4 Mn2+, 5 SO4 2-).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np


from .. import tables
from ..fortran_io import record
from . import solver
from .constants import Constants, Template, f32

HEADER = (" time Voltage Current mAh/g State C_Zn(avg) C_Mn(avg) C_SO4(avg) pH_phreeq(avg) pH_ZHS(avg) "
          "pH_standard(avg)  KMn8O16(avg) ZMC(avg) ZMC_max(avg) ZHS(avg) intercalation_Zn_to_Mn_ratio(avg)  "
          "R1_OCP(avg) R2_OCP(avg) R3_OCP(avg)  R1_Eta(avg) R2_Eta(avg) R3_Eta(avg)  R1_Exi(avg) R2_Exi(avg) "
          "R3_Exi(avg)  R1_I(avg) R2_I(avg) R3_I(Avg)  Area_KMn8O16 Area_ZMC Area_ZMC_max Area_ZHS  "
          "mAh/g_sim dom_react")

F32 = np.float32
EIGHT_THIRDS = f32(8.0 / 3.0)          # (8.0/3) and (8.0/3.0): single-precision constant folding
TWO_THIRDS = f32(2.0 / 3.0)
R2_U0 = float(F32(1.78) + F32(0.76))   # (1.78+0.76) folded in single precision
PH_SLOPE = float(F32(0.0592) * F32(4.0))  # 0.0592*(8.0/2.0) folded in single precision
C03_REF, C04_REF, C05_REF = f32(0.002), f32(0.0001), f32(0.0021)
ZN_REF = f32(0.001)
SOLID_LIMIT = f32(0.1)
ONE_SIXTH_SP = f32(1.0 / 6.0)
FD = f32(0.001)


def _exp(x: float) -> float:
    """exp with IEEE semantics (libm's exp; +inf on overflow, as in Fortran)."""
    try:
        return math.exp(x)
    except OverflowError:
        return math.inf


def _log(x: float) -> float:
    """log with IEEE semantics: -inf at 0, NaN for negative or NaN arguments."""
    if x > 0.0:
        return math.log(x)
    if x == 0.0:
        return -math.inf
    return math.nan


def _log10(x: float) -> float:
    if x > 0.0:
        return math.log10(x)
    if x == 0.0:
        return -math.inf
    return math.nan


def _pow(x: float, y: float) -> float:
    """x**y for real y with IEEE semantics (libm's pow; inf/NaN instead of exceptions or complex)."""
    try:
        return math.pow(x, y)
    except OverflowError:
        return math.inf if (x > 0 or float(y).is_integer() and int(y) % 2 == 0) else -math.inf
    except ValueError:
        if x == 0.0 and y < 0:
            return math.inf
        return math.nan


class Stop(Exception):
    """The original program executed STOP (or EXIT) here; the message is what it printed."""


def powi(x, n: int):
    """x**n for a constant integer n as gcc evaluates it (repeated squaring)."""
    r = 1.0 if not isinstance(x, np.floating) else x.dtype.type(1.0)
    while n:
        if n & 1:
            r = r * x
        x = x * x
        n >>= 1
    return r


# ---------------------------------------------------------------------------------------- pH spline
_ZN_PTS, _MN_PTS, _COEF = tables.ph_spline()
_OCP_THETA, _OCP_COEF = tables.r3_ocp()
_OCP_T = [float(x) for x in _OCP_THETA]
_OCP_C = [[float(x) for x in row] for row in _OCP_COEF]


def _log10f(x: np.float32) -> np.float32:
    return F32(_log10(float(x)))


def _basis(i: int, k: int, x: np.float32, t) -> np.float32:
    """Recursive B-spline basis (1-based knot index i), in single precision."""
    if i + k > len(t):                            # the original reads past the knot array here
        if x != x:
            return F32("nan")                     # NaN input: NaN whatever it reads
        raise Stop(f"pH spline knot read out of bounds (t({i + k})): not reproducible")
    if k == 1:
        return F32(1.0) if (t[i - 1] <= x and x < t[i]) else F32(0.0)
    a = (x - t[i - 1]) / (t[i + k - 2] - t[i - 1]) if t[i + k - 2] != t[i - 1] else F32(0.0)
    b = (t[i + k - 1] - x) / (t[i + k - 1] - t[i]) if t[i + k - 1] != t[i] else F32(0.0)
    return a * _basis(i, k - 1, x, t) + b * _basis(i + 1, k - 1, x, t)


def _bisect(x, arr) -> int:
    low, high = 1, len(arr)
    while low < high:
        mid = (low + high) // 2
        if x < arr[mid - 1]:
            high = mid
        else:
            low = mid + 1
    return low


def _find_interval(x, knots) -> int:
    n = len(knots)
    return max(min(_bisect(x, knots) - 1, n - 2), 0)


def manual_spline(conc1: np.float32, conc2: np.float32) -> np.float32:
    zn, mn = conc2, conc1                     # flipped on purpose in the original
    lg_zn, lg_mn = _log10f(zn), _log10f(mn)
    kx = ky = 3
    ix, iy = _find_interval(lg_zn, _ZN_PTS), _find_interval(lg_mn, _MN_PTS)
    bx = [(_basis(ix - kx + i, kx + 1, lg_zn, _ZN_PTS) if 1 <= ix - kx + i <= 153 - kx else F32(0.0)) for i in range(kx + 1)]
    by = [(_basis(iy - ky + i, ky + 1, lg_mn, _MN_PTS) if 1 <= iy - ky + i <= 153 - ky else F32(0.0)) for i in range(ky + 1)]
    ph = F32(0.0)
    for i in range(1, 5):
        for j in range(1, 5):
            ci, cj = ix - kx + i - 1, iy - ky + j - 1
            if 1 <= ci <= 149 and 1 <= cj <= 149:
                ph = ph + _COEF[ci - 1, cj - 1] * bx[i - 1] * by[j - 1]
    return ph


class ChargeModel:
    def __init__(self, template: Template):
        self.k = Constants().derive(template)
        k, v = self.k, self.k.values
        self.v = v
        N, NJ, S = k.N, k.NJ, k.SEP_NODE
        self.N, self.NJ, self.S = N, NJ, S
        self.RT = k.Rigc * k.Temp
        self._initial_condition()
        # run state
        self.time = 0.0
        self.delT = k.delT_nominal
        self.current = 0.0
        self.c_density = 0.0
        self.c_specific = 0.0
        self.mAhg = 0.0
        self.state = "D" if v["current"] >= 0.0 else "C"
        self.anode_pot = 0.0
        self.ramp_on = True
        self.ramp_count = 0
        self.ramp_initial = 0.0
        self.ramp_target = v["ramp_current_target"]
        self.last_write_time = 0                  # implicitly INTEGER in the original (D-4 pattern)
        self.delC = np.zeros((NJ, N))
        self.rows: list[str] = [HEADER]

    # ------------------------------------------------------------------ setup
    def _initial_condition(self):
        k, v = self.k, self.v
        N, NJ, S = self.N, self.NJ, self.S
        len_sep, xmax = v["len_sep"], v["xmax"]
        h_sep = len_sep / float(S - 2)
        h_cath = xmax / float(NJ - S - 1)
        xx, delx = np.zeros(NJ), np.zeros(NJ)
        for j in range(1, NJ + 1):
            if j == 1:
                x = 0.0
            elif j < S:
                x = h_sep * float(j - 1) - h_sep / 2.0
            elif j == S:
                x = len_sep
            elif j == NJ:
                x = xmax + len_sep
            else:
                x = len_sep + h_cath * float(j - S) - h_cath / 2.0
            xx[j - 1] = x
        for j in range(2, NJ):
            if j < S:
                delx[j - 1] = h_sep
            elif j > S:
                delx[j - 1] = h_cath
        delx[0] = delx[S - 1] = delx[NJ - 1] = 0.0
        self.xx, self.delx = xx, delx
        c = np.zeros((NJ, N))
        z = dict(KMn8O16=np.zeros(NJ), ZHS=np.zeros(NJ), ZMCx=np.zeros(NJ), ZMC_max=np.zeros(NJ),
                 por=np.zeros(NJ), tort=np.zeros(NJ), a_K=np.zeros(NJ), a_ZHS=np.zeros(NJ), a_ZMCx=np.zeros(NJ),
                 a_ZMC_max=np.zeros(NJ), znm=np.zeros(NJ), mnm=np.zeros(NJ), ratio=np.zeros(NJ), MW=np.zeros(NJ),
                 pH_phreeqc=np.full(NJ, k.pH_init), pH_ZHS=np.full(NJ, k.pH_init), pH_standard=np.full(NJ, k.pH_init))
        ci = v["c_initial"]
        for j in range(1, NJ + 1):
            i = j - 1
            c[i] = [v["Phi_1_init"], 0.0, ci[2], ci[3], ci[4]]
            if j < S:
                z["por"][i] = v["eps_sep"]
                z["tort"][i] = 2 * z["por"][i] ** (-0.5)
            else:
                z["KMn8O16"][i] = v["KMn8O16_init"]
                z["ZHS"][i] = v["ZHS_init"]
                z["ZMCx"][i] = v["ZMCx_init"]
                z["ZMC_max"][i] = v["ZMC_max_init"]
                z["por"][i] = v["eps"]
                z["tort"][i] = 2 * z["por"][i] ** (-0.5)
                z["a_K"][i] = 3.0 * z["KMn8O16"][i] / (k.density_KMn8O16 * k.xmax_c_KMn8O16)
                z["a_ZHS"][i] = 3.0 * z["ZHS"][i] / (k.density_ZHS * k.xmax_c_ZHS)
                z["a_ZMCx"][i] = 3.0 * z["ZMCx"][i] / (k.density_ZMC * k.xmax_c_ZMC)
                z["a_ZMC_max"][i] = 3.0 * z["ZMC_max"][i] / (k.density_ZMC * k.xmax_c_ZMC)
                z["ratio"][i] = v["ratio_initial"]
                z["MW"][i] = f32(65.38) * z["ratio"][i] + f32(54.938) + float(F32(15.999) * F32(2))
                z["znm"][i] = z["ZMCx"][i] * z["ratio"][i] / z["MW"][i]
                z["mnm"][i] = z["ZMCx"][i] / z["MW"][i]
        self.c, self.s = c, z
        self.diff_term = np.zeros((NJ, N))
        self.mig_term = np.zeros((NJ, N))
        self._transport_terms(range(1, NJ + 1))

    def _transport_terms(self, nodes):
        k, v = self.k, self.v
        sc = v["resevoir_scaling"]
        for j in nodes:
            i = j - 1
            p, t = self.s["por"][i], self.s["tort"][i]
            for ic in range(3, 6):
                d, zz = v["diff_ion"][ic - 1], v["z_ion"][ic - 1]
                if j < self.S:
                    self.diff_term[i, ic - 1] = (sc ** 2) * p * d / t
                    self.mig_term[i, ic - 1] = (sc ** 2) * p * zz * d * k.Fconst / (k.Rigc * k.Temp * t)
                else:
                    self.diff_term[i, ic - 1] = p * d / t
                    self.mig_term[i, ic - 1] = p * zz * d * k.Fconst / (k.Rigc * k.Temp * t)

    # ------------------------------------------------------------------ chemistry
    def zhs_pH(self, c_zn, c_so4):
        zn, so4 = c_zn * 1000.0, c_so4 * 1000.0
        c_h = _pow(powi(zn, 4) * so4 / self.v["K_sp"], 1.0 / 6.0)
        return -_log10(c_h)

    def eval_pH(self, zn: np.float32, mn: np.float32, j: int) -> np.float32:
        zn_m, mn_m = zn * F32(1000.0), mn * F32(1000.0)
        c_so4 = zn_m + mn_m
        ph_ksp = F32(-_log10(_pow(float(powi(zn_m, 4) * c_so4) / self.v["K_sp"], ONE_SIXTH_SP)))
        ph = manual_spline(zn_m, mn_m)
        if j > self.S and ph_ksp < ph:
            ph = ph_ksp
        return ph

    def r3_ocp_spline(self, theta: float) -> float:
        ind = 1
        while ind <= 51:
            if theta <= _OCP_T[ind - 1]:
                break
            ind += 1
        ind -= 1
        if theta != theta:                        # NaN: every comparison fails, the read runs off the
            return float("nan")                   # table, and the result is NaN whatever it reads
        if ind < 1 or ind > 50:
            raise Stop(f"R3 OCP spline read out of bounds (row {ind}, theta={theta!r}): not reproducible (M-17)")
        p1, p2, p3, p4 = _OCP_C[ind - 1]
        td = theta - _OCP_T[ind - 1]
        vs = p1 * powi(td, 3) + p2 * powi(td, 2) + p3 * td + p4
        return ((self.k.V_at_Zmin - self.k.V_at_Zmax) * vs) + self.k.V_at_Zmax

    def reaction_2(self, p1, p2, c3, c4, c5, j):
        k, v, s, i = self.k, self.v, self.s, j - 1
        out = [0.0] * 12
        xreact = v["Zmax"]
        ph = self.zhs_pH(c3, c5)
        aa = 0.5
        ac = 1.0 - aa
        rk = v["Rxn2_K"]
        n_e = 2.0 * (1.0 - xreact)
        rt = k.Rigc * k.Temp / (2 * k.Fconst)
        ocp = R2_U0 + (rt * (1.0 * (-1) * _log(c3 / C03_REF))) + (rt * (1.0 * (-2.0) * _log(c4 / C04_REF))) \
            - PH_SLOPE * ph
        exi = k.Fconst * rk \
            * _pow(c3 / C03_REF, 1.0 * (EIGHT_THIRDS - xreact) * aa / n_e) \
            * _pow(c4 / C04_REF, 1.0 * -1.0 * ac / n_e) \
            * _pow(c5 / C05_REF, 1.0 * TWO_THIRDS * aa / n_e)
        eta = p1 - p2 - ocp
        area = s["a_ZHS"][i] if self.state == "C" else s["a_ZMC_max"][i]
        irxn = area * exi * (_exp(aa * k.Fconst * eta / self.RT) - _exp(-ac * k.Fconst * eta / self.RT))
        nF = n_e * k.Fconst
        M = v["molar_mass_ZMC_max"]
        dt = self.delT

        def rates(irxn, with_h=True, h=None):
            return (((EIGHT_THIRDS - xreact) * irxn / nF), (-irxn / nF), (4 * irxn / nF if with_h else h),
                    (irxn / nF), (TWO_THIRDS * irxn / nF), (-TWO_THIRDS * irxn / nF))

        dzn, dmn, dh, dzmc, dso4, dzhs = rates(irxn)
        bl = k.BL_thickness
        if (-1.0 * dmn / area) > (c4 * k.diff_Mn / bl):
            irxn = ((c4 * k.diff_Mn / bl) / ((-1.0 * dmn) / area)) * irxn
            dzn, dmn, dh, dzmc, dso4, dzhs = rates(irxn)
        if (-1.0 * dzn / area) > (c3 * k.diff_Zn / bl):
            irxn = ((c3 * k.diff_Zn / bl) / ((-1.0 * dzn) / area)) * irxn
            dzn, dmn, dh, dzmc, dso4, dzhs = rates(irxn)
        if (-1.0 * dso4 / area) > (c5 * k.diff_SO4 / bl):
            irxn = ((c5 * k.diff_SO4 / bl) / ((-1.0 * dso4) / area)) * irxn
            dzn, dmn, dh, dzmc, dso4, dzhs = rates(irxn)
        if (-1.0 * (dzmc * M * dt) / s["ZMC_max"][i]) >= SOLID_LIMIT:
            irxn = irxn * abs(s["ZMC_max"][i] / (dzmc * M * dt)) * SOLID_LIMIT
            dzn, dmn, dh, dzmc, dso4, dzhs = rates(irxn, with_h=False, h=dh)
        if (-1.0 * (dzhs * M * dt) / s["ZHS"][i]) >= SOLID_LIMIT:
            irxn = irxn * abs(s["ZHS"][i] / (dzhs * M * dt)) * SOLID_LIMIT
            dzn, dmn, dh, dzmc, dso4, dzhs = rates(irxn, with_h=False, h=dh)
        if (s["ZMC_max"][i] + (dzmc * M * dt)) <= 0.0:
            irxn = irxn * abs(s["ZMC_max"][i] / (dzmc * M * dt)) * SOLID_LIMIT
            dzn, dmn, dh, dzmc, dso4, dzhs = rates(irxn)
        if (s["ZHS"][i] + (dzhs * M * dt)) <= 0.0:
            irxn = irxn * abs(s["ZHS"][i] / (dzhs * M * dt)) * SOLID_LIMIT
            dzn, dmn, dh, dzmc, dso4, dzhs = rates(irxn)
        out[0], out[1], out[2], out[3] = ocp, eta, exi, irxn
        out[4], out[5], out[6], out[7] = dzn, dmn, dso4, dh
        out[9], out[10] = dzhs, dzmc
        return out

    def reaction_3(self, p1, p2, c3, c4, c5, j):
        k, v, s, i = self.k, self.v, self.s, j - 1
        out = [0.0] * 12
        self.zhs_pH(c3, c5)                       # computed and unused in the original
        aa = 0.5
        ac = 1.0 - aa
        rk = v["Rxn3_K"]
        n_e = 2.0
        ratio = s["ratio"][i]
        theta = (ratio - v["Zmin"]) / (v["Zmax"] - v["Zmin"])
        if theta <= f32(-0.1):
            vint = f32(1.8)
        else:
            vint = f32(self.r3_ocp_spline(theta) - f32(0.762))   # Vint: implicitly REAL(4)
        nernst = ((k.Rigc * k.Temp / (n_e * k.Fconst)) * (1.0 * _log(c3 / C03_REF)))
        ocp = nernst + vint - self.anode_pot
        if ratio < v["Zmax"]:
            exi = k.Fconst * rk * _pow(c3 / C03_REF, 0.5 * ac) * _pow(ratio, aa) * _pow(v["Zmax"] - ratio, ac)
        else:
            exi = 1.0e-15
        eta = p1 - p2 - ocp
        area = s["a_ZMCx"][i]
        irxn = area * exi * (_exp(aa * k.Fconst * eta / self.RT) - _exp(-ac * k.Fconst * eta / self.RT))
        dzn = irxn / (n_e * k.Fconst)
        if (-1.0 * dzn / area) > (c3 * k.diff_Zn / k.BL_thickness):
            irxn = ((c5 * k.diff_Zn / k.BL_thickness) / abs(dzn)) * irxn        # c5: original typo (M-6)
            dzn = irxn / (n_e * k.Fconst)
        if (s["znm"][i] + (-1.0 * dzn * self.delT)) <= 0.0:
            irxn = irxn * abs(s["znm"][i] / (-1.0 * dzn * self.delT)) * SOLID_LIMIT
            dzn = irxn / (n_e * k.Fconst)
        out[0], out[1], out[2], out[3], out[4] = ocp, eta, exi, irxn, dzn
        out[11] = -1.0 * dzn
        return out

    def react_tot(self, p1, p2, c3, c4, c5, j):
        r2 = self.reaction_2(p1, p2, c3, c4, c5, j)
        r3 = self.reaction_3(p1, p2, c3, c4, c5, j)
        tot = [0.0] * 12
        for n in range(3, 12):
            acc = 0.0
            acc = 0.0 + acc                       # R1 (off)
            acc = r2[n] + acc
            acc = r3[n] + acc
            tot[n] = acc
        return tot

    def drx_dc(self, p1, p2, c3, c4, c5, j):
        d = [[0.0] * 4 for _ in range(5)]
        sp1 = p1 * FD
        sp2 = p1 * FD                              # the phi2 step uses phi1 (M-7)
        args = [p1, p2, c3, c4, c5]
        steps = [sp1, sp2, c3 * FD, c4 * FD, c5 * FD]
        for m in range(5):
            st = steps[m]
            if m >= 2 and args[m] <= st:
                a = list(args); a[m] = args[m] + st
                t1 = self.react_tot(*a, j)
                t2 = self.react_tot(*args, j)
                d[m] = [(t1[q] - t2[q]) / st for q in range(3, 7)]
            else:
                a = list(args); a[m] = args[m] + st
                b = list(args); b[m] = args[m] - st
                t1 = self.react_tot(*a, j)
                t2 = self.react_tot(*b, j)
                d[m] = [(t1[q] - t2[q]) / (2.0 * st) for q in range(3, 7)]
        return d

    def zn_anode_pot(self, c_zn):
        k = self.k
        return f32(-0.762) + (k.Rigc * k.Temp / (2.0 * k.Fconst)) * _log(c_zn / ZN_REF)

    # ------------------------------------------------------------------ one linearized step
    def assemble(self):
        """fillmat + ABDGXY for every node: the blocks of A delta_{j-1} + B delta_j + D delta_{j+1} = G."""
        k, v = self.k, self.v
        N, NJ, S = self.N, self.NJ, self.S
        A = np.zeros((NJ, N, N)); B = np.zeros((NJ, N, N)); D = np.zeros((NJ, N, N)); G = np.zeros((NJ, N))
        c, delx = self.c, self.delx
        z = v["z_ion"]
        alphaE = betaE = alphaW = betaW = 0.0
        for j in range(1, NJ + 1):
            i = j - 1
            dE = np.zeros((N, N)); dW = np.zeros((N, N)); fE = np.zeros((N, N)); fW = np.zeros((N, N))
            rj = np.zeros((N, N)); g = np.zeros(N)
            p1, p2, c3, c4, c5 = (float(x) for x in c[i])
            cE = np.zeros(N); cW = np.zeros(N); dcdxE = np.zeros(N); dcdxW = np.zeros(N)
            if j == 1:
                alphaE = delx[i] / (delx[i + 1] + delx[i])
                betaE = 2.0 / (delx[i] + delx[i + 1])
                for ic in range(N):
                    cE[ic] = alphaE * c[i + 1, ic] + (1.0 - alphaE) * c[i, ic]
                    dcdxE[ic] = betaE * (c[i + 1, ic] - c[i, ic])
                dE[0, 0] = -1.0
                g[0] = -(dE[0, 0] * dcdxE[0])
                g[1] = -(z[2] * c3) - (z[3] * c4) - (z[4] * c5)
                rj[1, 2], rj[1, 3], rj[1, 4] = z[2], z[3], z[4]
                dE[2, 2] = -1.0 * self.diff_term[i, 2]
                dE[2, 1] = -1.0 * self.mig_term[i, 2] * cE[2]
                fE[2, 2] = -1.0 * self.mig_term[i, 2] * dcdxE[1]
                g[2] = -1.0 * self.c_density / (z[2] * k.Fconst) + (dE[2, 2] * dcdxE[2] + fE[2, 2] * cE[2])
                for ic in range(4, N):            # Fortran ic = 4 .. N-1
                    q = ic - 1
                    dE[q, q] = -1.0 * self.diff_term[i, q]
                    dE[q, 1] = -1.0 * self.mig_term[i, q] * cE[q]
                    fE[q, q] = -1.0 * self.mig_term[i, q] * dcdxE[1]
                    g[q] = 0.0 + (dE[q, q] * dcdxE[q] + fE[q, q] * cE[q])
                g[4] = 0.0 - c[i, 1]
                rj[4, 1] = 1.0
                B[i] = rj - (1.0 - alphaE) * fE + betaE * dE
                D[i] = -alphaE * fE - betaE * dE
                G[i] = g
                continue
            alphaW = delx[i - 1] / (delx[i - 1] + delx[i])
            betaW = 2.0 / (delx[i - 1] + delx[i])
            if j < NJ:
                alphaE = delx[i] / (delx[i + 1] + delx[i])
                betaE = 2.0 / (delx[i] + delx[i + 1])
            for ic in range(N):
                cW[ic] = alphaW * c[i, ic] + (1.0 - alphaW) * c[i - 1, ic]
                dcdxW[ic] = betaW * (c[i, ic] - c[i - 1, ic])
                if j < NJ:
                    cE[ic] = alphaE * c[i + 1, ic] + (1.0 - alphaE) * c[i, ic]
                    dcdxE[ic] = betaE * (c[i + 1, ic] - c[i, ic])
            if j == S:
                dE[0, 0] = -1.0
                g[0] = -(dE[0, 0] * dcdxE[0])
                g[1] = -z[2] * c3 - z[3] * c4 - z[4] * c5
                rj[1, 2], rj[1, 3], rj[1, 4] = z[2], z[3], z[4]
                for q in range(2, N):
                    dW[q, q] = -1 * self.diff_term[i - 1, q]
                    dE[q, q] = -1 * self.diff_term[i + 1, q]
                    fW[q, q] = -1.0 * self.mig_term[i - 1, q] * dcdxW[1]
                    fE[q, q] = -1 * self.mig_term[i + 1, q] * dcdxE[1]
                    dW[q, 1] = -1 * self.mig_term[i - 1, q] * cW[q]
                    dE[q, 1] = -1 * self.mig_term[i + 1, q] * cE[q]
                    g[q] = 0.0 - (fW[q, q] * cW[q] + dW[q, q] * dcdxW[q]) + (fE[q, q] * cE[q] + dE[q, q] * dcdxE[q])
            elif j == NJ:
                dW[0, 0] = -(1.0 - self.s["por"][i]) * v["sigma"]
                g[0] = (self.c_density) - dW[0, 0] * dcdxW[0]
                g[1] = -z[2] * c3 - z[3] * c4 - z[4] * c5
                rj[1, 2], rj[1, 3], rj[1, 4] = z[2], z[3], z[4]
                for q in range(2, N):
                    dW[q, q] = -1.0 * self.diff_term[i, q]
                    dW[q, 1] = -1.0 * self.mig_term[i, q] * cW[q]
                    fW[q, q] = -1.0 * self.mig_term[i, q] * dcdxW[1]
                    g[q] = 0.0 - (dW[q, q] * dcdxW[q] + fW[q, q] * cW[q])
                A[i] = (1.0 - alphaW) * fW - betaW * dW
                B[i] = rj + betaW * dW + alphaW * fW
                G[i] = g
                continue
            elif j < S:
                por = self.s["por"][i]
                dE[0, 0] = -(1.0 - por) * k.sigma_sep
                dW[0, 0] = -(1.0 - por) * k.sigma_sep
                g[0] = 0.0 - (fW[0, 0] * cW[0] + dW[0, 0] * dcdxW[0]) + (fE[0, 0] * cE[0] + dE[0, 0] * dcdxE[0])
                g[1] = -z[2] * c3 - z[3] * c4 - z[4] * c5
                rj[1, 2], rj[1, 3], rj[1, 4] = z[2], z[3], z[4]
                for q in range(2, N):
                    dW[q, q] = -1 * self.diff_term[i, q]
                    dE[q, q] = -1 * self.diff_term[i, q]
                    fW[q, q] = -1 * self.mig_term[i, q] * dcdxW[1]
                    fE[q, q] = -1 * self.mig_term[i, q] * dcdxE[1]
                    dW[q, 1] = -1 * self.mig_term[i, q] * cW[q]
                    dE[q, 1] = -1 * self.mig_term[i, q] * cE[q]
                    g[q] = - (fW[q, q] * cW[q] + dW[q, q] * dcdxW[q]) + (fE[q, q] * cE[q] + dE[q, q] * dcdxE[q])
                for qe in range(2, N):
                    for qv in range(N):
                        rj[qe, qv] = -(por / self.delT) * delx[i] if qe == qv else 0.0
            else:                                  # cathode interior
                rxn = self.react_tot(p1, p2, c3, c4, c5, j)
                drx = self.drx_dc(p1, p2, c3, c4, c5, j)
                por = self.s["por"][i]
                dE[0, 0] = -(1.0 - por) * v["sigma"]
                dW[0, 0] = -(1.0 - por) * v["sigma"]
                for q in range(N):
                    rj[0, q] = -drx[q][0] * delx[i]
                g[0] = rxn[3] * delx[i] - (fW[0, 0] * cW[0] + dW[0, 0] * dcdxW[0]) + (fE[0, 0] * cE[0] + dE[0, 0] * dcdxE[0])
                g[1] = -z[2] * c3 - z[3] * c4 - z[4] * c5
                rj[1, 2], rj[1, 3], rj[1, 4] = z[2], z[3], z[4]
                for q in range(2, N):
                    dW[q, q] = -1 * self.diff_term[i, q]
                    dE[q, q] = -1 * self.diff_term[i, q]
                    fW[q, q] = -1 * self.mig_term[i, q] * dcdxW[1]
                    fE[q, q] = -1 * self.mig_term[i, q] * dcdxE[1]
                    dW[q, 1] = -1 * self.mig_term[i, q] * cW[q]
                    dE[q, 1] = -1 * self.mig_term[i, q] * cE[q]
                    g[q] = -(rxn[q + 2]) * delx[i] - (fW[q, q] * cW[q] + dW[q, q] * dcdxW[q]) + (fE[q, q] * cE[q] + dE[q, q] * dcdxE[q])
                for qe in range(2, N):
                    for qv in range(N):
                        if qe == qv:
                            rj[qe, qv] = drx[qv][qe - 1] * delx[i] - (por / self.delT) * delx[i]
                        else:
                            rj[qe, qv] = drx[qv][qe - 1] * delx[i]
            A[i] = (1.0 - alphaW) * fW - betaW * dW
            B[i] = rj + betaW * dW + alphaW * fW - (1.0 - alphaE) * fE + betaE * dE
            D[i] = -alphaE * fE - betaE * dE
            G[i] = g
        self.last_coeffs = (dE, dW, fE, fW, rj, g)
        return A, B, D, G

    # ------------------------------------------------------------------ after the solve
    def update_other_variables(self, it: int):
        k, v, s = self.k, self.v, self.s
        NJ, S = self.NJ, self.S
        M = v["molar_mass_ZMC_max"]
        cutoff_theta = f32(0.985)
        c, dC = self.c, self.delC
        for j in range(S + 1, NJ):
            i = j - 1
            r = self.react_tot(*(float(c[i, q] - (dC[i, q] / 2)) for q in range(5)), j)
            s["KMn8O16"][i] = s["KMn8O16"][i] + (k.molar_mass_KMn8O16 * r[8] * self.delT)
            s["ZHS"][i] = s["ZHS"][i] + (k.molar_mass_ZHS * r[9] * self.delT)
            s["ZMC_max"][i] = s["ZMC_max"][i] + (M * r[10] * self.delT)
            s["znm"][i] = s["znm"][i] + (r[11] * self.delT)
            s["ratio"][i] = s["znm"][i] / s["mnm"][i]
            if s["ratio"][i] >= cutoff_theta * v["Zmax"]:
                zn_t, mn_t = s["znm"][i], s["mnm"][i]
                mn1 = (zn_t - v["Zmax"] * mn_t) / (cutoff_theta * v["Zmax"] - v["Zmax"])
                mn2 = mn_t - mn1
                zn1 = cutoff_theta * v["Zmax"] * mn1
                s["znm"][i], s["mnm"][i] = zn1, mn1
                s["ratio"][i] = s["znm"][i] / s["mnm"][i]
                s["ZMC_max"][i] = s["ZMC_max"][i] + (mn2 * M)
            s["MW"][i] = f32(65.38) * s["ratio"][i] + f32(54.93) + float(F32(15.999) * F32(2))
            s["ZMCx"][i] = s["MW"][i] * s["mnm"][i]
            if s["KMn8O16"][i] < 0.0:
                raise Stop("KMn8O16 is negative")
            if s["ZHS"][i] < 0.0:
                s["ZHS"][i] = 1.0e-50
            if s["ZMCx"][i] < 0.0:
                s["ZMCx"][i] = 1.0e-50
            if s["ZMC_max"][i] < 0.0:
                s["ZMC_max"][i] = 1.0e-50
            if s["znm"][i] < 0.0:
                s["znm"][i] = 1.7e-50
            if s["mnm"][i] < 0.0:
                s["mnm"][i] = 1.0e-50
            s["a_K"][i] = 3.0 * s["KMn8O16"][i] / (k.density_KMn8O16 * k.xmax_c_KMn8O16)
            s["a_ZMCx"][i] = 3.0 * s["ZMCx"][i] / (k.density_ZMC * k.xmax_c_ZMC)
            s["a_ZMC_max"][i] = 3.0 * s["ZMC_max"][i] / (k.density_ZMC * k.xmax_c_ZMC)
            s["a_ZHS"][i] = 3.0 * s["ZHS"][i] / (k.density_ZHS * k.xmax_c_ZHS)
            por_past = s["por"][i]
            for q in (2, 3, 4):
                c[i, q] = c[i, q] * por_past / s["por"][i]
            s["pH_phreeqc"][i] = float(self.eval_pH(F32(c[i, 2]), F32(c[i, 3]), j))
            s["pH_ZHS"][i] = self.zhs_pH(float(c[i, 2]), float(c[i, 4]))
            # cprev(6,j) with N = 5 reads cprev(1,j+1) (M-14); single-precision log10
            s["pH_standard"][i] = float(F32(-1.0) * _log10f(F32(1000 * c[i + 1, 0])))
        # after the loop the index j is NJ (the original's leftover loop variable, M-16)
        for key in ("KMn8O16", "ZMCx", "ZMC_max", "ZHS", "znm", "mnm", "ratio", "MW", "por", "tort",
                    "a_K", "a_ZMCx", "a_ZMC_max", "a_ZHS"):
            s[key][S - 1] = s[key][S]
        for key in ("pH_phreeqc", "pH_ZHS", "pH_standard"):
            s[key][NJ - 1] = s[key][S]
        for key in ("KMn8O16", "ZMCx", "ZMC_max", "ZHS", "znm", "mnm", "ratio", "MW", "por", "tort",
                    "a_K", "a_ZMCx", "a_ZMC_max", "a_ZHS"):
            s[key][NJ - 1] = s[key][NJ - 2]
        for key in ("pH_phreeqc", "pH_ZHS", "pH_standard"):
            s[key][NJ - 1] = s[key][NJ - 2]
        self.anode_pot = self.zn_anode_pot(float(c[0, 2]))

    # ------------------------------------------------------------------ output
    def write_row(self):
        NJ, S = self.NJ, self.S
        n = NJ - S - 1
        idx = range(S, NJ - 1)                    # Fortran SEP_NODE+1 .. NJ-1
        r_all = {2: [[0.0] * 4 for _ in range(NJ)], 3: [[0.0] * 4 for _ in range(NJ)]}
        for j in range(S + 1, NJ + 1):
            i = j - 1
            if j > S:
                args = [float(x) for x in self.c[i]]
                r2 = self.reaction_2(*args, j)
                r3 = self.reaction_3(*args, j)
                r_all[2][i] = r2[:4]
                r_all[3][i] = r3[:4]

        def mean_div(arr):                        # SUM(x(a:b)/n): divide each, then add in order
            acc = 0.0
            for i in idx:
                acc = acc + arr[i] / n
            return acc

        def sum_div(vals):                        # SUM(x(a:b))/n
            acc = 0.0
            for x in vals:
                acc = acc + x
            return acc / n

        s, c = self.s, self.c
        r1 = [sum_div([0.0 for _ in idx])] * 4
        r2 = [sum_div([r_all[2][i][q] for i in idx]) for q in range(4)]
        r3 = [sum_div([r_all[3][i][q] for i in idx]) for q in range(4)]
        temp = [abs(r1[3]), abs(r2[3]), abs(r3[3]), abs(sum_div([0.0 for _ in idx]))]
        dom = max(range(4), key=lambda q: (temp[q], -q)) + 1
        items = [self.time, float(c[NJ - 1, 0]), self.current, self.mAhg, self.state,
                 mean_div(c[:, 2]), mean_div(c[:, 3]), mean_div(c[:, 4]),
                 mean_div(s["pH_phreeqc"]), mean_div(s["pH_ZHS"]), mean_div(s["pH_standard"]),
                 mean_div(s["KMn8O16"]), mean_div(s["ZMCx"]), mean_div(s["ZMC_max"]), mean_div(s["ZHS"]),
                 mean_div(s["ratio"]),
                 r1[0], r2[0], r3[0], r1[1], r2[1], r3[1], r1[2], r2[2], r3[2], r1[3], r2[3], r3[3],
                 sum_div([s["a_K"][i] for i in idx]), sum_div([s["a_ZMCx"][i] for i in idx]),
                 sum_div([s["a_ZMC_max"][i] for i in idx]), sum_div([s["a_ZHS"][i] for i in idx]),
                 self.mAhg * self.v["AM_Grams"] / self.v["AM_Grams_sim"], dom]
        self.rows.append(record(*items))

    # ------------------------------------------------------------------ main loop
    def current_ramp(self):
        if self.ramp_on:
            if self.ramp_count == 0:
                self.ramp_initial = self.current
            self.delT = self.k.ramp_delT
            self.current = (float(self.ramp_count) * (self.ramp_target - self.ramp_initial) / float(self.k.ramp_iters)) \
                + self.ramp_initial
            self.ramp_count += 1
            if self.ramp_count == self.k.ramp_iters:
                self.current = self.ramp_target
                self.ramp_on = False
                self.ramp_count = 0
                self.delT = self.k.delT_nominal
        self.c_density = self.current / self.v["Area_CS"]
        self.c_specific = self.current / self.v["AM_Grams"]

    def _maybe_write(self, it):
        if it >= self.k.ramp_iters:
            if (self.time - self.k.write_time_density) >= self.last_write_time:
                self.write_row()
                self.last_write_time = int(self.time)
        elif it == 1:
            self.write_row()
            self.last_write_time = int(self.time)

    def run(self, max_steps: int | None = None, on_step=None) -> str:
        """Run to an exit; returns the reason the original printed."""
        NJ = self.NJ
        it = 0
        while True:
            it += 1
            if max_steps is not None and it > max_steps:
                return "max_steps"
            c = self.c
            if (c[NJ - 1, 1] >= 99.0) and self.state == "C":
                return "EXIT BECAUSE END OF CHARGE"
            if math.isnan(self.delC[0, 0]):
                return "EXIT BECAUSE delC ISNAN"
            if self.time >= 99.0 * 3600.0:
                return "EXIT BECAUSE END OF SIMULATION TIME"
            if (c[NJ - 1, 0] - c[NJ - 1, 1]) <= f32(0.8):
                return "EXIT BECAUSE LOWER VOLTAGE CUTOFF"
            if c[NJ - 1, 0] >= f32(1.8):
                return "EXIT BECAUSE Upper VOLTAGE CUTOFF"
            self.current_ramp()
            if self.current > 1.0e-10:
                self.state = "D"
            elif self.current <= -1.0e-10:
                self.state = "C"
            else:
                self.state = "R"
            if it >= self.k.ramp_iters:
                if (self.time - self.k.write_time_density) >= self.last_write_time:
                    self.write_row()
                    self.last_write_time = int(self.time)
            elif it == 1:
                self.write_row()
                self.last_write_time = int(self.time)
            A, B, D, G = self.assemble()
            self.delC = solver.solve(A, B, D, G)
            self.c = self.c + self.delC
            try:
                self.update_other_variables(it)
            except Stop as e:
                return str(e)
            if np.isnan(self.c).any():
                return "NaN in cprev"
            if np.isnan(self.delC).any():
                return "NaN in delC"
            if any(np.isnan(x).any() for x in self.last_coeffs):
                return "NaN in coefficients"
            if (self.time - self.k.write_time_density) >= self.last_write_time:
                self.write_row()
                self.last_write_time = int(self.time)
            self.time = self.time + self.delT
            if self.state == "D":
                self.mAhg = self.mAhg + 1000.0 * self.c_specific * self.delT / 3600.0
            elif self.state == "C":
                self.mAhg = self.mAhg - 1000.0 * self.c_specific * self.delT / 3600.0
            if on_step is not None:
                on_step(it, self)

    def write(self, path):
        Path(path).write_text("\n".join(self.rows) + "\n")
