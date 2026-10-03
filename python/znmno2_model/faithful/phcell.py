"""Faithful port of the original pH-cell / GITT line (ZnMn02_v2.1_GITT_no_probe, N = 6, H+ transported).

Same conventions as charge.py: a close transliteration of the Fortran, evaluated operation by operation
in the original's order, with single-precision literals and implicitly single variables rounded with
f32(). Unknowns k = 1..6: 1 phi1, 2 phi2, 3 Zn2+, 4 Mn2+, 5 SO4 2-, 6 H+.

The archived runs switch on R5 only (MnO2 dissolution with the ZHS pH correction), so R1-R4 are not
ported; a run that switches them on raises NotImplementedError. GITT cycling is off in these runs.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..fortran_io import record
from . import solver
from .charge import F32, Stop, _exp, _log, _log10, _log10f, _pow, manual_spline, powi
from .constants import f32, f32_mul, f32_pow10

HEADER = (" time Voltage Current mAh/g State  C_Zn(avg) C_Mn(avg) C_SO4(avg) C_H(avg)  pH_phreeq(avg) "
          "pH_ZHS(avg) pH_standard(avg)  KMn8O16(avg) ZMC(avg) ZMC_max(avg) ZHS(avg) MnO2(avg) "
          "intercalation_Zn_to_Mn_ratio(avg)  R1_OCP(avg) R2_OCP(avg) R3_OCP(avg) R4_OCP(avg) R5_OCP(avg)  "
          "R1_Eta(avg) R2_Eta(avg) R3_Eta(avg) R4_Eta(avg) R5_Eta(avg)  R1_Exi(avg) R2_Exi(avg) R3_Exi(avg) "
          "R4_Exi(avg) R5_Exi(avg)  R1_I(avg) R2_I(avg) R3_I(Avg) R4_I(avg) R5_I(avg)  Area_KMn8O16 Area_ZMC "
          "Area_ZMC_max Area_ZHS Area_MnO2  mAh/g_sim dom_react")

ONE_SIXTH_SP = f32(1.0 / 6.0)
FD = f32(0.001)
ZN_REF = f32(0.001)
CONC_FLOOR = 1.0e-20
R5_C04_REF = 0.0001                 # 0.0001d0: double in Reaction_5
ONE_POINT_ONE = f32(1.1)


@dataclass
class PhTemplate:
    """One row of the run generator's parameter file (values as the text that was substituted)."""
    RXNK_5: str
    FRACTION_KMNO2: str


@dataclass
class PhConstants:
    N: int = 6
    NJ: int = 134
    SEP_NODE: int = 63
    Rigc: float = f32(8.314)
    Temp: float = 298.0
    Fconst: float = 96485.0
    diff_Zn: float = 7.15e-6
    diff_Mn: float = 6.88e-6
    diff_H: float = 9.0e-5
    diff_SO4: float = 1.07e-5
    molar_mass_KMn8O16: float = f32(734.59)
    density_KMn8O16: float = f32(5.03)
    molar_mass_MnO2: float = f32(86.9368)
    density_MnO2: float = f32(5.03)
    molar_mass_ZHS: float = f32(549.819)
    density_ZHS: float = f32(2.67)
    density_ZMC: float = 5.0
    xmax_c: float = 200.0e-5           # the same for every solid
    K_sp: float = 7.0e-26
    BL_thickness: float = 0.5 * 1.0e-4
    sigma_sep: float = 1.0e-20
    ramp_iters: int = 1000
    ramp_delT: float = 1.0e-9
    delT_nominal: float = 1.0
    write_time_density: float = 10.0

    def derive(self, t: PhTemplate) -> dict:
        v: dict = {}
        v["Rxn5_K"] = f32_pow10(t.RXNK_5)
        v["Rxn2_Mn_correction"] = 1.0
        v["Zmin"] = f32(0.35)
        v["Zmax"] = f32(0.65)
        v["ratio_initial"] = v["Zmax"] * f32(0.9999)
        cZn, cMn, cH = f32(0.002), f32(0.00005), f32(0.000543)
        ph_init = 5.5
        if (-1.0 * _log10(1000 * cH)) <= ph_init:
            ph_init = -1.0 * _log10(1000 * cH)
        v["pH_init"] = ph_init
        v["sigma"] = f32(0.01)
        v["xmax"] = f32(0.0218)
        v["Area_CS"] = f32(0.178134094)
        v["len_sep"] = 600.0 * 1.0e-4
        v["eps_sep"] = f32(0.9)
        v["current"] = float(F32(0.107) / F32(1000.0))
        v["eps"] = f32(0.815)
        v["volfrac_KMn8O16"] = f32(0.00001)
        v["volfrac_MnO2"] = f32_mul("0.0616", t.FRACTION_KMNO2)
        v["volfrac_ZHS"] = f32(0.000001)
        v["volfrac_ZMC_max"] = f32(0.00001)
        v["volfrac_ZMCx"] = v["volfrac_ZMC_max"] * f32(0.01)
        v["volfrac_inert"] = (1.0 - v["eps"] - v["volfrac_KMn8O16"] - v["volfrac_ZHS"] - v["volfrac_ZMCx"]
                              - v["volfrac_ZMC_max"] - v["volfrac_MnO2"])
        v["AM_Grams"] = f32(0.00108)
        v["AM_Grams_sim"] = v["volfrac_ZMCx"] * v["Area_CS"] * v["xmax"] * self.density_ZMC
        v["KMn8O16_init"] = v["volfrac_KMn8O16"] * self.density_KMn8O16
        v["MnO2_init"] = v["volfrac_MnO2"] * self.density_MnO2
        v["ZHS_init"] = v["volfrac_ZHS"] * self.density_ZHS
        v["ZMCx_init"] = v["volfrac_ZMCx"] * self.density_ZMC
        v["ZMC_max_init"] = v["volfrac_ZMC_max"] * self.density_ZMC
        v["Phi_1_init"] = f32(1.1)
        v["molar_mass_ZMC_max"] = f32(65.38) * v["Zmax"] + f32(54.93) + f32_mul(15.999, 2)
        z = [0.0, 0.0, 2.0, 2.0, -2.0, 1.0]
        c = [0.0, 0.0, cZn, cMn, 0.0, cH]
        c[4] = -1 * ((z[2] * c[2]) + (z[3] * c[3]) + (z[5] * c[5])) / z[4]
        v["z_ion"], v["c_initial"] = z, c
        v["diff_ion"] = [0.0, 0.0, self.diff_Zn, self.diff_Mn, self.diff_SO4, self.diff_H]
        v["ramp_current_target"] = v["current"]
        # reservoir: resevoir_vol = 0 gives a negative scaling, which the program resets to 1.0
        v["resevoir_scaling"] = 1.0
        return v


class PhCellModel:
    def __init__(self, template: PhTemplate, reactions_on=(5,)):
        if tuple(reactions_on) != (5,):
            raise NotImplementedError("only R5 is ported (the archived pH-cell runs switch on R5 only)")
        self.k = k = PhConstants()
        self.v = v = k.derive(template)
        self.N, self.NJ, self.S = k.N, k.NJ, k.SEP_NODE
        self.RT = k.Rigc * k.Temp
        self._initial_condition()
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
        self.last_write_time = 0                  # implicitly INTEGER in the original
        self.delC = np.zeros((self.NJ, self.N))
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
        keys = ("KMn8O16", "MnO2", "ZHS", "ZMCx", "ZMC_max", "por", "tort", "a_K", "a_MnO2", "a_ZHS", "a_ZMCx",
                "a_ZMC_max", "znm", "mnm", "ratio", "MW")
        s = {key: np.zeros(NJ) for key in keys}
        for key in ("pH_phreeqc", "pH_ZHS", "pH_standard"):
            s[key] = np.full(NJ, v["pH_init"])
        ci = v["c_initial"]
        xc = k.xmax_c
        for j in range(1, NJ + 1):
            i = j - 1
            c[i] = [v["Phi_1_init"], 0.0, ci[2], ci[3], ci[4], ci[5]]
            if j < S:
                s["por"][i] = v["eps_sep"]
                s["tort"][i] = 2 * _pow(s["por"][i], -0.5)
            else:
                s["KMn8O16"][i] = v["KMn8O16_init"]
                s["MnO2"][i] = v["MnO2_init"]
                s["ZHS"][i] = v["ZHS_init"]
                s["ZMCx"][i] = v["ZMCx_init"]
                s["ZMC_max"][i] = v["ZMC_max_init"]
                s["por"][i] = v["eps"]
                s["tort"][i] = 2 * _pow(s["por"][i], -0.5)
                s["a_K"][i] = 3.0 * s["KMn8O16"][i] / (k.density_KMn8O16 * xc)
                s["a_MnO2"][i] = 3.0 * s["MnO2"][i] / (k.density_MnO2 * xc)
                s["a_ZHS"][i] = 3.0 * s["ZHS"][i] / (k.density_ZHS * xc)
                s["a_ZMCx"][i] = 3.0 * s["ZMCx"][i] / (k.density_ZMC * xc)
                s["a_ZMC_max"][i] = 3.0 * s["ZMC_max"][i] / (k.density_ZMC * xc)
                s["ratio"][i] = v["ratio_initial"]
                s["MW"][i] = f32(65.38) * s["ratio"][i] + f32(54.938) + float(F32(15.999) * F32(2))
                s["znm"][i] = s["ZMCx"][i] * s["ratio"][i] / s["MW"][i]
                s["mnm"][i] = s["ZMCx"][i] / s["MW"][i]
        self.c, self.s = c, s
        self.diff_term = np.zeros((NJ, N))
        self.mig_term = np.zeros((NJ, N))
        self._transport_terms(range(1, NJ + 1))

    def _transport_terms(self, nodes):
        """diff_term, mig_term at the given nodes. The reservoir scaling is 1.0, so the original's two
        branches (inside/outside the probe region) give the same bits."""
        k, v = self.k, self.v
        for j in nodes:
            i = j - 1
            p, t = self.s["por"][i], self.s["tort"][i]
            for ic in range(3, self.N + 1):
                d, zz = v["diff_ion"][ic - 1], v["z_ion"][ic - 1]
                self.diff_term[i, ic - 1] = p * d / t
                self.mig_term[i, ic - 1] = p * zz * d * k.Fconst / (k.Rigc * k.Temp * t)

    # ------------------------------------------------------------------ chemistry
    def zhs_pH(self, c_zn, c_so4):
        zn, so4 = c_zn * 1000.0, c_so4 * 1000.0
        c_h = _pow(powi(zn, 4) * so4 / self.k.K_sp, -1.0 / 6.0)
        return -_log10(c_h)

    def eval_pH(self, zn: np.float32, mn: np.float32, j: int) -> np.float32:
        zn_m, mn_m = zn * F32(1000.0), mn * F32(1000.0)
        c_so4 = zn_m + mn_m
        ph_ksp = F32(-_log10(_pow(float(powi(zn_m, 4) * c_so4) / self.k.K_sp, ONE_SIXTH_SP)))
        ph = manual_spline(zn_m, mn_m)
        if j > self.S and ph_ksp < ph:
            ph = ph_ksp
        return ph

    def zn_anode_pot(self, c_zn):
        k = self.k
        return f32(-0.762) + (k.Rigc * k.Temp / (2.0 * k.Fconst)) * _log(c_zn / ZN_REF)

    def reaction_5(self, p1, p2, c3, c4, c5, c6, j):
        k, v, s, i = self.k, self.v, self.s, j - 1
        out = [0.0] * 13
        F, RT, dt, bl = k.Fconst, self.RT, self.delT, k.BL_thickness
        ph_std = -_log10(1000.0 * c6)
        ph_zhs = self.zhs_pH(c3, c5)
        aa = 0.5
        ac = 1.0 - aa
        uref = 1.2225
        uan = self.anode_pot
        rk = v["Rxn5_K"]
        n_e = 2.0
        nF = n_e * F
        area = s["a_MnO2"][i]
        M_MnO2 = k.molar_mass_MnO2
        limit = 0.1
        mn_corr = v["Rxn2_Mn_correction"]

        def ocp_eta_irxn(ph):
            ocp = (uref - uan) + (k.Rigc * k.Temp / nF) * (-_log(c4 / R5_C04_REF)) - 0.0592 * (4.0 / n_e) * ph
            exi = F * rk * _pow(c4 / R5_C04_REF, -ac * mn_corr / n_e)
            eta = p1 - p2 - ocp
            irxn = area * exi * (_exp(aa * F * eta / RT) - _exp(-ac * F * eta / RT))
            return ocp, exi, eta, irxn

        if ph_std < ph_zhs:
            ocp, exi, eta, irxn = ocp_eta_irxn(ph_std)
            dmn, dh, dmno2 = -irxn / nF, 4.0 * irxn / nF, irxn / nF
            dzhs = dzn = dso4 = 0.0
            if (-dmn / area) > (c4 * k.diff_Mn / bl):
                irxn = ((c4 * k.diff_Mn / bl) / ((-dmn) / area)) * irxn
                dmn, dh, dmno2 = -irxn / nF, 4.0 * irxn / nF, irxn / nF
            if (-dmno2 * M_MnO2 * dt) / s["MnO2"][i] >= limit:
                irxn = irxn * abs(s["MnO2"][i] / (dmno2 * M_MnO2 * dt)) * limit
                dmn, dh, dmno2 = -irxn / nF, 4.0 * irxn / nF, irxn / nF
            h_post = (c6 * s["por"][i] * self.delx[i] * v["Area_CS"]) + (dh * dt)
            zhs_h = _pow(10.0, -ph_zhs) / 1000.0
            if (h_post < zhs_h) and (ONE_POINT_ONE * ph_std >= ph_zhs):
                h_needed = zhs_h - h_post
                h_rate = h_needed / dt
                r_zhs = h_rate / 6.0
                dh = dh + h_rate
                dzn = dzn - 4.0 * r_zhs
                dso4 = dso4 - 1.0 * r_zhs
                dzhs = dzhs + r_zhs
        else:
            ocp, exi, eta, irxn = ocp_eta_irxn(ph_zhs)

            def rates(irxn):
                return (-irxn / nF, irxn / nF, 0.0, -(4.0 * irxn / nF) / 6.0, (8.0 / 3.0) * irxn / nF,
                        (2.0 / 3.0) * irxn / nF)

            dmn, dmno2, dh, dzhs, dzn, dso4 = rates(irxn)
            if (-dzn / area) > (c3 * k.diff_Zn / bl):
                irxn = ((c3 * k.diff_Zn / bl) / ((-dzn) / area)) * irxn
                dmn, dmno2, dh, dzhs, dzn, dso4 = rates(irxn)
            if (-dso4 / area) > (c5 * k.diff_SO4 / bl):
                irxn = ((c5 * k.diff_SO4 / bl) / ((-dso4) / area)) * irxn
                dmn, dmno2, dh, dzhs, dzn, dso4 = rates(irxn)
            if (-dmno2 * M_MnO2 * dt) / s["MnO2"][i] >= limit:
                irxn = irxn * abs(s["MnO2"][i] / (dmno2 * M_MnO2 * dt)) * limit
                dmn, dmno2, dh, dzhs, dzn, dso4 = rates(irxn)
            if (-dzhs * M_MnO2 * dt) / s["ZHS"][i] >= limit:
                irxn = irxn * abs(s["ZHS"][i] / (dzhs * M_MnO2 * dt)) * limit
                dmn, dmno2, dh, dzhs, dzn, dso4 = rates(irxn)
        out[0], out[1], out[2], out[3] = ocp, eta, exi, irxn
        out[4], out[5], out[6], out[7] = dzn, dmn, dso4, dh
        out[9], out[12] = dzhs, dmno2
        return out

    def react_tot(self, p1, p2, c3, c4, c5, c6, j):
        r5 = self.reaction_5(p1, p2, c3, c4, c5, c6, j)
        tot = [0.0] * 13
        for n in range(3, 13):                    # R_Matrix rows 1..9 added in order; only row 5 is set
            acc = 0.0
            for _ in range(4):
                acc = 0.0 + acc
            acc = r5[n] + acc
            for _ in range(4):
                acc = 0.0 + acc
            tot[n] = acc
        return tot

    def drx_dc(self, p1, p2, c3, c4, c5, c6, j):
        args = [p1, p2, c3, c4, c5, c6]
        steps = [x * FD for x in args]
        steps = [1.0e-6 if st == 0.0 else st for st in steps]
        d = []
        for m in range(6):
            st = steps[m]
            a = list(args); a[m] = args[m] + st
            if m >= 2 and args[m] <= st:
                t1 = self.react_tot(*a, j)
                t2 = self.react_tot(*args, j)
                d.append([(t1[q] - t2[q]) / st for q in range(3, 8)])
            else:
                b = list(args); b[m] = args[m] - st
                t1 = self.react_tot(*a, j)
                t2 = self.react_tot(*b, j)
                d.append([(t1[q] - t2[q]) / (2.0 * st) for q in range(3, 8)])
        return d

    # ------------------------------------------------------------------ one linearized step
    def assemble(self):
        """fillmat + ABDGXY for every node: the blocks of A delta_{j-1} + B delta_j + D delta_{j+1} = G."""
        k, v = self.k, self.v
        N, NJ, S = self.N, self.NJ, self.S
        A = np.zeros((NJ, N, N)); B = np.zeros((NJ, N, N)); D = np.zeros((NJ, N, N)); G = np.zeros((NJ, N))
        c, delx, por = self.c, self.delx, self.s["por"]
        z = v["z_ion"]
        dt_, mt_ = self.diff_term, self.mig_term
        alphaE = betaE = alphaW = betaW = 0.0

        def charge_row(g, rj, cc):
            g[1] = -(z[2] * cc[2]) - (z[3] * cc[3]) - (z[4] * cc[4]) - (z[5] * cc[5])
            rj[1, 2], rj[1, 3], rj[1, 4], rj[1, 5] = z[2], z[3], z[4], z[5]

        for j in range(1, NJ + 1):
            i = j - 1
            dE = np.zeros((N, N)); dW = np.zeros((N, N)); fE = np.zeros((N, N)); fW = np.zeros((N, N))
            rj = np.zeros((N, N)); g = np.zeros(N)
            cc = [float(x) for x in c[i]]
            cE = np.zeros(N); cW = np.zeros(N); dcdxE = np.zeros(N); dcdxW = np.zeros(N)
            if j == 1:
                alphaE = delx[i] / (delx[i + 1] + delx[i])
                betaE = 2.0 / (delx[i] + delx[i + 1])
                for ic in range(N):
                    cE[ic] = alphaE * c[i + 1, ic] + (1.0 - alphaE) * c[i, ic]
                    dcdxE[ic] = betaE * (c[i + 1, ic] - c[i, ic])
                dE[0, 0] = -1.0
                g[0] = -(dE[0, 0] * dcdxE[0])
                g[1] = -(z[2] * cc[2]) - (z[3] * cc[3]) - (z[4] * cc[4]) - (z[5] * cc[5])
                rj[1, 2], rj[1, 3], rj[1, 4], rj[1, 5] = z[2], z[3], z[4], z[5]
                dE[2, 2] = -1.0 * dt_[i, 2]
                dE[2, 1] = -1.0 * mt_[i, 2] * cE[2]
                fE[2, 2] = -1.0 * mt_[i, 2] * dcdxE[1]
                g[2] = -1.0 * self.c_density / (z[2] * k.Fconst) + (dE[2, 2] * dcdxE[2] + fE[2, 2] * cE[2])
                for q in range(3, N - 1):          # Fortran ic = 4 .. N-1
                    dE[q, q] = -1.0 * dt_[i, q]
                    dE[q, 1] = -1.0 * mt_[i, q] * cE[q]
                    fE[q, q] = -1.0 * mt_[i, q] * dcdxE[1]
                    g[q] = 0.0 + (dE[q, q] * dcdxE[q] + fE[q, q] * cE[q])
                g[N - 1] = 0.0 - c[i, 1]
                rj[N - 1, 1] = 1.0
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
                charge_row(g, rj, cc)
                for q in range(2, N):
                    dW[q, q] = -1 * dt_[i - 1, q]
                    dE[q, q] = -1 * dt_[i + 1, q]
                    fW[q, q] = -1.0 * mt_[i - 1, q] * dcdxW[1]
                    fE[q, q] = -1 * mt_[i + 1, q] * dcdxE[1]
                    dW[q, 1] = -1 * mt_[i - 1, q] * cW[q]
                    dE[q, 1] = -1 * mt_[i + 1, q] * cE[q]
                    g[q] = 0.0 - (fW[q, q] * cW[q] + dW[q, q] * dcdxW[q]) + (fE[q, q] * cE[q] + dE[q, q] * dcdxE[q])
            elif j == NJ:
                dW[0, 0] = -(1.0 - por[i]) * v["sigma"]
                g[0] = (self.c_density) - dW[0, 0] * dcdxW[0]
                charge_row(g, rj, cc)
                for q in range(2, N):
                    dW[q, q] = -1.0 * dt_[i, q]
                    dW[q, 1] = -1.0 * mt_[i, q] * cW[q]
                    fW[q, q] = -1.0 * mt_[i, q] * dcdxW[1]
                    g[q] = 0.0 - (dW[q, q] * dcdxW[q] + fW[q, q] * cW[q])
                A[i] = (1.0 - alphaW) * fW - betaW * dW
                B[i] = rj + betaW * dW + alphaW * fW
                G[i] = g
                continue
            elif j < S:
                p = por[i]
                dE[0, 0] = -(1.0 - p) * k.sigma_sep
                dW[0, 0] = -(1.0 - p) * k.sigma_sep
                g[0] = 0.0 - (fW[0, 0] * cW[0] + dW[0, 0] * dcdxW[0]) + (fE[0, 0] * cE[0] + dE[0, 0] * dcdxE[0])
                charge_row(g, rj, cc)
                for q in range(2, N):
                    dW[q, q] = -1 * dt_[i, q]
                    dE[q, q] = -1 * dt_[i, q]
                    fW[q, q] = -1 * mt_[i, q] * dcdxW[1]
                    fE[q, q] = -1 * mt_[i, q] * dcdxE[1]
                    dW[q, 1] = -1 * mt_[i, q] * cW[q]
                    dE[q, 1] = -1 * mt_[i, q] * cE[q]
                    g[q] = - (fW[q, q] * cW[q] + dW[q, q] * dcdxW[q]) + (fE[q, q] * cE[q] + dE[q, q] * dcdxE[q])
                for qe in range(2, N):
                    for qv in range(N):
                        rj[qe, qv] = -(p / self.delT) * delx[i] if qe == qv else 0.0
            else:                                  # cathode interior
                rxn = self.react_tot(*cc, j)
                drx = self.drx_dc(*cc, j)
                # POROSITY_UPDATE_ON: face values interpolated
                por_W = alphaW * por[i] + (1.0 - alphaW) * por[i - 1]
                por_E = alphaE * por[i + 1] + (1.0 - alphaE) * por[i]
                mig_W = alphaW * mt_[i] + (1.0 - alphaW) * mt_[i - 1]
                mig_E = alphaE * mt_[i + 1] + (1.0 - alphaE) * mt_[i]
                dif_W = alphaW * dt_[i] + (1.0 - alphaW) * dt_[i - 1]
                dif_E = alphaE * dt_[i + 1] + (1.0 - alphaE) * dt_[i]
                dE[0, 0] = -(1.0 - por_E) * v["sigma"]
                dW[0, 0] = -(1.0 - por_W) * v["sigma"]
                for q in range(N):
                    rj[0, q] = -drx[q][0] * delx[i]
                g[0] = rxn[3] * delx[i] - (fW[0, 0] * cW[0] + dW[0, 0] * dcdxW[0]) + (fE[0, 0] * cE[0] + dE[0, 0] * dcdxE[0])
                charge_row(g, rj, cc)
                for q in range(2, N):
                    dW[q, q] = -1 * dif_W[q]
                    dE[q, q] = -1 * dif_E[q]
                    fW[q, q] = -1 * mig_W[q] * dcdxW[1]
                    fE[q, q] = -1 * mig_E[q] * dcdxE[1]
                    dW[q, 1] = -1 * mig_W[q] * cW[q]
                    dE[q, 1] = -1 * mig_E[q] * cE[q]
                    g[q] = -(rxn[q + 2]) * delx[i] - (fW[q, q] * cW[q] + dW[q, q] * dcdxW[q]) + (fE[q, q] * cE[q] + dE[q, q] * dcdxE[q])
                p = por[i]
                for qe in range(2, N):
                    for qv in range(N):
                        if qe == qv:
                            rj[qe, qv] = drx[qv][qe - 1] * delx[i] - (p / self.delT) * delx[i]
                        else:
                            rj[qe, qv] = drx[qv][qe - 1] * delx[i]
            A[i] = (1.0 - alphaW) * fW - betaW * dW
            B[i] = rj + betaW * dW + alphaW * fW - (1.0 - alphaE) * fE + betaE * dE
            D[i] = -alphaE * fE - betaE * dE
            G[i] = g
        self.last_coeffs = (dE, dW, fE, fW, rj, g)
        return A, B, D, G

    # ------------------------------------------------------------------ after the solve
    def update_band_variables(self):
        c = self.c + self.delC
        for i in range(self.NJ):
            for q in range(2, self.N):
                x = c[i, q]
                c[i, q] = x if x >= CONC_FLOOR else (x if x != x else CONC_FLOOR)
        self.c = c

    def update_other_variables(self, it: int):
        k, v, s = self.k, self.v, self.s
        NJ, S = self.NJ, self.S
        M = v["molar_mass_ZMC_max"]
        xc = k.xmax_c
        cutoff_theta = f32(0.985)
        c, dC = self.c, self.delC
        dt = self.delT
        for j in range(S + 1, NJ):
            i = j - 1
            r = self.react_tot(*(float(c[i, q] - (dC[i, q] / 2)) for q in range(6)), j)
            s["KMn8O16"][i] = s["KMn8O16"][i] + (k.molar_mass_KMn8O16 * r[8] * dt)
            s["ZHS"][i] = s["ZHS"][i] + (k.molar_mass_ZHS * r[9] * dt)
            s["ZMC_max"][i] = s["ZMC_max"][i] + (M * r[10] * dt)
            s["MnO2"][i] = s["MnO2"][i] + (k.molar_mass_MnO2 * r[12] * dt)
            s["znm"][i] = s["znm"][i] + (r[11] * dt)
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
            for key, floor in (("ZHS", 1.0e-50), ("ZMCx", 1.0e-50), ("ZMC_max", 1.0e-50), ("MnO2", 1.0e-50),
                               ("znm", 1.7e-50), ("mnm", 1.0e-50)):
                if s[key][i] < 0.0:
                    s[key][i] = floor
            s["a_K"][i] = 3.0 * s["KMn8O16"][i] / (k.density_KMn8O16 * xc)
            s["a_ZMCx"][i] = 3.0 * s["ZMCx"][i] / (k.density_ZMC * xc)
            s["a_ZMC_max"][i] = 3.0 * s["ZMC_max"][i] / (k.density_ZMC * xc)
            s["a_ZHS"][i] = 3.0 * s["ZHS"][i] / (k.density_ZHS * xc)
            s["a_MnO2"][i] = 3.0 * s["MnO2"][i] / (k.density_MnO2 * xc)
            por_past = s["por"][i]
            s["por"][i] = (1.0 - v["volfrac_inert"] - (s["KMn8O16"][i] / k.density_KMn8O16)
                           - (s["ZMCx"][i] / k.density_ZMC) - (s["ZHS"][i] / k.density_ZHS)
                           - (s["ZMC_max"][i] / k.density_ZMC) - (s["MnO2"][i] / k.density_MnO2))
            s["tort"][i] = 2.0 * _pow(s["por"][i], -0.5)
            self._transport_terms((j,))
            for q in (2, 3, 4, 5):
                c[i, q] = c[i, q] * por_past / s["por"][i]
            s["pH_phreeqc"][i] = float(self.eval_pH(F32(c[i, 2]), F32(c[i, 3]), j))
            s["pH_ZHS"][i] = self.zhs_pH(float(c[i, 2]), float(c[i, 4]))
            s["pH_standard"][i] = float(F32(-1.0) * _log10f(F32(1000 * c[i, 5])))
        # copies to the boundary nodes; the pH copies use the leftover loop index j = NJ (M-16)
        keys = ("KMn8O16", "ZMCx", "ZMC_max", "ZHS", "MnO2", "znm", "mnm", "ratio", "MW", "por", "tort",
                "a_K", "a_ZMCx", "a_ZMC_max", "a_ZHS", "a_MnO2")
        for key in keys:
            s[key][S - 1] = s[key][S]
        for key in ("pH_phreeqc", "pH_ZHS", "pH_standard"):
            s[key][NJ - 1] = s[key][S]
        for key in keys:
            s[key][NJ - 1] = s[key][NJ - 2]
        for key in ("pH_phreeqc", "pH_ZHS", "pH_standard"):
            s[key][NJ - 1] = s[key][NJ - 2]
        self.anode_pot = self.zn_anode_pot(float(c[0, 2]))

    # ------------------------------------------------------------------ output
    def write_row(self):
        NJ, S = self.NJ, self.S
        n = NJ - S - 1
        idx = range(S, NJ - 1)                    # Fortran SEP_NODE+1 .. NJ-1
        r5_all = [[0.0] * 4 for _ in range(NJ)]
        for i in idx:
            r5_all[i] = self.reaction_5(*(float(x) for x in self.c[i]), i + 1)[:4]

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
        zero = sum_div([0.0 for _ in idx])
        r5 = [sum_div([r5_all[i][q] for i in idx]) for q in range(4)]
        dom = 1                                   # maxloc over |R1..R4 currents|, all zero
        items = [self.time, float(c[NJ - 1, 0]), self.current, self.mAhg, self.state,
                 mean_div(c[:, 2]), mean_div(c[:, 3]), mean_div(c[:, 4]), mean_div(c[:, 5]),
                 mean_div(s["pH_phreeqc"]), mean_div(s["pH_ZHS"]), mean_div(s["pH_standard"]),
                 mean_div(s["KMn8O16"]), mean_div(s["ZMCx"]), mean_div(s["ZMC_max"]), mean_div(s["ZHS"]),
                 mean_div(s["MnO2"]), mean_div(s["ratio"])]
        for q in range(4):
            items += [zero, zero, zero, zero, r5[q]]
        items += [sum_div([s[key][i] for i in idx]) for key in ("a_K", "a_ZMCx", "a_ZMC_max", "a_ZHS", "a_MnO2")]
        items += [self.mAhg * self.v["AM_Grams"] / self.v["AM_Grams_sim"], dom]
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
            if c[NJ - 1, 0] <= 1.0:
                return "EXIT BECAUSE LOWER VOLTAGE CUTOFF"
            if c[NJ - 1, 0] >= 2.0:
                return "EXIT BECAUSE Upper VOLTAGE CUTOFF"
            self.current_ramp()
            if self.current > 1.0e-10:
                self.state = "D"
            elif self.current <= -1.0e-10:
                self.state = "C"
            else:
                self.state = "R"
            if self.mAhg >= 89.0:
                self.delT = f32(0.1)              # single-precision literal
            if self.mAhg >= 100.0:
                self.delT = 1.0
            if it >= self.k.ramp_iters:
                if (self.time - self.k.write_time_density) >= self.last_write_time:
                    self.write_row()
                    self.last_write_time = int(self.time)
            elif it == 1:
                self.write_row()
                self.last_write_time = int(self.time)
            A, B, D, G = self.assemble()
            self.delC = solver.solve(A, B, D, G)
            self.update_band_variables()
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
