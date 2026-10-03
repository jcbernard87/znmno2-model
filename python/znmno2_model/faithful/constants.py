"""Constants of the original charge-line program (ZnMn02_v3), evaluated as gfortran evaluates them.

Template placeholders are substituted as text, so each becomes a default-kind (single-precision)
literal: its value is the float32 rounding of the decimal, and products of two such literals are
folded in single precision. `10.0**(x)` with a single-precision literal is folded to the correctly
rounded float32 result.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def f32(x) -> float:
    """Round to float32 and return a Python float."""
    return float(np.float32(x))


def f32_mul(a, b) -> float:
    return float(np.float32(a) * np.float32(b))


def f32_pow10(x) -> float:
    """10.0**(x) for a single-precision literal x, folded at compile time (correctly rounded)."""
    return f32(10.0 ** f32(x))


@dataclass
class Template:
    """One row of the run generator's parameter file (values as the text that was substituted)."""
    RXNK_2: str
    RXNK_3: str
    FRAC_ZMCX: str
    FRAC_ZMCMAX: str
    XMAX: str
    APPLIED_CURRENT: str
    POROSITY: str
    VOLFRAC_MNO2: str
    STATED_MASS_LOADING: str
    ZHS_KSP: str


@dataclass
class Constants:
    N: int = 5
    NJ: int = 122
    SEP_NODE: int = 51
    # general
    Rigc: float = f32(8.314)
    Temp: float = 298.0
    Fconst: float = 96485.0
    # electrolyte
    diff_Zn: float = 7.15e-6
    diff_Mn: float = 6.88e-6
    diff_SO4: float = 1.07e-5
    pH_init: float = 5.5
    z_Zn: float = 2.0
    z_Mn: float = 2.0
    z_SO4: float = -2.0
    # solids
    molar_mass_KMn8O16: float = f32(734.59)
    density_KMn8O16: float = f32(5.03)
    molar_mass_ZHS: float = f32(549.819)
    density_ZHS: float = f32(2.67)
    density_ZMC: float = 5.0
    xmax_c_KMn8O16: float = 200.0e-5
    xmax_c_ZHS: float = 200.0e-5
    xmax_c_ZMC: float = 200.0e-5
    V_at_Zmin: float = 1.75
    V_at_Zmax: float = f32(1.45)
    BL_thickness: float = 0.5 * 1.0e-4
    sigma_sep: float = 1.0e-20
    ramp_iters: int = 1000
    ramp_delT: float = 1.0e-9
    delT_nominal: float = 1.0
    write_time_density: float = 10.0
    # set from the template in derive()
    values: dict = field(default_factory=dict)

    def derive(self, t: Template) -> "Constants":
        v = self.values
        v["Rxn1_K"] = f32(10.0 ** -10.0)
        v["Rxn2_K"] = f32_pow10(t.RXNK_2)
        v["Rxn3_K"] = f32_pow10(t.RXNK_3)
        v["Rxn4F_K"] = f32(1.0e-6)
        v["Rxn4B_K"] = f32(1.0e-6)
        v["K_sp"] = f32(10.0 ** int(t.ZHS_KSP)) if t.ZHS_KSP.lstrip("-").isdigit() else f32_pow10(t.ZHS_KSP)
        v["Zmin"] = f32(0.2)
        v["Zmax"] = 0.5
        v["ratio_initial"] = v["Zmin"] * f32(1.001)
        v["cbulk_Zn"] = f32(0.002)
        v["cbulk_Mn"] = f32(0.00005)
        v["sigma"] = f32(0.1)
        v["xmax"] = f32(t.XMAX)
        v["Area_CS"] = f32(0.178134094)
        len_sep = 600.0 * 1.0e-4
        v["eps_sep"] = f32(0.9)
        v["current"] = f32(t.APPLIED_CURRENT)
        v["eps"] = f32(t.POROSITY)
        v["tortuosity"] = 2 * v["eps"] ** (-0.5)
        v["volfrac_KMn8O16"] = f32(0.0001)
        v["volfrac_ZHS"] = f32(0.0001)
        v["volfrac_ZMCx"] = f32_mul(t.VOLFRAC_MNO2, t.FRAC_ZMCX)
        v["volfrac_ZMC_max"] = f32_mul(t.VOLFRAC_MNO2, t.FRAC_ZMCMAX)
        v["volfrac_inert"] = 1.0 - v["eps"] - v["volfrac_KMn8O16"] - v["volfrac_ZHS"] - v["volfrac_ZMCx"] - v["volfrac_ZMC_max"]
        v["AM_Grams"] = f32(t.STATED_MASS_LOADING)
        v["AM_Grams_sim"] = v["volfrac_ZMCx"] * v["Area_CS"] * v["xmax"] * self.density_ZMC
        v["KMn8O16_init"] = v["volfrac_KMn8O16"] * self.density_KMn8O16
        v["ZHS_init"] = v["volfrac_ZHS"] * self.density_ZHS
        v["ZMCx_init"] = v["volfrac_ZMCx"] * self.density_ZMC
        v["ZMC_max_init"] = v["volfrac_ZMC_max"] * self.density_ZMC
        v["Phi_1_init"] = f32(1.7)
        # molar_mass_ZMC_max = 65.38*Zmax + 54.93 + (15.999*2): single literals; 15.999*2 folded in single
        v["molar_mass_ZMC_max"] = f32(65.38) * v["Zmax"] + f32(54.93) + f32_mul(15.999, 2)
        # ions 3..5 (Zn, Mn, SO4); z_ion(4) = z_Mn0: an undeclared, implicitly real variable = 0 (M-1)
        z = [0.0, 0.0, self.z_Zn, 0.0, self.z_SO4]
        c = [0.0, 0.0, v["cbulk_Zn"], v["cbulk_Mn"], 0.0]
        c[4] = -1 * ((z[2] * c[2]) + (z[3] * c[3])) / z[4]
        v["z_ion"], v["c_initial"] = z, c
        v["diff_ion"] = [0.0, 0.0, self.diff_Zn, self.diff_Mn, self.diff_SO4]
        v["ramp_current_target"] = v["current"]
        v["resevoir_vol"] = f32(0.0001)
        sc = (v["resevoir_vol"] + (len_sep * v["Area_CS"] * v["eps_sep"]) + (v["xmax"] * v["Area_CS"] * v["eps"])) / (
            len_sep * v["Area_CS"] * v["eps_sep"])
        v["resevoir_scaling"] = sc
        v["len_sep"] = len_sep * sc
        return self
