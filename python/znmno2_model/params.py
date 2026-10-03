"""Parameters of the corrected model (units: cm, s, mol, A, V; electrolyte inputs in mol/L).

Every length, porosity and area is an input. Defaults are illustrative values from the original
programs (converted where the formulation changed); they are not fitted to any cell.
"""
from __future__ import annotations

from dataclasses import dataclass, fields, replace


@dataclass(frozen=True)
class Params:
    # --- geometry: Zn anode | probe region | separator | cathode | collector ---
    A_cell: float = 0.178134094     # cross-section of the cell [cm2], the same in every region (1-D)
    L_probe: float = 0.1            # probe region length [cm]; 0 = no probe region
    eps_probe: float = 1.0          # probe region: open electrolyte
    tau_probe: float = 1.0          # probe region tortuosity
    L_sep: float = 0.06             # separator thickness [cm]
    eps_sep: float = 0.9            # separator porosity
    tau_factor_sep: float = 2.0     # tortuosity = tau_factor * eps**bruggeman
    bruggeman_sep: float = -0.5     # separator Bruggeman exponent
    L_cath: float = 0.0218          # cathode thickness [cm]
    eps_cath: float = 0.815         # initial cathode porosity
    tau_factor_cath: float = 2.0    # cathode tortuosity = tau_factor_cath * eps**bruggeman_cath
    bruggeman_cath: float = -0.5    # cathode Bruggeman exponent
    n_probe: int = 20               # finite-volume cells per region
    n_sep: int = 30                 # separator cells
    n_cath: int = 40                # cathode cells
    sigma: float = 0.1              # cathode solid conductivity [S/cm]; effective sigma (1 - eps)
    # --- cathode solids: initial volume fractions, molar masses [g/mol], densities [g/cm3], radii [cm] ---
    vf_MnO2: float = 1.0e-4         # pristine MnO2 (R1)
    vf_ZMO: float = 0.01            # Zn_z MnO2, dissolution/deposition (R2)
    vf_host: float = 0.03           # Zn-insertion host (R3), fixed
    vf_ZHS: float = 0.0             # zinc hydroxide sulfate
    M_MnO2: float = 86.937          # MnO2 molar mass [g/mol]
    rho_MnO2: float = 5.03          # MnO2 density [g/cm3]
    rho_ZMO: float = 5.0            # Zn_z MnO2 density [g/cm3]
    rho_host: float = 5.0           # insertion host density [g/cm3]
    M_ZHS: float = 549.819          # Zn4SO4(OH)6 . 5 H2O
    rho_ZHS: float = 2.67           # ZHS density [g/cm3]
    r_MnO2: float = 2.0e-3          # MnO2 particle radius [cm]
    r_ZMO: float = 2.0e-3           # ZMO particle radius [cm]
    r_host: float = 2.0e-3          # host particle radius [cm]
    r_ZHS: float = 2.0e-3           # ZHS (and ZnO, Zn(OH)2) particle radius [cm]
    z_ZMO: float = 0.5              # Zn per Mn in the dissolving phase
    zmin: float = 0.2               # insertion range of the host (Zn per Mn)
    zmax: float = 0.5               # Zn per Mn of the full host
    theta0: float = 1.0e-3          # initial insertion fraction (0 = zmin, charged)
    mass_AM: float = 0.0            # active mass [g] for mAh/g; 0 = MnO2-equivalent mass of the solids
    # --- electrolyte [mol/L] and transport [cm2/s] ---
    c_ZnSO4: float = 2.0            # initial ZnSO4 [mol/L]
    c_MnSO4: float = 0.05           # initial MnSO4 [mol/L]
    c_H2SO4: float = 0.0            # initial H2SO4 [mol/L]
    D_Zn: float = 7.15e-6           # Zn2+ diffusion coefficient [cm2/s]
    D_Mn: float = 6.88e-6           # Mn2+ diffusion coefficient [cm2/s]
    D_SO4: float = 1.07e-5          # SO4 2- diffusion coefficient [cm2/s]
    D_H: float = 9.0e-5             # H+ diffusion coefficient [cm2/s]
    # --- reactions (potentials vs Zn/Zn2+ at 1 M; free-ion concentrations in mol/L) ---
    U1: float = 1.986               # R1 MnO2 + 4H+ + 2e -> Mn2+ + 2H2O
    k1: float = 1.0e-10             # [mol/cm2/s]
    alpha1: float = 0.5             # R1 transfer coefficient
    U2: float = 2.49                # R2 Zn_z MnO2 + 4H+ + (2-2z)e <-> z Zn2+ + Mn2+ + 2H2O
    k2: float = 1.0e-9              # R2 rate constant [mol/cm2/s]
    alpha2: float = 0.5             # R2 transfer coefficient
    a_seed_R2: float = 10.0         # deposition area besides ZMO and ZHS [cm2/cm3]
    k3: float = 1.0e-9              # R3 insertion
    alpha3: float = 0.5             # R3 transfer coefficient
    V_at_zmin: float = 1.75         # empirical OCP spline end points [V]
    V_at_zmax: float = 1.45         # R3 OCP spline value at z_max [V]
    c_ref3: float = 2.0             # reference Zn2+ concentration of the R3 OCP [mol/L]
    k_an: float = 1.0e-6            # Zn anode, i0 = F k_an sqrt(c_Zn2+)
    alpha_an: float = 0.5           # anode transfer coefficient
    logK_ZHS: float = 28.4          # log10([Zn2+]^4 [SO4 2-] / [H+]^6) at saturation (Herrmann et al.)
    k_ZHS: float = 1.0e-7           # precipitation/dissolution rate constant [mol/cm2/s]
    a_seed_ZHS: float = 10.0        # precipitation area besides existing ZHS [cm2/cm3]
    # --- model options (docs: the option matrix; defaults are the corrected model) ---
    ph_mode: str = "speciation"     # pH in the reactions: speciation | zhs_equilibrium | fixed | spline
    pH_fixed: float = 4.8           # ph_mode = fixed (R2Fixed of Bernard et al. 2025)
    species: str = "with_H"         # with_H: H_T transported | no_H: H_T held at its initial value (charge line)
    transport: str = "ions"         # ions: each total moves with its ion's D | quasi: species fluxes summed
    basis: str = "free"             # Zn, Mn, SO4 in the Nernst and rate terms: free (speciation) | totals
    R1_on: bool = True              # R1 (pristine MnO2 dissolution) on
    R2_on: bool = True              # R2 (ZMO dissolution/deposition) on
    R3_on: bool = True              # R3 (Zn insertion) on
    zhs: str = "kinetic"            # kinetic | equilibrium (instantaneous) | lumped (into R1/R2) | off
    zhs_nucleation: float = 1.0     # Zn supersaturation c_Zn/c_sat needed for growth on the seed area (Herrmann: 1.05)
    ZnO_on: bool = False            # extra precipitates (kinetic, same law as ZHS)
    ZnOH2_on: bool = False          # Zn(OH)2 precipitation on
    logK_ZnO: float = 11.17         # log10([Zn2+]/[H+]^2) at saturation (Herrmann & Horstmann 2024, Table 1)
    logK_ZnOH2: float = 12.45       # log10([Zn2+]/[H+]^2) at Zn(OH)2 saturation
    k_ZnO: float = 1.0e-7           # ZnO precipitation rate constant [mol/cm2/s]
    k_ZnOH2: float = 1.0e-7         # Zn(OH)2 precipitation rate constant [mol/cm2/s]
    M_ZnO: float = 81.38            # ZnO molar mass [g/mol]
    rho_ZnO: float = 5.61           # ZnO density [g/cm3]
    M_ZnOH2: float = 99.42          # Zn(OH)2 molar mass [g/mol]
    rho_ZnOH2: float = 3.05         # Zn(OH)2 density [g/cm3]
    r3_ocp: str = "spline_nernst"   # spline_nernst: empirical OCP with Nernstian ends | spline | nernst
    r3_end_width: float = 0.01      # spline_nernst: the end terms act within about this fraction of theta = 0 and 1
    U3_nernst: float = 1.55         # Herrmann et al. (2024), U_ins,Zn
    logk_file: str = ""             # log K overrides ('name value' per line); empty: literature values
    equilibria_db: str = ""         # PHREEQC database for the equilibria; empty: Herrmann et al. (2023) Table S1
    D_OH: float = 5.27e-5           # species diffusion coefficients for transport = quasi [cm2/s]
    D_HSO4: float = 1.33e-5         # HSO4- diffusion coefficient (transport = quasi) [cm2/s]
    D_complex: float = 5.0e-6       # every other complex (assumed; plan Q-1)
    # --- constants ---
    R: float = 8.314462618          # gas constant [J/mol/K]
    T: float = 298.15               # temperature [K]
    F: float = 96485.33212          # Faraday constant [C/mol]
    # --- operation and numerics ---
    steps: str = "cc I=100 Vmin=1.0"   # protocol (docs/protocol.md); I in mA/g, positive = discharge
    cycles: int = 1                 # number of times the protocol is repeated
    end_on_cutoff: bool = True      # a voltage cutoff ends the protocol (GITT); False: go to the next step
    V_min: float = 1.0              # default lower cutoff [V]
    V_max: float = 1.9              # default upper cutoff [V]
    dt: float = 10.0                # [s]
    write_interval: float = 60.0    # [s]
    newton_tol: float = 1.0e-9      # Newton update tolerance (scaled)
    newton_max_iter: int = 30       # Newton iterations per (sub-)step

    def __post_init__(self):
        choices = dict(ph_mode=("speciation", "zhs_equilibrium", "fixed", "spline"), species=("with_H", "no_H"),
                       transport=("ions", "quasi"), basis=("free", "totals"),
                       zhs=("kinetic", "equilibrium", "lumped", "off"), r3_ocp=("spline_nernst", "spline", "nernst"))
        for name, allowed in choices.items():
            if getattr(self, name) not in allowed:
                raise ValueError(f"{name} must be one of {allowed}, got {getattr(self, name)!r}")

    def with_(self, **changes) -> "Params":
        return replace(self, **changes)

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    # --- derived ---
    @property
    def M_ZMO(self) -> float:
        return 65.38 * self.z_ZMO + 54.938 + 2 * 15.999

    @property
    def M_host(self) -> float:
        return 65.38 * self.zmin + 54.938 + 2 * 15.999

    @property
    def mass(self) -> float:
        """Active mass [g]: mass_AM, or the solids' MnO2-equivalent mass (mol Mn x M_MnO2)."""
        if self.mass_AM > 0:
            return self.mass_AM
        mol_mn = (self.vf_MnO2 * self.rho_MnO2 / self.M_MnO2 + self.vf_ZMO * self.rho_ZMO / self.M_ZMO
                  + self.vf_host * self.rho_host / self.M_host)
        return mol_mn * self.M_MnO2 * self.A_cell * self.L_cath
