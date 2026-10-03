# Zn/MnO₂ model: equations and discretization

This document describes the **corrected model** (`python/znmno2_model/model.py`). Its options are listed in [options.md](options.md). The two **faithful** ports, which reproduce the original research programs, are summarized in §9; their defects and how the corrected model treats each one are in [deviations.md](deviations.md).

Units: cm, s, mol, A, V. Electrolyte concentrations are inputs in mol/L and unknowns in mol/cm³. Potentials are referred to the Zn anode metal (0 V), and equilibrium potentials to Zn/Zn²⁺ at 1 M.

## 1. Cell and domain

```
x = 0       L_probe              L_probe + L_sep                    L_probe + L_sep + L_cath
 |-- probe region --|---- separator ----|---------- porous MnO2 cathode ----------|
Zn anode          (open electrolyte)                                        current collector
```

- **Probe region**: open electrolyte between the anode and the separator, where a pH probe sits in the experimental cell. Porosity `eps_probe` (1), tortuosity `tau_probe` (1), length `L_probe`. `L_probe = 0` removes it.
- **Separator**: porosity `eps_sep`, tortuosity `tau_factor_sep · ε^bruggeman_sep`.
- **Cathode**: porosity from the solids (§4), tortuosity `tau_factor_cath · ε^bruggeman_cath`.
- One cross-section `A_cell` throughout (one-dimensional). Every length and porosity is an input.
- The electrolyte starts uniform with the given ZnSO₄, MnSO₄ and H₂SO₄; a salt given as 0 starts as a 10⁻⁹ M trace (the equilibrium potentials need ln c).
- **Mesh**: cell-centred finite volumes, `n_probe`, `n_sep`, `n_cath` uniform cells per region.

## 2. Unknowns (per cell, N = 12)

| Column | Symbol | Meaning | Unit |
|---|---|---|---|
| 0 | φ₁ | solid potential (cathode; set to 0 elsewhere) | V |
| 1 | φ₂ | electrolyte potential | V |
| 2–5 | Zn_T, Mn_T, S_T, H_T | totals of the master species Zn²⁺, Mn²⁺, SO₄²⁻, H⁺ in the electrolyte | mol/cm³ |
| 6 | n_MnO₂ | pristine MnO₂ (R1) | mol/cm³ of electrode |
| 7 | n_ZMO | Zn_zMnO₂, the dissolving/depositing phase (R2) | mol Mn/cm³ |
| 8 | s | insertion log-odds, s = ln(θ/(1−θ)) (R3) | – |
| 9 | n_ZHS | zinc hydroxide sulfate, Zn₄SO₄(OH)₆·5H₂O | mol/cm³ |
| 10, 11 | n_ZnO, n_Zn(OH)₂ | optional precipitates | mol/cm³ |

H_T counts master protons (H⁺, HSO₄⁻, 2·H₂SO₄, minus OH⁻ and the hydroxo complexes). It can be negative.

## 3. Speciation

At every cell the free concentrations of H⁺, Zn²⁺, Mn²⁺ and SO₄²⁻ (and so the pH) follow from the four totals by mass action with ideal activities (concentrations in mol/L):

- species s: `c_s = K_s ∏_m c_m^ν_sm` over the masters m;
- balances: `Σ_s ν_sm c_s = T_m` for m = H, Zn, Mn, S.

The equilibria are the 21 of Herrmann et al. (2023), SI Table S1 (`speciation.HERRMANN_2023`): water, HSO₄⁻, H₂SO₄, seven Zn hydroxo, four Zn sulfate, six Mn hydroxo complexes and MnSO₄⁰. Any log K can be overridden from a file, or the equilibria read from a PHREEQC database ([options.md](options.md)). The solver (`eqchem.py`) is Newton in log₁₀ of the free concentrations, warm-started from the previous state; where Newton fails it brackets the pH (the proton balance is monotonic in it) with the other balances solved at each trial pH. Its sensitivities ∂c/∂T come from the implicit-function theorem.

## 4. Solids and porosity

Cathode porosity: `ε = 1 − ε_fixed − Σ_k V_k n_k` over MnO₂, ZMO, ZHS, ZnO and Zn(OH)₂, with molar volumes `V_k = M_k/ρ_k`. `ε_fixed` (inert material plus the insertion host) is fixed by the initial porosity and volume fractions.

Specific areas of spheres: `a_k = 3 V_k n_k / r_k`. The host area is fixed: `a_host = 3 vf_host / r_host`. The host holds `n_host = vf_host ρ_host / M_host` mol Mn/cm³; its inserted Zn is `n_host (z_min + θ (z_max − z_min))`.

## 5. Reactions

Currents are per electrode volume (A/cm³), positive anodic; `f = F/RT`. `c` denotes the free concentrations from §3 (or the totals, option `basis = totals`); `c_H` follows option `ph_mode`.

**R1, pristine MnO₂ dissolution:** MnO₂ + 4H⁺ + 2e⁻ → Mn²⁺ + 2H₂O, irreversible.
- `U₁ = U1 − (1/2f)(ln c_Mn − 4 ln c_H)`, `η₁ = φ₁ − φ₂ − U₁`;
- `i₁ = −a_MnO₂ F k1 g(−b₁)` with `b₁ = exp(2α₁fη₁) − exp(−2(1−α₁)fη₁)`.

**R2, ZMO dissolution and deposition:** Zn_zMnO₂ + 4H⁺ + (2−2z)e⁻ ⇌ zZn²⁺ + Mn²⁺ + 2H₂O, with z = `z_ZMO` and n₂ = 2 − 2z.
- `U₂ = U2 − (1/(n₂f))(z ln c_Zn + ln c_Mn − 4 ln c_H)`;
- `i₂ = F k2 [a_ZMO b₂ + (a_ZHS + a_seed_R2) g(b₂)]`, with `b₂ = exp(α₂n₂fη₂) − exp(−(1−α₂)n₂fη₂)`.

Both directions occur on ZMO. Deposition can also grow on ZHS and on a seed area (carbon, current collector), but only while it is thermodynamically favoured.

**R3, Zn insertion:** Zn²⁺ + 2e⁻ + host ⇌ Zn(host), θ ∈ (0, 1).
- `U₃` from option `r3_ocp`:
  - `spline_nernst` (default): the empirical 51-point OCP spline, scaled between `V_at_zmin` and `V_at_zmax`, plus `(1/2f) ln(c_Zn/c_ref3)`, plus Nernstian end terms −(1/2f) ln(θ/(1−θ)) applied smoothly within about `r3_end_width` of θ = 0 and 1;
  - `spline`: the same without the end terms;
  - `nernst`: `U3_nernst + (1/2f)(ln c_Zn − s)`.
- `i₃ = a_host F k3 √(c_Zn θ(1−θ)) [exp(2α₃fη₃) − exp(−2(1−α₃)fη₃)]`.

**Precipitation** of ZHS (4Zn²⁺ + SO₄²⁻ + 6H₂O ⇌ ZHS + 6H⁺), ZnO and Zn(OH)₂ (Zn²⁺ + H₂O ⇌ ZnO + 2H⁺, …):
- `w = (Q/K)^(1/ν_H) = c_H,eq/c_H`, with `Q = c_Zn^ν_Zn c_SO4^ν_S / c_H^ν_H` and `log K` from `logK_ZHS`, `logK_ZnO`, `logK_ZnOH2`;
- `r = k [a (w − 1) + a_seed_ZHS g(w − w_nuc)]` (mol/cm³/s, positive = forms), with `w_nuc` the nucleation supersaturation `zhs_nucleation` expressed in w.
- Other ZHS treatments (option `zhs`): instantaneous equilibrium, lumped into R1/R2 as in Bernard et al. (2024), or off.

**g, the one-sided term:** `g(z) = z²/(z + δ)` for z > 0 and 0 otherwise (δ = 10⁻²). It is C¹ and exactly zero at equilibrium. Every rate therefore vanishes at equilibrium, and no current circulates in a cell at rest (a test checks this).

**Sources:** with `ξ₁ = −i₁/2F`, `ξ₂ = −i₂/(n₂F)`, `ξ₃ = −i₃/2F`:

| Quantity | Source per electrode volume |
|---|---|
| Zn_T | `z ξ₂ − ξ₃ − 4r_ZHS − r_ZnO − r_Zn(OH)₂` |
| Mn_T | `ξ₁ + ξ₂` |
| S_T | `−r_ZHS` |
| H_T | `−4ξ₁ − 4ξ₂ + 6r_ZHS + 2r_ZnO + 2r_Zn(OH)₂` |
| n_MnO₂, n_ZMO | `−ξ₁`, `−ξ₂` |
| inserted Zn | `+ξ₃` |
| precipitates | `+r` |

The sources conserve Zn, Mn, S and H exactly; charge is conserved because electroneutrality holds (§6).

## 6. Balances

- **Species** (k = Zn, Mn, S, H): `∂(ε c_k)/∂t = −∂N_k/∂x + source_k`, with dilute-solution diffusion and migration: `N_k = −(ε D_k/τ)(∂c_k/∂x + z_k c_k f ∂φ₂/∂x)`. Charges are 2, 2, −2 and +1.
- **Electroneutrality:** `2 Zn_T + 2 Mn_T + H_T − 2 S_T = 0`.
- **Solid current** (cathode): `∂/∂x[(1−ε) σ ∂φ₁/∂x] = Σ_k i_k`.
- **Solids:** `∂n/∂t = source`.
- **Option `species = no_H`:** H_T is held at its initial value and not transported, as in the original charge line.
- **Option `transport = quasi`:** each total moves with the sum of its species' Nernst–Planck fluxes (Herrmann & Horstmann 2024) instead of with one ion's diffusivity.

## 7. Boundaries

- **Zn anode (x = 0):** Butler–Volmer, `i_an = F k_an √c_Zn [exp(2α f η) − exp(−2(1−α) f η)]` with `η = 0 − φ₂ − (1/2f) ln c_Zn`, using the first cell's values. This is first order in the cell width: 3 µV at the default mesh (see [validation.md](validation.md)). Zn²⁺ enters at `i_an A/(2F)`; the other species have no flux. The anode metal is the potential reference, so the cell voltage is φ₁ at the collector face.
- **Current collector:** the solid current `I` (total, A) leaves; no species flux.
- **Separator/cathode face:** no solid current.

## 8. Discretization and solution

- **Fluxes:** each face has the geometric conductance `G = 1/(h_L τ_L/(2ε_L A) + h_R τ_R/(2ε_R A))`, exact for diffusion across a porosity or tortuosity step. Diffusion plus migration use exponential fitting (Scharfetter–Gummel): `flow = D G [B(a) c_L − B(−a) c_R]`, with `a = z f (φ₂,R − φ₂,L)` and `B(a) = a/(eᵃ − 1)`.
- **Time:** backward Euler. Each step is solved by Newton's method with bandsolver's block-tridiagonal solver (partial pivoting), rows and columns equilibrated.
- **Jacobian:** analytic for storage, electroneutrality and transport (including the quasi-particle fluxes, through the speciation sensitivities). The local reaction and speciation terms use finite differences, with steps relative to each value, so a depleted species keeps its precision.
- **Insertion unknown:** stored as s, solved for in θ where θ > ½ and in s below. θ and 1 − θ are kept separately to full precision.
- **Convergence:** a full Newton step whose update is below `newton_tol` (scaled by typical sizes, the insertion unknown by its change in θ), with every residual below 10⁻⁶ of the 1 mA/g equivalent of its row or its round-off floor.
- **Initial state:** the given compositions, with potentials settled at open circuit (a 10⁻⁴ s step at zero current).
- **Protocol driver:** halves the sub-step on a Newton failure, locates voltage cutoffs to 10⁻⁴ V, finds constant-voltage currents by bracketing, and reports physical limits (host full or empty, Zn²⁺ or Mn²⁺ depleted, MnO₂ exhausted, pores clogged) ([protocol.md](protocol.md)).

## 9. The faithful ports

`python/znmno2_model/faithful/` reproduces two original programs byte for byte (their output files), defects included:

| Line | Original | Unknowns | Notes |
|---|---|---|---|
| `charge.py` | ZnMn02_v3 (the charge line, used for parameter sweeps) | φ₁, φ₂, Zn²⁺, Mn²⁺, SO₄²⁻ (N = 5) | R2 and R3; pH from a closed-form ZHS solubility; Mn²⁺ treated as uncharged (M-1) |
| `phcell.py` | ZnMn02_v2.1_GITT_no_probe (the pH-cell line) | + H⁺ (N = 6) | R5 only (MnO₂ dissolution with a ZHS correction) |

Both use one linearized solve per time step, explicit solid updates, single-precision constants where the original had them, the PHREEQC pH spline evaluated in single precision, and Newman's BAND/MATINV: byte-for-byte with a private transliteration of the original code, otherwise bandsolver's legacy pivot (the same results to within 10⁻⁹ relative; see [validation.md](validation.md)).

## References

- J. C. Bernard et al., J. Electrochem. Soc. 171, 050502 (2024), doi:10.1149/1945-7111/ad3f55 (the base model).
- J. C. Bernard et al., J. Electrochem. Soc. 172, 030517 (2025), doi:10.1149/1945-7111/adbc25 (GITT; pH treatments).
- J. C. Bernard et al., J. Electrochem. Soc. 172, 100509 (2025) (electrolyte engineering: pH and ZnSO₄).
- N. J. Herrmann, H. Euchner, A. Groß, B. Horstmann, Adv. Energy Mater. 14, 2302553 (2024) (equilibria, quasi-particle transport, ZHS).
- N. J. Herrmann, B. Horstmann, Energy Storage Mater. 70, 103437 (2024) (electrolyte design; precipitates).
