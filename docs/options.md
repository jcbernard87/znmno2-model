# Model options

The original research went through several versions of this model: with and without the pH-probe region, different ways of getting the pH, different reaction and ZHS treatments, with and without H⁺ as a transported species. Here each of those differences is an **option** of one corrected model, so the versions can be run side by side and compared on the same cell (the comparison notebooks do this). Defaults are in **bold**. All options are set in the `&options` group of the input file (and the geometry ones in `&cell`); [parameters.md](parameters.md) lists every value.

Every option keeps the model's conservation of Zn, Mn, S and H and its charge balance; the tests check each one ([validation.md](validation.md)).

## Geometry

| Option | Values | Meaning |
|---|---|---|
| `L_probe` | **0.1 cm**, or 0 | The open-electrolyte region between the Zn anode and the separator where the pH probe sits. 0 removes it (the cell without a probe). |

## How the pH is obtained

| `ph_mode` | Meaning | Corresponds to |
|---|---|---|
| **`speciation`** | The pH (and the free Zn²⁺, Mn²⁺, SO₄²⁻) from the 21 solution equilibria at every point, from the transported totals. | Herrmann et al. (2024); the question "what is the pH in the cell at any time, given the solvation structures" |
| `zhs_equilibrium` | The pH at which ZHS is saturated: c_H = (c_Zn⁴ c_SO₄ / K_sp)^(1/6). | *R2Precipitation* of Bernard et al. (2025, GITT); the original charge line's `ZHS_pH` |
| `fixed` | A constant pH, `pH_fixed` (4.8). | *R2Fixed* of Bernard et al. (2025, GITT) |
| `spline` | The original programs' PHREEQC surrogate pH(Zn_T, Mn_T), a bicubic spline. Needs the spline data file. | *R2PHREEQC* of Bernard et al. (2025, GITT); the original pH-cell line |

`ph_mode` sets the H⁺ concentration used in the reaction potentials and rates. The ZHS precipitation rate always uses the solution's own pH from the speciation.

| `basis` | Meaning |
|---|---|
| **`free`** | Zn²⁺, Mn²⁺, SO₄²⁻ in the Nernst and rate terms are the free ions from the speciation. |
| `totals` | They are the totals, as in the original programs. |

## Species handling

| `species` | Meaning | Corresponds to |
|---|---|---|
| **`with_H`** | H⁺ (the master-proton total H_T) is transported and consumed or released by the reactions. | the original pH-cell line |
| `no_H` | H_T is held at its initial value; only Zn²⁺, Mn²⁺ and SO₄²⁻ are transported. | the original charge line (N = 5) |

## Transport

| `transport` | Meaning | Corresponds to |
|---|---|---|
| **`ions`** | Each total moves with one ion's diffusion coefficient and charge (`D_Zn`, `D_Mn`, `D_SO4`, `D_H`), as dilute-solution theory for free ions. | the original programs |
| `quasi` | Each total moves with the sum of the Nernst–Planck fluxes of all the species it is part of (free ions and complexes), each with its own charge and diffusion coefficient (`D_OH`, `D_HSO4`, `D_complex` for the others). The "_T" quasi-particle transport. | Herrmann & Horstmann (2024), Eqs. 12–16 |

## Reactions

| Option | Values | Meaning |
|---|---|---|
| `R1_on`, `R2_on`, `R3_on` | **true** / false | Switch each reaction: R1 pristine MnO₂ dissolution, R2 Zn_zMnO₂ dissolution and deposition, R3 Zn insertion. |
| `zhs` | **`kinetic`** | ZHS forms and dissolves at a finite rate, driven by the local supersaturation. |
| | `equilibrium` | ZHS is always at equilibrium where it exists (instantaneous; a complementarity condition). |
| | `lumped` | ZHS forms in proportion to R1 and R2, so they consume no net H⁺ (Bernard et al. 2024, Eqs. 4b and 5b). Deposition then needs ZHS present (the original's charge-area switch). |
| | `off` | No ZHS. |
| `zhs_nucleation` | **1.0**, e.g. 1.05 | Zn supersaturation needed before ZHS grows on the seed area (Herrmann et al. 2024 use 105 %). |
| `ZnO_on`, `ZnOH2_on` | **false** / true | ZnO and Zn(OH)₂ precipitation, with the same rate law as ZHS. |
| `r3_ocp` | **`spline_nernst`** | The empirical insertion OCP (spline) with Nernstian end terms. |
| | `spline` | The empirical OCP as in the original programs. |
| | `nernst` | A Nernstian insertion OCP, `U3_nernst + RT/2F (ln c_Zn − ln θ/(1−θ))` (Herrmann et al. 2024). |

**Why the insertion OCP needs Nernstian ends.** The same material reacts by insertion (R3) and by dissolution (R2): a discharge goes through one and then the other, and the cell must be able to go back. With an OCP that stays finite at θ = 0 and 1, and an exchange current ∝ √(θ(1−θ)), the host fills or empties *exactly* in a finite time, and then has no exchange current: it can no longer react. A charge after a full discharge then recovers only what R2 deposition can supply. The logarithmic end terms make the host approach full or empty only asymptotically, so both reactions reverse. With `spline_nernst` the middle of the empirical curve is unchanged (within 1 mV for 0.1 ≤ θ ≤ 0.9).

## Chemistry

| Option | Values | Meaning |
|---|---|---|
| `logk_file` | **''** or a path | Override log K values: one `species_name value` pair per line (e.g. fitted values). Empty: the literature values. |
| `equilibria_db` | **''** or a path | Read the equilibria from a PHREEQC database (species made of H⁺, Zn²⁺, Mn²⁺, SO₄²⁻ and H₂O) instead of Herrmann et al. (2023) Table S1. |

Activities are ideal (concentrations). Non-ideal activity models (PHREEQC's WATEQ/Davies, Pitzer) are not implemented: Davies-type models are outside their range at the cell's ionic strength (about 8 mol/L).

## Electrolyte

`c_ZnSO4`, `c_MnSO4` and `c_H2SO4` set the initial electrolyte (mol/L). Acid adds a high-voltage first plateau (R2 is driven by the low pH).
