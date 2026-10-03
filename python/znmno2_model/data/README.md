# Data files

Two tables from the original research programs, used by the model as inputs. Plain text, read by the Python, Fortran and C++ implementations; the values are float32 (as in the original programs) written with 9 significant digits, which reproduces them exactly.

## `r3_ocp_spline.txt`: the empirical Zn-insertion open-circuit potential

- `theta` (51 points, 0 to 1) and `coeffs` (50 × 4): a piecewise cubic, V_s(θ) = p₁d³ + p₂d² + p₃d + p₄ with d = θ − θ_k on [θ_k, θ_{k+1}], normalized so that V_s(0) = 1.
- The model scales it as V = (V_at_zmin − V_at_zmax)·V_s + V_at_zmax, so the shape is fixed and the end voltages are inputs.
- Derived from measured cell data in the original research (Bernard et al., J. Electrochem. Soc. 171, 050502, 2024). Copied unchanged from the original program.
- Used by `r3_ocp = spline_nernst` (the default) and `spline`.

## `ph_spline_phreeqc.txt`: the original programs' pH surrogate

- `zn_knots`, `mn_knots` (153 each) and `coeffs` (149 × 149): a bicubic B-spline in log₁₀ of the total Zn and Mn (mol/L), giving pH.
- Fitted (scipy `RectBivariateSpline`) to PHREEQC solutions of ZnSO₄ + MnSO₄ on a 149 × 149 logarithmic grid (10⁻¹⁰ to 10 mol/kg water). The equilibria were those of Herrmann et al. (2023) added to `phreeqc.dat`, with log K values adjusted to operando pH-probe measurements (Bernard et al., J. Electrochem. Soc. 172, 030517, 2025, SI Table S1).
- Copied unchanged from the original program. Note the original evaluates it with its first knot axis at log₁₀(Mn) and the second at log₁₀(Zn); the arguments are swapped on purpose, to undo a transposed coefficient read.
- Used only by `ph_mode = spline` and by the faithful ports. The default model computes the pH from the equilibria instead.
