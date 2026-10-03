# Validation

What is checked, how, and what is not. The tests are in `tests/`; run them with `python -m pytest` (RuntimeWarnings are errors).

## Faithful ports: byte-identical to the original programs

| Line | Reference | Result |
|---|---|---|
| charge line (`faithful/charge.py`) | the original program, rebuilt from its source and run (gfortran, -O0), five parameter sets | output files **byte-identical** (618, 463, 705, 614 and 252 lines; exits on NaN or the cutoff, as the original) |
| pH-cell line (`faithful/phcell.py`) | the original program's archived outputs, four parameter sets | output files **byte-identical** (1491, 1681, 1492 and 1676 lines) |

**The solver.** Byte-identity needs the original BAND/MATINV operation by operation. That code follows the listing in Newman's *Electrochemical Systems* (Appendix C) and is not distributed: the faithful ports use a private transliteration when it is available (`ZNMNO2_ORACLE`), and otherwise bandsolver's legacy pivot. With the public solver, full runs give the same rows and exits as the originals, with values within 5×10⁻¹¹ (charge line) and 10⁻⁹ (pH-cell line, about 20,000 steps) relative. The difference comes only from fused multiply-adds in the compiled bandsolver: built with `-ffp-contract=off`, bandsolver's legacy pivot (any backend or kernel) reproduces both lines byte for byte (checked on the first 1100 steps of each). The distributed bandsolver is left as it is (contraction off costs 1–15 % in solve time); for byte-identical faithful output, use the Fortran or C++ program, which are built without contraction.

The references (the original sources, their parameter files and outputs) are not distributed. `tests/test_faithful_oracle.py` compares the first 1100 steps of two runs of each line when `ZNMNO2_ORACLE` points to them, and is skipped otherwise. Exactness needed single-precision constants where the original had them, the original's single-precision spline evaluation, IEEE semantics for overflow and NaN, the original BAND/MATINV operation order (a compiled solver differs at 10⁻¹⁴ per solve because of fused multiply-adds), and the original's integer write timer.

## Corrected model

**Speciation** (`tests/test_eqchem.py`):
- agrees with the general solver `speciation.py` (ideal activities) to 10⁻⁹ in pH and 10⁻⁸ relative in every species, on salt and base solutions; `speciation.py` itself reproduces PHREEQC (same equilibria and activity model) to 1.3×10⁻⁵ pH for ionic strengths up to 1 mol/kg;
- every balance closes to 10⁻¹² and the charge of the species equals that of the totals (acid, salt and base solutions);
- reproduces Herrmann et al. (2023), Fig. S1: 2 M ZnSO₄ + 0.5 M MnSO₄ at pH 4 has 72.7 % ZnSO₄⁰, 15.2 % Zn(SO₄)₂²⁻, 12.0 % Zn²⁺;
- sensitivities ∂x/∂T agree with central differences to 10⁻⁵.

**Conservation and balances** (`tests/test_model.py`, `tests/test_options.py`, `tests/test_cycling.py`):
- Mn, S and H conserved to 10⁻¹² relative over a discharge, for each of 18 option configurations; Zn changes by the anode charge / 2F to 10⁻⁸;
- over three discharge/charge cycles: Mn, S and H to 10⁻¹¹, Zn to the charge passed;
- the reaction currents sum to the applied current; electroneutrality to 10⁻¹⁵.

**Jacobian:** the assembled blocks agree with central differences of the whole residual (10⁻⁴ relative to each row's scale) for seven option configurations, including the quasi-particle transport.

**Limits and special cases:**
- a piecewise-linear profile with continuous flux has zero divergence across the probe/separator porosity step (exact);
- `fixed` and `zhs_equilibrium` pH modes give the stated pH; `spline` agrees with the original surrogate to 6×10⁻⁷ pH;
- `no_H` equals `with_H` (10⁻¹⁰ V) when no reaction consumes H⁺;
- `quasi` equals `ions` (10⁻⁹) when every complex is suppressed;
- the Nernstian OCP and the spline-with-Nernstian-ends OCP have the stated values and slopes;
- equilibrium ZHS: w ≤ 1 everywhere, and w = 1 where ZHS exists.

**Rest and equilibrium:** after a partial discharge and a long rest, the electrolyte is uniform to 10⁻⁶ and the net reaction rate in every cell vanishes. This test found an earlier rate law that let current circulate at rest.

**Protocols:** constant-current/constant-voltage charge (voltage held to 10⁻⁶ V, current decaying to `Imin`); a constant-voltage step never accepts a state off the set voltage, and a hold that cannot be held for a whole time step proceeds in sub-steps (`tests/test_driver.py`, with stand-in steppers, and a 1.9 V hold after a discharge in all three implementations); a discharge with acid is not stopped by the upper cutoff; a physical limit ends the step, not the protocol; a discharge followed by a charge recovers the capacity.

## Fortran and C++ programs (`tests/test_ports.py`)

| Check | Cases | Result |
|---|---|---|
| Fortran and C++ corrected model, against each other | every option configuration (18), CC-CV, three cycles, a long rest, a physical limit, GITT, a PHREEQC-format database (with and without quasi-particle transport), log K overrides, the full default mesh | output tables **byte-identical** |
| Fortran (and C++) against Python, corrected model | the same cases | equal to the printed 8 digits, except values near zero (e.g. a volume fraction of 10⁻²³) and the voltage on the row written at a physical limit (below) |
| faithful ports, against the original programs | every reference run (5 charge-line, 4 pH-cell), start to finish | output files **byte-identical**, Fortran and C++ (private, `ZNMNO2_ORACLE`) |
| faithful ports, against the Python faithful ports | the first 1,100 steps of each line | to 10⁻⁶ relative (the installed bandsolver uses fused multiply-adds) |
| parameters | every field | the same names, defaults and namelist groups in the three implementations |

The test suite runs shortened versions of the corrected-model cases (a 30-minute discharge and a rest, and a CC-CV cycle).

**How the ports stay identical.** C++ mirrors the Fortran operation by operation, and both are built without fused multiply-adds (`-ffp-contract=off`, bandsolver included). The C++ BAND solve mirrors bandsolver's partial-pivot kernel exactly, including its zero second-neighbour terms (they can flip the sign of a zero). `pow()` is called with its arguments hidden from the optimizer: clang rewrites `pow(x, 0.5)` as `sqrt(x)` and `pow(10, y)` as `exp10(y)` without fast-math, which changes the last bit. The faithful C++ output was checked with clang at -O0, -O2 and -O3 and with g++ 14 at -O2 and -O3.

**The voltage at a physical limit.** When a step fails because the host is full (`insertion_full`), the state written is the last converged sub-step. There the insertion rate no longer depends on 1 − θ (both terms of the rate vanish with it), so the Newton tolerance does not determine 1 − θ, and the voltage, which depends on ln(1 − θ), differs between implementations (by 87 mV in the test case). Everything else on that row, and every later row, agrees.

## Convergence (`scripts/convergence.py`)

Default discharge (100 mA/g to 1.0 V), voltage at 25, 50, 100, 140 and 170 mAh/g:

- **Time:** first order (backward Euler). At 170 mAh/g the error halves with Δt: 3.9×10⁻⁵, 1.9×10⁻⁵, 9.7×10⁻⁶ V for Δt = 10, 5, 2.5 s. The default Δt = 10 s is within about 8×10⁻⁵ V of the converged curve; the pH within 10⁻⁵.
- **Space:** observed first order, 2.0×10⁻⁶ → 1.0×10⁻⁶ V from mesh ×1 to ×2 to ×4. This comes entirely from the anode boundary, whose kinetics use the first cell's centre values: refining only the probe region changes the voltage by 3.1×10⁻⁶ V, refining the separator and cathode by 7×10⁻⁸ V (the interior scheme is second order).
- The capacity at the cutoff is the same for every mesh and time step (the cutoff is located to 0.1 mV).

## Not checked

- No comparison with measured cells here: the repository holds no data, and the default parameters are illustrative.
- Non-ideal activities (not implemented).
- The quasi-particle transport against an independent implementation; only its limits, its Jacobian and its conservation are tested.
- The Fortran and C++ programs with compilers other than gfortran 14 and Apple clang / g++ 14, or on Windows.
