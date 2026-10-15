# Changelog

## Unreleased

- **Configurable output** (all three implementations; #15): `columns` chooses the time-series columns; `profile_file` writes the variables of `Model.profiles` (now with `phi1`) against x in one long table, with every `profile_every`-th output row, at `profile_times` and at the end; `profile_vars` chooses the variables. Off by default: the output is unchanged.
- **Restart files** (all three implementations; #16): `restart_write` writes the final state, the reference state of the limits, time, capacity and current; `restart_read` starts a run from it. A continued run repeats the run it continues exactly; the implementations read each other's files.
- **Fortran program:** an error reading `&run`, `&output` or `&faithful` exited with status 0 (`stop 'message'`); it now prints `error: ...` and exits with status 1, like the model's other errors.
- **CI** (#20): every supported Python version (3.9-3.13), coverage in the run summary; the port tests' input builder refuses names it cannot place.

## 0.3.0 (2026-10-14)

- **Host-to-ZMO transfer** (`zmc_transfer`, off by default; all three implementations; #5): the insertion host's Mn is a new unknown (N = 13), so a nearly full host can turn into ZMO at a finite rate, smooth in θ around z = 0.985 z_max (the original programs' transfer, deviation M-9). The transferred Mn takes z_ZMO Zn along; the host keeps the rest; Zn, Mn, S and H stay conserved to round-off. The default output is unchanged.
- **D_complex:** its comment cited an internal note; it now says the value is assumed, with no cited source yet (#1).
- **Anode boundary at second order** (all three implementations; #6): the Zn anode kinetics use φ₂ and c_Zn²⁺ extrapolated to the anode face instead of the first cell's centre values, with the matching Jacobian terms. The voltage now converges at second order in space (it was first order); the default output changes by up to about 3 µV (more near the end of a discharge).
- **Adaptive time stepping** (all three implementations; #7): `dt_mode = 'adaptive'` chooses each time step from backward Euler's local error estimate (linear predictor), against `rtol` and `atol` (scaled by each variable's size), within `dt_min` and `dt_max`; steps land on the output times. Fortran and C++ take the same steps. The default stays `dt_mode = 'fixed'`, so existing inputs give the same output.

## 0.2.0 (2026-10-05)

- **Driver:** a run that ends at a full (or empty) insertion host writes the state where θ reaches 1 − 10⁻³ (or 10⁻³) instead of the last converged sub-step, whose voltage Newton did not determine (87 mV apart between implementations). Other runs are unchanged (#14).
- **Example inputs** in `input/examples/`: GITT, CC-CV, acid electrolyte with the probe region, and the two faithful lines; tested (#12).
- **Input checks:** the Fortran program rejects an unknown namelist group (it skipped it silently); all three implementations reject a group that appears twice (Fortran read the first, C++ merged them, Python kept the last); the C++ program names an unknown group as such.
- **Tests:** the programs' rejection of bad input (#10); quadratic Newton convergence (#9). `Model.newton_step` takes an optional `history` list.
- **CI:** a weekly workflow (also run on demand) executes every notebook (#13); tests read text files as UTF-8, so they pass on Windows; the version test also checks `CMakeLists.txt` and the newest CHANGELOG heading.

## 0.1.0 (2026-10-03)

First version.

- Corrected model (Python): speciation-coupled pH; R1–R3 with implicit solids; finite-rate ZHS (and optional ZnO, Zn(OH)₂); Zn anode with Butler–Volmer kinetics; probe region; protocol driver (cc, cv, rest, cycling, GITT) with physical exit reasons.
- Options for the model versions of the original research: pH handling, species handling, quasi-particle transport, ZHS treatment, reaction switches, insertion OCP form, log K overrides, PHREEQC databases.
- Faithful ports (Python) of the charge line and the pH-cell line, byte-identical to the original programs.
- Fortran program (with bandsolver) and self-contained C++ program: the corrected model with every option and the protocol driver, and the faithful ports; the same input file as the Python package. Fortran and C++ write byte-identical output; the faithful ports are byte-identical to the original programs.
- Driver: a constant-voltage step accepts a collapsed bracket only within 10⁻⁶ V of the set voltage, and proceeds in sub-steps when no current holds the voltage for a whole time step (learnings from lfp-model).
- CI on Linux, macOS and Windows (MinGW), with lint (ruff).
- Tests: speciation, conservation, Jacobian, limits, cycling and relaxation; the three implementations against each other; convergence study.
