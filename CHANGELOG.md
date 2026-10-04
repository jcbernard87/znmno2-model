# Changelog

## Unreleased

- **Input checks:** the Fortran program rejects an unknown namelist group (it skipped it silently); all three implementations reject a group that appears twice (Fortran read the first, C++ merged them, Python kept the last); the C++ program names an unknown group as such.
- **Tests:** the programs' rejection of bad input (#10); quadratic Newton convergence (#9). `Model.newton_step` takes an optional `history` list.

## 0.1.0 (2026-10-03)

First version.

- Corrected model (Python): speciation-coupled pH; R1–R3 with implicit solids; finite-rate ZHS (and optional ZnO, Zn(OH)₂); Zn anode with Butler–Volmer kinetics; probe region; protocol driver (cc, cv, rest, cycling, GITT) with physical exit reasons.
- Options for the model versions of the original research: pH handling, species handling, quasi-particle transport, ZHS treatment, reaction switches, insertion OCP form, log K overrides, PHREEQC databases.
- Faithful ports (Python) of the charge line and the pH-cell line, byte-identical to the original programs.
- Fortran program (with bandsolver) and self-contained C++ program: the corrected model with every option and the protocol driver, and the faithful ports; the same input file as the Python package. Fortran and C++ write byte-identical output; the faithful ports are byte-identical to the original programs.
- Driver: a constant-voltage step accepts a collapsed bracket only within 10⁻⁶ V of the set voltage, and proceeds in sub-steps when no current holds the voltage for a whole time step (learnings from lfp-model).
- CI on Linux, macOS and Windows (MinGW), with lint (ruff).
- Tests: speciation, conservation, Jacobian, limits, cycling and relaxation; the three implementations against each other; convergence study.
