# Cycling protocols and output

The cell is driven by a **protocol**: a list of steps that run in order, optionally repeated. It is given in the `&protocol` group of the input file:

```fortran
&protocol
  steps = 'cc I=100 Vmin=1.0; rest t=1800; cc I=-100 Vmax=1.85; cv V=1.85 Imin=5'
  cycles = 3, end_on_cutoff = .false.
/
```

## Steps

Steps are separated by `;`. Each is a keyword followed by `key=value` settings (case-insensitive, any order):

| Step | Settings | Ends when |
|---|---|---|
| `cc` | `I` = specific current [mA/g of active material], **positive = discharge, negative = charge** (required); `t` = maximum duration [s]; `Vmin`, `Vmax` [V] (default: the global `V_min`, `V_max`) | a discharge reaches `Vmin`, a charge reaches `Vmax`, or `t` has elapsed |
| `cv` | `V` = held voltage [V] (required); `t` = maximum duration [s]; `Imin` = current magnitude [mA/g] below which the step ends | \|I\| ≤ `Imin`, or `t` has elapsed |
| `rest` | `t` = duration [s] (required) | `t` has elapsed |

A discharge is not stopped by `Vmax`, nor a charge by `Vmin`: with acid in the electrolyte a discharge starts above the usual upper cutoff (a high-voltage plateau). `cycles = n` runs the whole list n times.

- **`end_on_cutoff = .true.`** (default): a voltage cutoff or a physical limit ends the whole protocol. Use this for a single discharge, or for **GITT**, e.g. `steps = 'cc I=50 t=720; rest t=1800'` with `cycles = 1000`, which then runs until the cutoff.
- **`end_on_cutoff = .false.`**: a cutoff or a physical limit ends only the current step, and the protocol continues (cycling).

**Specific current.** mA/g refers to the active mass, `mass_AM`, or by default the MnO₂-equivalent mass of the cathode solids (moles of Mn × M_MnO₂).

## How steps are solved

- **Time step** `dt` (10 s by default). A step's last time step is shortened so the step ends exactly at `t`.
- **Adaptive time step** (`dt_mode = 'adaptive'`; the default is `'fixed'`). The time step is chosen from an estimate of backward Euler's local error: the difference between the new state and the linear predictor through the two states before it, times h/(h + h_prev), in the RMS norm weighted by `atol·scale + rtol·|x|`. The scale is 1 for the potentials [V] and the log-odds of θ; for a concentration or an amount it is that variable's largest magnitude on the mesh, at least 10⁻³ of the initial sulfate concentration, so `atol` is relative to each variable's size and variables near zero are controlled absolutely. A step whose error exceeds 1 is redone with h·max(0.2, 0.9/√err); after an accepted step h grows by up to 2×, bounded by `dt_min` and `dt_max` (0: `write_interval`). Steps are shortened to land on the output times and on the end of a step. Each step of the protocol starts with `dt` and no predictor, since a change of current or voltage is a discontinuity. Newton failures, cutoffs, constant voltage and a full or empty host are handled inside each time step as below. On a 200 mA/g discharge to 1.0 V (coarse mesh), the default tolerances (`rtol = atol = 10⁻⁴`) take 443 steps against 378 fixed 10 s steps, with a third of their voltage error (0.13 mV against 0.36 mV, from fixed 1 s steps); `rtol = atol = 10⁻³` matches the fixed steps' accuracy (0.44 mV) in 164 steps.
- **Newton** at every time step ([model.md](model.md) §8). A Newton failure halves the sub-step (down to 10⁻¹⁰ s); each success doubles it again, up to `dt`. After 200 failures in one time step the step stops.
- **Cutoffs** are located inside a time step by halving the sub-step until the voltage is within 0.1 mV of the cutoff.
- **Constant voltage:** each time step finds the current at which the cell voltage equals `V`, by bracketing and regula falsi around the full Newton solve (|V − V_set| ≤ 10⁻⁹ V). If the bracket collapses without reaching the voltage (a jump in V(I)), the state is accepted only within 10⁻⁶ V of `V`; otherwise the step has failed. If no current holds `V` for a whole time step (e.g. a hold right after a discharge), the step is retried over half the time, down to 10⁻¹⁰ s, and the hold continues from there.
- **Cell voltage:** φ₁ at the collector face against the Zn anode metal (0 V). The anode's Butler–Volmer overpotential and its Nernst potential are part of the solution.

## Exit reasons

| Reason | Meaning |
|---|---|
| `cutoff_low`, `cutoff_high` | a discharge reached `Vmin`, or a charge `Vmax` (single-step protocols, or `end_on_cutoff`) |
| `duration`, `current_limit` | the last step ended on `t`, or on `Imin` (cv) |
| `end_of_protocol` | every step completed |
| `insertion_full`, `insertion_empty` | a time step could not be solved and the insertion host was full (on discharge) or empty (on charge); the step is redone and stopped where θ reaches 1 − 10⁻³ (or 10⁻³), and that state is written (docs/validation.md) |
| `zinc_depleted`, `manganese_depleted` | … and the electrolyte had run out of Zn²⁺ or Mn²⁺ somewhere |
| `dissolvable_mno2_exhausted` | … and no MnO₂ or Zn_zMnO₂ was left |
| `pores_clogged` | … and the cathode porosity had fallen below 10⁻³ |
| `solver_fail` | a time step could not be solved and none of the above applies |
| `nan`, `max_time` | non-finite state; 99 h of simulated time |

A physical limit is judged at the last converged sub-step, and the time and capacity up to it are kept.

## Output

`file` (default `znmno2_out.txt`) is a table with one header line and a row every `write_interval` seconds and at the end of every step:

| Column | Meaning |
|---|---|
| `t_h` | time [h] |
| `V` | cell voltage [V] |
| `I_mAg` | applied current [mA/g], positive on discharge |
| `mAhg` | capacity passed [mAh/g] (discharge positive) |
| `step` | 1-based index of the step in the expanded list (cycles × steps) |
| `pH_cath`, `pH_probe`, `pH_anode` | pH: from the cathode's mean free H⁺, at the middle of the probe region, and in the first cell |
| `Zn_cath_M`, `Mn_cath_M`, `S_cath_M` | mean totals in the cathode [mol/L] |
| `vf_MnO2`, `vf_ZMO`, `vf_ZHS` | mean volume fractions in the cathode |
| `theta` | mean insertion fraction of the host |
| `i_R1_mAg`, `i_R2_mAg`, `i_R3_mAg` | each reaction's current, summed over the cathode [mA/g] (positive = anodic) |
| `eps_cath` | mean cathode porosity |

## Restart files

`restart_write = '<file>'` (group `&output`) writes the state at the end of a run: every cell's unknowns, the reference state the physical limits are measured against (the first run's initial state), the time, the capacity and the last current, in plain text with 17 significant digits (an exact round trip). `restart_read = '<file>'` starts a run from it: the protocol runs from its first step at that state, time and capacity, and a first constant-voltage step starts its current search from the stored current. A run continued this way repeats the run it continues exactly (in Python, and byte for byte in the Fortran and C++ programs; `tests/test_restart.py`); its first output row repeats the last row of that run. The mesh must match the input. The three implementations read each other's files. Corrected mode only; both are empty (off) by default.

## Choosing the output

Both are in `&output`, corrected mode, and off by default (the output is then as above).

- **Columns** (`columns`): the time-series columns to write, by name, comma- or blank-separated, case-insensitive (e.g. `columns = 'V, pH_cath, theta'`). `t_h` is always written; the others keep their order and format. An unknown name stops the run.
- **Profiles** (`profile_file`): the state against position, in one long table: a header row, then one row per cell per profile time with `t_h` [h], `x_um` [µm] and the variables of `Model.profiles`: `pH`, `Zn_M`, `Mn_M`, `S_M`, `H_M` (totals [M]), `Zn_free_M` [M], `phi1`, `phi2` [V], `vf_ZHS`, `vf_ZMO`, `vf_MnO2`, `theta` and `eps`, each ES15.7 as in the time series. `profile_vars` chooses some (default: all). A value that does not exist in a cell (the solid potential and the solids outside the cathode) is NaN. A profile reuses its row's speciation, so it never changes the run.
- **When:** profiles are written on output rows, so they never change the time steps: with every `profile_every`-th row (counting the first as 0), at the first row at or after each time in `profile_times` [s], and with the last row.
