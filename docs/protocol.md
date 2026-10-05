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
