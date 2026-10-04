# znmno2-model

A one-dimensional model of an aqueous Zn/MnO₂ cell: Zn anode | electrolyte (pH-probe) region | separator | porous MnO₂ cathode | current collector, in a ZnSO₄/MnSO₄ (optionally H₂SO₄) electrolyte. It is solved with [bandsolver](https://github.com/jcbernard87/bandsolver), an implementation of Newman's BAND method.

The cathode reacts in three ways: dissolution of pristine MnO₂ (R1), dissolution and deposition of Zn_zMnO₂ (R2), and Zn insertion (R3). Zinc hydroxide sulfate (ZHS) and, optionally, ZnO and Zn(OH)₂ precipitate and dissolve. The pH and the free ions come from the solution equilibria (21 complexes) at every point and time. The model is discretized with finite volumes and advanced by backward Euler, with Newton's method at every step, and runs constant-current, constant-voltage, rest, cycling and GITT protocols.

This repository shares the **model only**. It contains no simulation results or experimental data; run the model to generate results.

Three implementations share one input file: the Python package (`python/znmno2_model`), a Fortran program (`fortran/`, using bandsolver) and a self-contained C++ program (`cpp/znmno2.cpp`). The Fortran and C++ programs write byte-identical output, and agree with Python to the printed digits; they run 10–15 times faster.

## Quick start

```sh
pip install "bandsolver @ git+https://github.com/jcbernard87/bandsolver@v0.1.2"
pip install -e ".[notebooks]"
python -m znmno2_model input/default.nml      # writes znmno2_out.txt
```

```python
from znmno2_model.params import Params
from znmno2_model.simulate import run

r = run(Params(steps="cc I=100 Vmin=1.0; cc I=-100 Vmax=1.85", end_on_cutoff=False))
capacity, volts = r.column("mAhg"), r.column("V")
```

### Fortran and C++ programs

```sh
cmake -S . -B build          # fetches bandsolver v0.1.2; or -DFETCHCONTENT_SOURCE_DIR_BANDSOLVER=/path/to/bandsolver
cmake --build build
build/fortran/znmno2_f input/default.nml
build/cpp/znmno2_cpp input/default.nml
```

Run them from the repository root, or set `data_dir` in the `&run` group to the folder of the spline tables (`python/znmno2_model/data`). The C++ program also builds on its own: `c++ -std=c++17 -O2 -ffp-contract=off cpp/znmno2.cpp -o znmno2_cpp`. `&run mode` selects the corrected model (default) or a faithful port (`faithful_charge`, `faithful_phcell`, with the original run inputs in a `&faithful` group: `rxnk_2, rxnk_3, frac_zmcx, frac_zmcmax, xmax_t, applied_current, porosity, volfrac_mno2, stated_mass_loading, zhs_ksp` for the charge line, `rxnk_5, fraction_kmno2` for the pH-cell line; the output file is `file` in `&output`).

## Model versions as options

The model went through several versions in the original research: with and without the pH-probe region, different ways of computing the pH, different reaction and ZHS treatments, with and without transported H⁺, and with "quasi-particle" transport of the totals. Each is an option of this model, so their effects can be compared on the same cell. See [docs/options.md](docs/options.md).

## Faithful ports of the original programs

Two of the original research programs are reproduced, defects included, for comparison with earlier work: the charge line and the pH-cell line (`znmno2_model.faithful`, and `mode = 'faithful_charge'` / `'faithful_phcell'` in the Fortran and C++ programs). The Fortran and C++ ports reproduce the originals' output files byte for byte (verified privately on every reference run; the references are not distributed). The Python ports do so with the original BAND/MATINV code (not distributed); with the installed bandsolver they agree to within 10⁻⁹ relative. The defects, and how the corrected model treats each, are in [docs/deviations.md](docs/deviations.md).

## Input

One Fortran-namelist file, [`input/default.nml`](input/default.nml), read by all three implementations, holds every parameter: geometry, solids, electrolyte, reactions, options, protocol, numerics and output. Every length, porosity and area is an input. [docs/parameters.md](docs/parameters.md) lists each one with its default and meaning.

[`input/examples/`](input/examples) has complete inputs for a GITT discharge (`gitt.nml`), a CC-CV cycle (`cccv.nml`), an acid electrolyte with the probe region (`acid_probe.nml`), and a faithful run of each original line (`faithful_charge.nml`, `faithful_phcell.nml`; Fortran and C++ programs).

## Documentation

- [docs/model.md](docs/model.md): equations, speciation, reactions, boundaries, discretization and solution
- [docs/options.md](docs/options.md): the model versions as options
- [docs/parameters.md](docs/parameters.md): every parameter
- [docs/protocol.md](docs/protocol.md): protocol steps, exit reasons and the output columns
- [docs/deviations.md](docs/deviations.md): defects in the original programs and the corrected model's treatment
- [docs/validation.md](docs/validation.md): what is tested and how, convergence, and what is not checked

## References

- J. C. Bernard et al., J. Electrochem. Soc. 171, 050502 (2024), the base model; 172, 030517 (2025), GITT and the pH treatments; 172, 100509 (2025), electrolyte engineering.
- N. J. Herrmann, H. Euchner, A. Groß, B. Horstmann, Adv. Energy Mater. 14, 2302553 (2024); N. J. Herrmann, B. Horstmann, Energy Storage Mater. 70, 103437 (2024).
- J. Newman, K. E. Thomas-Alyea, *Electrochemical Systems*, Appendix C (the BAND method).

## License

BSD 3-Clause (see [LICENSE](LICENSE)).
