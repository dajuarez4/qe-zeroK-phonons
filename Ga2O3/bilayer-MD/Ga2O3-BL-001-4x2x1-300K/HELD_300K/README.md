# Experimental HELD phonons: 80-atom Ga2O3 (001) bilayer

The installed HELD CLI officially supports monoatomic BCC/FCC/HCP
systems. `../run_held_phonons.py` uses the local experimental Ga/O
adapter created for this project, with species-aware symmetry and
separate Ga/O masses. The external HELD repository is unchanged.

## Current result

- Unit-cell atoms: 10
- Simulation-supercell atoms: 80
- Complete configurations: 30
- Simulated time represented: 29.028 fs
- Mean sampled temperature: 518.308 K
- Latest complete-frame temperature: 679.508 K
- HELD coefficients: 144
- Force equations per frame: 240
- Per-frame rank: 144 (full rank)
- Global design rank: 144
- Global condition number: 27.27
- Aggregation: mean of 30 per-frame fits
- Force-component RMSE: approximately 3.219 eV/Angstrom
- Force-vector RMSE: approximately 5.576 eV/Angstrom
- Frequency range: approximately -31.741 to 32.882 THz

The larger supercell gives a full-rank fit. However, this trajectory
heats from about 303 K to roughly 680--700 K in only 29 fs. It is not
an equilibrated 300 K ensemble, and the high force error and negative
branches make this an experimental diagnostic rather than a validated
HELD spectrum.

## Important files

- `Ga2O3-BL-001-4x2x1-HELD-experimental.png`: phonon plot
- `held_dispersion.dat`: 30 HELD branches in THz
- `held_coefficients.csv`: mean HELD coefficients
- `held_step_coefficients.csv`: coefficients for all complete frames
- `held_summary.json`: machine-readable diagnostics
