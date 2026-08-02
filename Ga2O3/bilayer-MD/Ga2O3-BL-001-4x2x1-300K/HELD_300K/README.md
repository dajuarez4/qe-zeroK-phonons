# Experimental HELD phonons: 80-atom Ga2O3 (001) bilayer

The installed HELD CLI officially supports monoatomic BCC/FCC/HCP
systems. `../run_held_phonons.py` uses the local experimental Ga/O
adapter created for this project, with species-aware symmetry and
separate Ga/O masses. The external HELD repository is unchanged.

## Current result

- Unit-cell atoms: 10
- Simulation-supercell atoms: 80
- Complete configurations: 126
- Simulated time represented: 121.918 fs
- Mean sampled temperature: 532.568 K
- Latest complete-frame temperature: 584.624 K
- HELD coefficients: 144
- Force equations per frame: 240
- Per-frame rank: 144 (full rank)
- Global design rank: 144
- Global condition number: 47.95
- Aggregation: mean of 126 per-frame fits
- Force-component RMSE: approximately 3.572 eV/Angstrom
- Force-vector RMSE: approximately 6.186 eV/Angstrom
- Frequency range: approximately -14.616 to 23.701 THz

The larger supercell gives a full-rank fit. However, this trajectory
remains far above 300 K and spans only 122 fs. It is not
an equilibrated 300 K ensemble, and the high force error and negative
branches make this an experimental diagnostic rather than a validated
HELD spectrum.

## Important files

- `Ga2O3-BL-001-4x2x1-HELD-experimental.png`: phonon plot
- `Ga2O3-BL-001-4x2x1-HELD-step-dashboard.gif`: 126-frame
  per-step HELD dashboard animation at 8 fps
- `Ga2O3-BL-001-4x2x1-HELD-step-dashboard-last-frame.png`: static
  preview of the final dashboard frame
- `held_step_dashboard_summary.json`: animation and per-step diagnostics
- `held_step_dashboard_cache.npz`: cached 126 per-step dispersions
- `held_dispersion.dat`: 30 HELD branches in THz
- `held_coefficients.csv`: mean HELD coefficients
- `held_step_coefficients.csv`: coefficients for all complete frames
- `held_summary.json`: machine-readable diagnostics

The dashboard force RMSE is evaluated against the same individual
frame used to fit each coefficient row. Its mean value (about
0.198 eV/Angstrom) is therefore an in-sample per-frame diagnostic and
must not be confused with the 3.572 eV/Angstrom residual obtained when
the mean coefficients are applied to the complete trajectory.
