# Experimental HELD phonons for the Ga2O3 (001) bilayer

The installed unified HELD implementation officially supports
single-element BCC, FCC, and HCP trajectories. It rejects this Ga/O
bilayer through its normal CLI.

`../run_held_phonons.py` is a local experimental adapter. It reuses
HELD's symmetry-constrained force-constant basis, acoustic sum rule,
least-squares design matrix, and dynamical-matrix implementation while
adding:

- species-aware spglib symmetry mapping;
- separate Ga and O masses;
- the custom reciprocal-fractional path `Gamma-X-S-Y-Gamma`.

The external HELD repository was not modified.

## Result from the current trajectory

- Complete MD frames: 20
- Mean sampled temperature: 417.740 K
- Latest complete-frame temperature: 616.859 K
- Species-aware symmetry operations: 1 (identity only)
- HELD coefficients: 144
- Force equations per frame: 30
- Rank per frame: 30
- Aggregation: mean of the per-frame minimum-norm fits, matching the
  standard HELD CLI default
- Force-component RMSE: 1.139 eV/Angstrom
- Force-vector RMSE: 1.973 eV/Angstrom
- Frequency range: -23.090 to 20.274 THz

Each individual fit is strongly underdetermined because it has 30
equations for 144 coefficients. A global fit is formally full rank over
all 20 frames, but its condition number is approximately `3.3e6` and it
produces nonphysical frequencies of order thousands of THz. The global
result is therefore not used.

This result must not be treated as a validated HELD prediction. The
trajectory is short and heats far above 300 K, the reference structure
has substantial residual forces, and the 10-atom primitive simulation
cell is too small for the 2.5 Angstrom interaction cutoff. Use the
80-atom `4 x 2 x 1` trajectory for a meaningful follow-up.

## Output files

- `Ga2O3-BL-001-HELD-experimental.png`: experimental dispersion plot
- `held_dispersion.dat`: path coordinate and 30 branches in THz
- `held_coefficients.csv`: mean HELD coefficients
- `held_step_coefficients.csv`: per-frame minimum-norm coefficients
- `held_summary.json`: machine-readable fit diagnostics
