# Experimental HELD diagnostic for the 4 Angstrom bilayer

This directory contains an experimental multi-species Ga/O HELD fit using the
594 valid position--force pairs in the separate `nraise=20` trajectory. The adapter
uses species-aware symmetry and separate Ga/O masses; the installed HELD CLI
itself officially targets monoatomic BCC/FCC/HCP systems.

## Result

- Complete configurations: 594
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 574.7544 fs
- Mean sampled temperature: 291.062 K
- Latest complete-frame temperature: 324.211 K
- Cutoff: 5.5 Angstrom
- Coefficients: 1,008
- Equations per frame: 240
- Per-frame rank: 240 (underdetermined relative to 1,008 coefficients)
- Global equations: 142,560
- Global design rank: 1,008 (full rank)
- Global condition number: 133.82
- Aggregation: global least squares
- Force-component RMSE: 0.2089 eV/Angstrom
- Force-vector RMSE: 0.3618 eV/Angstrom
- Frequency range: -1.206 to 24.873 THz

This is not a validated HELD prediction or a converged finite-temperature
spectrum. The full 574.8 fs sample includes the strong initial cooling
transient but now averages 291 K, while the newest complete frame remains above 300 K.
Imaginary branches must be treated as a
trajectory diagnostic rather than proof of physical instability.

Important outputs are `Ga2O3-BL-001-4A-HELD-experimental.png`,
`held_dispersion.dat`, `held_coefficients.csv`,
`held_step_coefficients.csv`, and `held_summary.json`.

## Frame-by-frame dashboard

`Ga2O3-BL-001-4A-HELD-step-dashboard.gif` contains all 594 frames at 8 fps.
At the 5.5 Angstrom cutoff each individual frame
has only 240 equations for 1,008 coefficients, so its per-frame HELD fit is
underdetermined. The near-zero per-frame force errors and the animated
per-frame spectra must not be interpreted as validated predictions. The
dashboard remains useful for temperature, structure, and qualitative
trajectory evolution.

Each frame contains its HELD dispersion, temperature history,
per-frame force-fit RMSE, and current atomic structure. The final dashboard is
also available as `Ga2O3-BL-001-4A-HELD-step-dashboard-last-frame.png`.

Across the updated dashboard the per-frame frequency range is -10.93 to 24.94 THz.
Cached per-frame dispersions and
machine-readable diagnostics are stored in `held_step_dashboard_cache.npz`
and `held_step_dashboard_summary.json`.
