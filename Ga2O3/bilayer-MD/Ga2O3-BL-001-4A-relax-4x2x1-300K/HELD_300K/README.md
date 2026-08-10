# Experimental HELD diagnostic for the 4 Angstrom bilayer

This directory contains an experimental multi-species Ga/O HELD fit using the
139 valid position--force pairs in the separate `nraise=20` trajectory. The adapter
uses species-aware symmetry and separate Ga/O masses; the installed HELD CLI
itself officially targets monoatomic BCC/FCC/HCP systems.

## Result

- Complete configurations: 139
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 134.4964 fs
- Mean sampled temperature: 232.888 K
- Coefficients: 144
- Equations per frame: 240
- Per-frame rank: 144 (full rank)
- Global equations: 33,360
- Global design rank: 144
- Global condition number: 16.04
- Aggregation: mean of the 139 per-frame fits
- Force-component RMSE: 0.2976 eV/Angstrom
- Force-vector RMSE: 0.5155 eV/Angstrom
- Frequency range: -2.226 to 23.861 THz

This is not a validated HELD prediction or a converged finite-temperature
spectrum. The full 134.5 fs sample includes the strong initial cooling
transient and averages 233 K despite later intervals closer to 300 K.
Imaginary branches must be treated as a
short-trajectory diagnostic rather than proof of physical instability.

Important outputs are `Ga2O3-BL-001-4A-HELD-experimental.png`,
`held_dispersion.dat`, `held_coefficients.csv`,
`held_step_coefficients.csv`, and `held_summary.json`.

## Frame-by-frame dashboard

`Ga2O3-BL-001-4A-HELD-step-dashboard.gif` animates all 139 complete frames
at 8 fps. Each frame contains its HELD dispersion, temperature history,
per-frame force-fit RMSE, and current atomic structure. The final dashboard is
also available as `Ga2O3-BL-001-4A-HELD-step-dashboard-last-frame.png`.

Across individual frames the frequency range is -12.77 to 28.48 THz, showing
that the per-step fits vary much more strongly than the mean HELD spectrum.
Its mean per-step in-sample force RMSE is 0.144 eV/Angstrom; this must not be
confused with the 0.298 eV/Angstrom residual obtained by applying the mean
coefficients to the complete trajectory. Cached per-frame dispersions and
machine-readable diagnostics are stored in `held_step_dashboard_cache.npz`
and `held_step_dashboard_summary.json`.
