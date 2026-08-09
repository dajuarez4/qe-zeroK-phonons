# Experimental HELD diagnostic for the 4 Angstrom bilayer

This directory contains an experimental multi-species Ga/O HELD fit using the
81 valid position--force pairs in the separate `nraise=20` trajectory. The adapter
uses species-aware symmetry and separate Ga/O masses; the installed HELD CLI
itself officially targets monoatomic BCC/FCC/HCP systems.

## Result

- Complete configurations: 81
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 78.3756 fs
- Mean sampled temperature: 207.300 K
- Coefficients: 144
- Equations per frame: 240
- Per-frame rank: 144 (full rank)
- Global equations: 19,440
- Global design rank: 144
- Global condition number: 17.62
- Aggregation: mean of the 81 per-frame fits
- Force-component RMSE: 0.2485 eV/Angstrom
- Force-vector RMSE: 0.4304 eV/Angstrom
- Frequency range: -4.354 to 24.028 THz

This is not a validated HELD prediction or a converged finite-temperature
spectrum. Recent frames have reached the 300 K target, but the full 78.4 fs
sample includes the strong initial cooling transient and averages 207 K.
Imaginary branches must be treated as a
short-trajectory diagnostic rather than proof of physical instability.

Important outputs are `Ga2O3-BL-001-4A-HELD-experimental.png`,
`held_dispersion.dat`, `held_coefficients.csv`,
`held_step_coefficients.csv`, and `held_summary.json`.

## Frame-by-frame dashboard

`Ga2O3-BL-001-4A-HELD-step-dashboard.gif` animates all 81 complete MD
frames at 8 fps. Each frame contains its HELD dispersion, temperature history,
per-frame force-fit RMSE, and current atomic structure. The final dashboard is
also available as `Ga2O3-BL-001-4A-HELD-step-dashboard-last-frame.png`.

Across individual frames the frequency range is -12.50 to 28.48 THz, showing
that the per-step fits vary much more strongly than the mean HELD spectrum.
The mean per-step in-sample force RMSE is 0.128 eV/Angstrom; this must not be
confused with the 0.248 eV/Angstrom residual obtained by applying the mean
coefficients to the complete trajectory. Cached per-frame dispersions and
machine-readable diagnostics are stored in `held_step_dashboard_cache.npz`
and `held_step_dashboard_summary.json`.
