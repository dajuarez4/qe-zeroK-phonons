# Experimental HELD diagnostic for the 4 Angstrom bilayer

This directory contains an experimental multi-species Ga/O HELD fit using the
18 valid position--force pairs in the separate `nraise=20` trajectory. The adapter
uses species-aware symmetry and separate Ga/O masses; the installed HELD CLI
itself officially targets monoatomic BCC/FCC/HCP systems.

## Result

- Complete configurations: 18
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 17.4168 fs
- Mean sampled temperature: 221.381 K
- Coefficients: 144
- Equations per frame: 240
- Per-frame rank: 144 (full rank)
- Global equations: 4,320
- Global design rank: 144
- Global condition number: 25.81
- Aggregation: mean of the 18 per-frame fits
- Force-component RMSE: 0.0914 eV/Angstrom
- Force-vector RMSE: 0.1583 eV/Angstrom
- Frequency range: -8.095 to 26.249 THz

This is not a validated HELD prediction or a converged finite-temperature
spectrum. The sample is only 17.4 fs, is not equilibrated, and cools
strongly during the available frames. Imaginary branches must be treated as a
short-trajectory diagnostic rather than proof of physical instability.

Important outputs are `Ga2O3-BL-001-4A-HELD-experimental.png`,
`held_dispersion.dat`, `held_coefficients.csv`,
`held_step_coefficients.csv`, and `held_summary.json`.
