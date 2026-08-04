# Experimental HELD diagnostic for the 4 Angstrom bilayer

This directory contains an experimental multi-species Ga/O HELD fit using the
nine valid position--force pairs in the interrupted trajectory. The adapter
uses species-aware symmetry and separate Ga/O masses; the installed HELD CLI
itself officially targets monoatomic BCC/FCC/HCP systems.

## Result

- Complete configurations: 9
- Represented trajectory time: 8.7084 fs
- Mean sampled temperature: 265.872 K
- Coefficients: 144
- Equations per frame: 240
- Per-frame rank: 144 (full rank)
- Global equations: 2,160
- Global design rank: 144
- Global condition number: 21.32
- Aggregation: mean of the nine per-frame fits
- Force-component RMSE: 0.0453 eV/Angstrom
- Force-vector RMSE: 0.0784 eV/Angstrom
- Frequency range: -5.779 to 25.473 THz

This is not a validated HELD prediction or a converged finite-temperature
spectrum. The sample is shorter than 9 fs, is not equilibrated, and cools
strongly during the available frames. Imaginary branches must be treated as a
short-trajectory diagnostic rather than proof of physical instability.

Important outputs are `Ga2O3-BL-001-4A-HELD-experimental.png`,
`held_dispersion.dat`, `held_coefficients.csv`,
`held_step_coefficients.csv`, and `held_summary.json`.
