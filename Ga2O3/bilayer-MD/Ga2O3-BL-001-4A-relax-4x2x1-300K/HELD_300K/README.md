# Experimental HELD diagnostic for the 4 Angstrom bilayer

This directory contains an experimental multi-species Ga/O HELD fit using the
30 valid position--force pairs in the interrupted trajectory. The adapter
uses species-aware symmetry and separate Ga/O masses; the installed HELD CLI
itself officially targets monoatomic BCC/FCC/HCP systems.

## Result

- Complete configurations: 30
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 29.028 fs
- Mean sampled temperature: 191.147 K
- Coefficients: 144
- Equations per frame: 240
- Per-frame rank: 144 (full rank)
- Global equations: 7,200
- Global design rank: 144
- Global condition number: 14.54
- Aggregation: mean of the 30 per-frame fits
- Force-component RMSE: 0.1255 eV/Angstrom
- Force-vector RMSE: 0.2174 eV/Angstrom
- Frequency range: -4.197 to 25.256 THz

This is not a validated HELD prediction or a converged finite-temperature
spectrum. The sample is only 29.0 fs, is not equilibrated, and cools
strongly during the available frames. Imaginary branches must be treated as a
short-trajectory diagnostic rather than proof of physical instability.

Important outputs are `Ga2O3-BL-001-4A-HELD-experimental.png`,
`held_dispersion.dat`, `held_coefficients.csv`,
`held_step_coefficients.csv`, and `held_summary.json`.
