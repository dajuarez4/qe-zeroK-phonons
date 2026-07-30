# Diagnostic TDEP phonons: 80-atom Ga2O3 (001) bilayer

TDEP successfully mapped the 80-atom `4 x 2 x 1` simulation supercell
onto the 10-atom phonon unit cell and fitted second-order force
constants with a 2.5 Angstrom cutoff.

## Current result

- Complete configurations: 30
- Simulated time represented: 29.028 fs
- Mean sampled temperature: 518.308 K
- Latest complete-frame temperature: 679.508 K
- Incomplete position blocks excluded: 1
- Force equations: 7,200
- First- plus second-order parameters: 171
- Overdetermination ratio: 42.1
- Total residual-force RMSE: approximately 0.667 eV/Angstrom
- Frequency range: approximately -4.163 to 13.871 THz
- Requested cutoff: 2.5 Angstrom
- Maximum unaliased supercell cutoff reported by TDEP: 5.671 Angstrom

This is a pipeline diagnostic, not a converged 300 K spectrum. The
trajectory heats from about 303 K to roughly 680--700 K and spans only
29 fs, so it is neither equilibrated nor representative of a 300 K
thermal ensemble. Negative branches remain and must not yet be
interpreted as proof of physical instability.

Ga2O3 is polar, but this fit does not include a dielectric tensor or
Born effective charges, so non-analytic LO-TO corrections are absent.

## Important files

- `Ga2O3-BL-001-4x2x1-TDEP-diagnostic.png`: phonon plot
- `tdep_summary.json`: machine-readable dispersion diagnostics
- `outfile.forceconstant`: fitted second-order force constants
- `outfile.dispersion_relations`: 30 branches in THz
- `extract_forceconstants.log`: fit diagnostics
- `phonon_dispersion_relations.log`: dispersion log
