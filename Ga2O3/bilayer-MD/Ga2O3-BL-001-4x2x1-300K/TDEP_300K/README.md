# Diagnostic TDEP phonons: 80-atom Ga2O3 (001) bilayer

TDEP successfully mapped the 80-atom `4 x 2 x 1` simulation supercell
onto the 10-atom phonon unit cell and fitted second-order force
constants with a 2.5 Angstrom cutoff.

## Current result

- Complete configurations: 136
- Simulated time represented: 131.594 fs
- Mean sampled temperature: 532.732 K
- Latest complete-frame temperature: 460.184 K
- Incomplete position blocks excluded: 1
- Force equations: 32,640
- First- plus second-order parameters: 171
- Overdetermination ratio: 190.9
- First-order reference-force RMSE: approximately 0.232 eV/Angstrom
- Harmonic force-fit residual RMS: approximately 0.909 eV/Angstrom
- Harmonic force-fit residual standard deviation: 0.525 eV/Angstrom
- Frequency range: approximately -1.467 to 12.925 THz
- Requested cutoff: 2.5 Angstrom
- Maximum unaliased supercell cutoff reported by TDEP: 5.671 Angstrom

This is a pipeline diagnostic, not a converged 300 K spectrum. The
trajectory is far above its 300 K target and spans only 132 fs, so it
is neither equilibrated nor representative of a 300 K
thermal ensemble. Negative branches remain and must not yet be
interpreted as proof of physical instability.

`WORKFLOW_AUDIT.md` documents an additional critical issue: the cell
and supercell mapping is exact, but the reference coordinates are not
relaxed. TDEP predicts a 3.904 Angstrom RMS reference correction, which
is far too large for a trustworthy equilibrium force-constant fit.

Ga2O3 is polar, but this fit does not include a dielectric tensor or
Born effective charges, so non-analytic LO-TO corrections are absent.

## Important files

- `Ga2O3-BL-001-4x2x1-TDEP-diagnostic.png`: phonon plot
- `tdep_summary.json`: machine-readable dispersion diagnostics
- `WORKFLOW_AUDIT.md`: lattice mapping and reference-structure audit
- `anharmonicity_convergence/Ga2O3-BL-001-TDEP-anharmonicity-convergence.png`:
  cumulative convergence of the TDEP force-residual anharmonicity
- `anharmonicity_convergence/tdep_anharmonicity_convergence.csv`: values
  at 20, 40, 60, 80, 100, 120, and 136 frames
- `anharmonicity_convergence/tdep_anharmonicity_summary.json`:
  machine-readable anharmonicity diagnostics

## TDEP anharmonicity measure

For the 136-frame fit, TDEP reports
`sigma_A = std(F_DFT - F_harmonic) / std(F_DFT) = 0.800358`.
The residual-force standard deviation is 0.524768 eV/Angstrom and the
corresponding residual R-squared is 0.359427. The cumulative value is
roughly 0.79--0.82 after 40 frames, but this apparent numerical
stability does not make it a 300 K material property: the reference is
unrelaxed and the sampled mean temperature is 532.732 K.
- `outfile.forceconstant`: fitted second-order force constants
- `outfile.dispersion_relations`: 30 branches in THz
- `extract_forceconstants.log`: fit diagnostics
- `phonon_dispersion_relations.log`: dispersion log
