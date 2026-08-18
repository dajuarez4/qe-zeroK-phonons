# TDEP phonon diagnostic for the 4 Angstrom bilayer

This directory contains a completed second-order TDEP fit using every valid
position--force pair in the interrupted Quantum ESPRESSO trajectory.

## Result

- Complete configurations: 413
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 399.6188 fs
- Mean sampled temperature: 282.418 K
- Latest complete-frame temperature: 341.657 K
- QE termination: incomplete (`JOB DONE.` is absent)
- Reference mapping: exact 4 x 2 x 1 replication to numerical precision
- Second-order cutoff: 5.5 Angstrom
- Force equations: 99,120
- First- plus second-order parameters: 1,125
- Overdetermination ratio: 88.1
- Harmonic force-fit residual RMS: 0.330 eV/Angstrom
- Force residual R-squared: 0.902
- TDEP anharmonicity measure: 0.313
- Frequency range: -1.439 to 25.897 THz

This is a pipeline diagnostic, not a converged finite-temperature phonon
spectrum. This separate `nraise=20` run has 413 adjacent complete MD frames
covering 399.6 fs. The full trajectory retains the initial cooling transient
but now averages 282 K, while the newest complete frames remain above 300 K.
The negative branches therefore cannot be
used alone to establish a physical instability. A production result requires
an equilibrated trajectory with hundreds or preferably thousands of
decorrelated configurations.

The fit omits Born effective charges and the dielectric tensor, so non-analytic
LO--TO corrections are absent for polar Ga2O3.

Important outputs are `Ga2O3-BL-001-4A-TDEP-diagnostic.png`,
`outfile.dispersion_relations`, `outfile.forceconstant`, `tdep_summary.json`,
and `extract_forceconstants.log`.

## Cumulative anharmonicity

`anharmonicity_convergence/Ga2O3-BL-001-4A-TDEP-anharmonicity-convergence.png`
is the previously generated 266-frame cumulative snapshot at the 5.5 Angstrom
cutoff. The normalized residual
`sigma_A = std(F_DFT - F_harmonic) / std(F_DFT)` rises from 0.006 at 10
frames to 0.254 at 266 frames. The residual-force standard deviation reaches
0.149 eV/Angstrom. This snapshot predates the current 413-frame fit, so it does
not establish convergence of the present trajectory. Numerical values are in
`anharmonicity_convergence/tdep_anharmonicity_convergence.csv` and
`anharmonicity_convergence/tdep_anharmonicity_summary.json`.

The actual relaxed reference printed by QE was reconstructed from the initial
`Crystallographic axes` block. This was necessary because the removed MD input
stored in Git contains different, pre-relaxation coordinates.
