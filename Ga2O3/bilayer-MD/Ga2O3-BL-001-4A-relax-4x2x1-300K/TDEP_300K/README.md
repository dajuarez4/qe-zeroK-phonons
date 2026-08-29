# TDEP phonon diagnostic for the 4 Angstrom bilayer

This directory contains a completed second-order TDEP fit using every valid
position--force pair in the interrupted Quantum ESPRESSO trajectory.

## Result

- Complete configurations: 711
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 687.9636 fs
- Mean sampled temperature: 295.909 K
- Latest complete-frame temperature: 311.717 K
- QE termination: incomplete (`JOB DONE.` is absent)
- Reference mapping: exact 4 x 2 x 1 replication to numerical precision
- Second-order cutoff: 5.5 Angstrom
- Force equations: 170,640
- First- plus second-order parameters: 1,125
- Overdetermination ratio: 151.7
- Harmonic force-fit residual RMS: 0.372 eV/Angstrom
- Force residual R-squared: 0.879
- TDEP anharmonicity measure: 0.348
- Frequency range: -0.805 to 24.583 THz

This is a pipeline diagnostic, not a converged finite-temperature phonon
spectrum. This separate `nraise=20` run has 711 adjacent complete MD frames
covering 688.0 fs. The full trajectory retains the initial cooling transient
but now averages 296 K, while the newest complete frame remains above 300 K.
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
contains the updated 618-frame cumulative series at the 5.5 Angstrom cutoff.
The normalized residual
`sigma_A = std(F_DFT - F_harmonic) / std(F_DFT)` rises from 0.006 at 10
frames to 0.343 at 618 frames. The residual-force standard deviation reaches
0.211 eV/Angstrom. Relative to the 594-frame fit, the imaginary minimum becomes
slightly shallower, from -0.748 to -0.710 THz, while the force residual and
anharmonicity remain essentially unchanged. Numerical values are in
`anharmonicity_convergence/tdep_anharmonicity_convergence.csv` and
`anharmonicity_convergence/tdep_anharmonicity_summary.json`.

The actual relaxed reference printed by QE was reconstructed from the initial
`Crystallographic axes` block. This was necessary because the removed MD input
stored in Git contains different, pre-relaxation coordinates.
