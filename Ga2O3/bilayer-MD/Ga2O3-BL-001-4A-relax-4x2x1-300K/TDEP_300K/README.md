# TDEP phonon diagnostic for the 4 Angstrom bilayer

This directory contains a completed second-order TDEP fit using every valid
position--force pair in the interrupted Quantum ESPRESSO trajectory.

## Result

- Complete configurations: 203
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 196.4228 fs
- Mean sampled temperature: 249.236 K
- Latest complete-frame temperature: 236.040 K
- QE termination: incomplete (`JOB DONE.` is absent)
- Reference mapping: exact 4 x 2 x 1 replication to numerical precision
- Second-order cutoff: 2.5 Angstrom
- Force equations: 48,720
- First- plus second-order parameters: 171
- Overdetermination ratio: 284.9
- Harmonic force-fit residual RMS: 0.662 eV/Angstrom
- Force residual R-squared: 0.551
- TDEP anharmonicity measure: 0.670
- Frequency range: -3.039 to 19.498 THz

This is a pipeline diagnostic, not a converged finite-temperature phonon
spectrum. This separate `nraise=20` run has 203 adjacent complete MD frames
covering 196.4 fs. The full trajectory retains the initial cooling transient
and averages 249 K despite later intervals closer to 300 K. The negative
branches therefore cannot be
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
shows cumulative fits at 10-frame intervals through the final 203 frames. The
normalized residual
`sigma_A = std(F_DFT - F_harmonic) / std(F_DFT)` rises from 0.301 at 10
frames to 0.670 at 203 frames. The residual-force standard deviation likewise
rises from 0.103 to 0.382 eV/Angstrom. This is not convergence to a stable
anharmonicity value: later frames broaden the sampled force distribution while
the cumulative sample still mixes the cooling transient with recent 300 K
frames. Numerical values are in
`anharmonicity_convergence/tdep_anharmonicity_convergence.csv` and
`anharmonicity_convergence/tdep_anharmonicity_summary.json`.

The actual relaxed reference printed by QE was reconstructed from the initial
`Crystallographic axes` block. This was necessary because the removed MD input
stored in Git contains different, pre-relaxation coordinates.
