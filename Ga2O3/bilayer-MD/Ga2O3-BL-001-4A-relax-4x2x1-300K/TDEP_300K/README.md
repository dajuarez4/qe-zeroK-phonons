# TDEP phonon diagnostic for the 4 Angstrom bilayer

This directory contains a completed second-order TDEP fit using every valid
position--force pair in the interrupted Quantum ESPRESSO trajectory.

## Result

- Complete configurations: 48
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 46.4448 fs
- Mean sampled temperature: 185.364 K
- Latest complete-frame temperature: 154.623 K
- QE termination: incomplete (`JOB DONE.` is absent)
- Reference mapping: exact 4 x 2 x 1 replication to numerical precision
- Second-order cutoff: 2.5 Angstrom
- Force equations: 11,520
- First- plus second-order parameters: 171
- Overdetermination ratio: 67.4
- Harmonic force-fit residual RMS: 0.340 eV/Angstrom
- Force residual R-squared: 0.829
- TDEP anharmonicity measure: 0.413
- Frequency range: -3.023 to 24.771 THz

This is a pipeline diagnostic, not a converged finite-temperature phonon
spectrum. The 48 adjacent complete MD frames cover only 46.4 fs and the
temperature falls from about 302 K to roughly 155 K. The negative branches therefore cannot be
used alone to establish a physical instability. A production result requires
an equilibrated trajectory with hundreds or preferably thousands of
decorrelated configurations.

The fit omits Born effective charges and the dielectric tensor, so non-analytic
LO--TO corrections are absent for polar Ga2O3.

Important outputs are `Ga2O3-BL-001-4A-TDEP-diagnostic.png`,
`outfile.dispersion_relations`, `outfile.forceconstant`, `tdep_summary.json`,
and `extract_forceconstants.log`.

The actual relaxed reference printed by QE was reconstructed from the initial
`Crystallographic axes` block. This was necessary because the removed MD input
stored in Git contains different, pre-relaxation coordinates.
