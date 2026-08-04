# TDEP phonon diagnostic for the 4 Angstrom bilayer

This directory contains a completed second-order TDEP fit using every valid
position--force pair in the interrupted Quantum ESPRESSO trajectory.

## Result

- Complete configurations: 9
- Represented trajectory time: 8.7084 fs
- Mean sampled temperature: 265.872 K
- Latest sampled temperature: 213.662 K
- QE termination: incomplete (`JOB DONE.` is absent)
- Reference mapping: exact 4 x 2 x 1 replication to numerical precision
- Second-order cutoff: 2.5 Angstrom
- Force equations: 2,160
- First- plus second-order parameters: 171
- Overdetermination ratio: 12.6
- Harmonic force-fit residual RMS: 0.162 eV/Angstrom
- Force residual R-squared: 0.933
- TDEP anharmonicity measure: 0.259
- Frequency range: -3.723 to 25.303 THz

This is a pipeline diagnostic, not a converged finite-temperature phonon
spectrum. Nine adjacent MD frames cover less than 9 fs and the temperature
falls from about 302 K to 214 K. The negative branches therefore cannot be
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
