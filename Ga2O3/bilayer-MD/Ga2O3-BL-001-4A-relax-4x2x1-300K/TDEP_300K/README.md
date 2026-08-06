# TDEP phonon diagnostic for the 4 Angstrom bilayer

This directory contains a completed second-order TDEP fit using every valid
position--force pair in the interrupted Quantum ESPRESSO trajectory.

## Result

- Complete configurations: 18
- Unfinished coordinate blocks excluded: 1
- Represented trajectory time: 17.4168 fs
- Mean sampled temperature: 221.381 K
- Latest complete-frame temperature: 174.833 K
- QE termination: incomplete (`JOB DONE.` is absent)
- Reference mapping: exact 4 x 2 x 1 replication to numerical precision
- Second-order cutoff: 2.5 Angstrom
- Force equations: 4,320
- First- plus second-order parameters: 171
- Overdetermination ratio: 25.3
- Harmonic force-fit residual RMS: 0.270 eV/Angstrom
- Force residual R-squared: 0.870
- TDEP anharmonicity measure: 0.361
- Frequency range: -3.493 to 26.820 THz

This is a pipeline diagnostic, not a converged finite-temperature phonon
spectrum. This separate `nraise=20` run has 18 adjacent complete MD frames
covering only 17.4 fs, and its temperature falls from about 302 K to roughly
175 K. The negative branches therefore cannot be
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
