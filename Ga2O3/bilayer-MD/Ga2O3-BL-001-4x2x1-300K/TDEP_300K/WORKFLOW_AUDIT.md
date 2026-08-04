# TDEP workflow audit: 136-frame Ga2O3 bilayer trajectory

## What is correct

- The 10-atom QE cell is copied exactly into `infile.ucposcar`.
- The 80-atom QE cell is copied exactly into `infile.ssposcar`.
- The 80-atom cell is an exact `4 x 2 x 1` replication of the 10-atom
  cell: the cell-volume ratio is 8 and the mapped atomic-position error
  is below `3.3e-15` Angstrom.
- TDEP reads 136 complete configurations, 32,640 force equations, and
  the intended 2.5 Angstrom second-order cutoff. The reported maximum
  unaliased cutoff is 5.671 Angstrom.
- The incomplete final QE frame is excluded.

## What is not physically satisfactory

The structure supplied as the TDEP ideal reference is
`../Ga2O3-BL-001-300K/Ga2O3-BL-001.md.in`. It is an MD starting input,
not a relaxed equilibrium structure. No relaxed bilayer output is
present in this workflow.

- Mean trajectory-site shift from the supplied reference: 0.376
  Angstrom RMS, 0.932 Angstrom maximum.
- Last-frame RMS displacement from the supplied reference: 0.630
  Angstrom.
- TDEP first-order predicted reference correction: 3.904 Angstrom RMS,
  8.179 Angstrom maximum.
- Mean sampled temperature: 532.732 K; only 4 of 136 complete frames
  lie between 270 and 330 K.

The very large predicted reference correction means that first-order
terms are trying to compensate for an unsuitable reference and a
non-equilibrated trajectory. `run_tdep.sh` then copies the fitted force
constants but continues to use the original `infile.ucposcar`; it does
not validate or adopt `outfile.new_ucposcar` as a new equilibrium
structure.

## Current diagnostic fit

- TDEP first-order reference-force RMSE: 0.232 eV/Angstrom.
- TDEP harmonic force-fit residual RMS: 0.909 eV/Angstrom.
- TDEP frequency range: -1.467 to 12.925 THz.
- TDEP normalized force-residual anharmonicity: 0.800358.
- TDEP residual-force standard deviation: 0.524768 eV/Angstrom.
- TDEP residual R-squared: 0.359427.
- HELD mean-coefficient force-component RMSE: 3.430 eV/Angstrom.
- HELD frequency range: -13.984 to 23.381 THz.

These are diagnostic spectra, not validated 300 K phonons.

## Required corrected workflow

1. Relax the 10-atom bilayer geometry while preserving the intended 2D
   boundary conditions and vacuum.
2. Build the 80-atom `4 x 2 x 1` cell from that relaxed structure.
3. Equilibrate the relaxed supercell near 300 K and discard the warm-up.
4. Collect a multi-picosecond production trajectory with a stationary
   temperature distribution near 300 K.
5. Use the relaxed 10-atom and 80-atom structures as TDEP's unit and
   simulation-supercell references.
6. Refit and verify that the first-order reference correction is small
   and that the dispersion is converged with trajectory length and
   cutoff.
