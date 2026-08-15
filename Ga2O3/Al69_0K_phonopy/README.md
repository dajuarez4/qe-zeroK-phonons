# 0 K phonons for the 69% Al beta-(Al,Ga)2O3 alloy

This directory converts `Bands_Alloy_69pct.fdf` from SIESTA to a reproducible
Quantum ESPRESSO + Phonopy finite-displacement workflow.

## Structure and composition

- Source cell: 40 atoms, copied unchanged as
  `Bands_Alloy_69pct.source.fdf` (SHA256
  `b3d38f03810a6971c76c964f62653a607107820125d24826584a3e23f0b5786f`).
- Counts: Al11 Ga5 O24, or Al fraction 11/16 = 0.6875 on the cation sites.
- Formula: `(Al0.6875 Ga0.3125)2O3` (eight formula units in the source cell).
- DFT: non-spin-polarized PBE, fixed occupations, `ecutwfc=80 Ry`,
  `ecutrho=640 Ry`.
- Pseudopotentials: the same PSLibrary PBE PAW files already used in this
  repository for Ga2O3 and Al2O3.

The supplied coordinates were relaxed using SIESTA. Because forces and stress
are code- and pseudopotential-dependent, a QE zero-pressure `vc-relax` is the
required first stage before computing harmonic force constants.

## Workflow

1. Relax the 40-atom alloy cell with QE:

   ```bash
   sbatch run_relax.sbatch
   ```

2. Once `01_vc-relax.out` contains `JOB DONE`, extract the final geometry and
   generate central plus/minus displacements:

   ```bash
   python3 extract_relaxed_structure.py
   bash generate_displacements.sh
   ```

   The source cell is already a 2x1x2 expansion of the underlying 10-atom
   beta-Ga2O3 cell. Phonopy therefore uses `DIM = 1 2 1`, producing an
   approximately isotropic 80-atom force supercell. The displaced-supercell
   k mesh is 6x6x6. A dry run on the supplied structure finds space group P1
   and generates 240 calculations (40 inequivalent atoms x 3 directions x
   plus/minus displacements). The final count is regenerated after QE relax.

3. `generate_displacements.sh` prints the exact Slurm array command. Submit
   that command; for example, if it generated N folders:

   ```bash
   sbatch --array=0-$((N-1))%4 run_displacements.sbatch
   ```

   The `%4` limits the run to four concurrent force calculations and can be
   adjusted for the allocation. The low-symmetry ordered alloy is
   computationally expensive; the script determines the count rather than
   hard-coding it in the Slurm file.

4. After every force calculation contains `JOB DONE`:

   ```bash
   bash collect_forces.sh
   ```

   This builds `FORCE_SETS`, `FORCE_CONSTANTS`, the band structure along the
   path supplied in the FDF (Gamma-A-Z-M-L-V-Gamma), total/projected DOS,
   0--1000 K harmonic thermal properties, and PNG/PDF plots.

## Important limitations

- The ordered 40-atom model represents one atomic arrangement at 68.75% Al;
  it is not a configurational average. Several symmetry-inequivalent alloy
  arrangements or an SQS are required for disorder uncertainty.
- This is a harmonic 0 K force-constant workflow. The reported thermal
  functions use those 0 K force constants; anharmonic lifetimes and lattice
  thermal conductivity require third-order force constants or MD/TDEP.
- Non-analytical LO-TO corrections are not enabled. Born effective charges
  and the dielectric tensor from a separate QE calculation are required for
  accurate polar modes near Gamma.
- Check convergence against a larger displacement supercell and denser k mesh
  before treating small imaginary frequencies as physical instabilities.
