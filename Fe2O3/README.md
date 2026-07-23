# Alpha-Fe2O3 (hematite): PBE+U phonon workflow

This directory is a new calculation and does not modify the completed
`Al2O3/` or `Ga2O3/` calculations.  It uses the same Jakar Slurm layout as
`Al2O3/symmetry_fixed/`.

## Physical model

- phase: alpha-Fe2O3 (hematite), corundum structure, space group R-3c
- cell: 10-atom rhombohedral primitive cell (`ibrav=5`)
- experimental starting structure: Maslen et al., Acta Cryst. B 50, 435 (1994)
- magnetic state: collinear antiferromagnetic `++--` order on the four Fe sites
- total magnetization constrained to zero, as required by QE 7.0 for fixed
  occupations with spin polarization
- functional: spin-polarized PBE+U
- starting Hubbard value: `U_eff(Fe-3d) = 4.0 eV`
- pseudopotentials: PSLibrary PBE PAW
- cutoffs: 80/640 Ry
- k mesh: 8x8x8 for relaxation and SCF
- q meshes: 2x2x2 test and 4x4x4 production

The value 4.0 eV is a defensible starting point, not a universal constant.
For publication-quality results, repeat at selected U values (for example
3, 4, and 5 eV), or calculate U consistently with `hp.x`.

The inputs use the QE 7.0 Hubbard syntax installed on Jakar:

```text
lda_plus_u = .true.
Hubbard_U(Fe1) = Hubbard_U(Fe2) = 4.0 eV
U_projection_type = 'ortho-atomic'
```

Do not replace this with the newer `HUBBARD` card unless the cluster QE
installation is upgraded to 7.1 or newer.

## Files

- `alpha_Fe2O3.relax.in`: fixed-cell ionic relaxation
- `alpha_Fe2O3.scf.in`: tightly converged SCF used by `ph.x`
- `update_scf_from_relax.py`: transfers the final relaxed coordinates into
  the SCF input; `run_relax.sbatch` calls it automatically
- `alpha_Fe2O3.ph_grid_2x2x2.in`: inexpensive phonon test
- `alpha_Fe2O3.ph_grid_4x4x4.in`: production q-grid input
- `alpha_Fe2O3.q2r.in`: real-space force constants
- `alpha_Fe2O3.path[1-4].matdyn.in`: disconnected rhombohedral band path
- `alpha_Fe2O3.phdos.in`: phonon DOS
- `plot_alpha_fe2o3_phonons.py`: dispersion and DOS plots
- `hematite_experimental.cif`: provenance copy of the starting structure

## Pseudopotentials

The folder `pseudo/` must contain:

```text
Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
O.pbe-n-kjpaw_psl.1.0.0.UPF
```

Both `Fe1` and `Fe2` deliberately use the same Fe file.  The different
species labels only let QE initialize opposite spin sublattices.

MD5 checksums of the included files are:

```text
fc81f059e5c5069939230b1155715ae8  Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
e99d9cef9b487d1ca56f5b95ecd0fd7a  O.pbe-n-kjpaw_psl.1.0.0.UPF
```

## Run order

Copy the complete `Fe2O3/` directory to your Jakar scratch area and enter it.
The inputs use `pseudo_dir='./pseudo'` and `outdir='./tmp'`, so the directory
is portable and no path needs to be edited.

### 1. Relax the ions

```bash
sbatch run_relax.sbatch
```

After completion, verify:

```bash
grep 'JOB DONE' alpha_Fe2O3.relax.out
grep -A12 'Begin final coordinates' alpha_Fe2O3.relax.out | tail -13
```

The batch script runs `update_scf_from_relax.py` only after a successful QE
job.  It backs up the original SCF input as `alpha_Fe2O3.scf.in.before_relax`
and inserts the final relaxed coordinates into `alpha_Fe2O3.scf.in`.

### 2. Run the final SCF

```bash
sbatch run_scf.sbatch
```

Before phonons, check:

```bash
grep 'JOB DONE' alpha_Fe2O3.scf.out
grep -E 'total magnetization|absolute magnetization' alpha_Fe2O3.scf.out | tail
grep 'Total force' alpha_Fe2O3.scf.out | tail
```

The total magnetization should be close to zero while the absolute
magnetization remains nonzero.  Also inspect the site-resolved Fe moments in
the verbose output and confirm the intended two-up/two-down AFM state did not
collapse into another solution.

### 3. Run the 2x2x2 phonon test

```bash
sbatch run_phonons.sbatch
```

This runs `ph.x`, `q2r.x`, the four dispersion segments, the DOS, and then
attempts plotting.  If matplotlib is unavailable on the compute node, the QE
results are still valid; run the plotting script later on your local machine.

Inspect the outputs before interpreting negative frequencies:

```bash
grep 'JOB DONE' alpha_Fe2O3.ph_grid_2x2x2.out alpha_Fe2O3.q2r.out alpha_Fe2O3.path*.matdyn.out alpha_Fe2O3.phdos.out
grep -iE 'error|warning|imaginary' alpha_Fe2O3.ph_grid_2x2x2.out alpha_Fe2O3.q2r.out alpha_Fe2O3.path*.matdyn.out
```

Small acoustic values close to Gamma can be interpolation or sum-rule error.
Large negative branches extending through the Brillouin zone require checking
the relaxation, magnetic state, U value, and numerical convergence.

### 4. Production q-grid

The `4x4x4` grid is much more expensive.  Run it only after the 2x2x2 test is
clean.  In `run_phonons.sbatch`, change these two names:

```text
alpha_Fe2O3.ph_grid_2x2x2.in  -> alpha_Fe2O3.ph_grid_4x4x4.in
alpha_Fe2O3.ph_grid_2x2x2.out -> alpha_Fe2O3.ph_grid_4x4x4.out
```

Archive the 2x2x2 `dyn*`, force-constant, frequency, and DOS files first,
because the production calculation intentionally uses the same output names.

## Band path

The four disconnected rhombohedral segments are

```text
Gamma-L-B1 | B-Z-Gamma-X | Q-F-P1-Z | L-P
```

The coordinates in the `matdyn` inputs were recalculated for the experimental
hematite rhombohedral metric; they were not copied numerically from Al2O3.

## Recommended convergence checks

Before using the spectrum in a paper, test at least:

1. `8x8x8` versus `10x10x10` k points.
2. `2x2x2` versus `4x4x4` q points.
3. Fe Hubbard U sensitivity.
4. Final force and stress tolerances.
5. AFM state and local moments after every relaxation/SCF calculation.

This workflow is scalar-relativistic and collinear.  That is appropriate for
the initial harmonic phonon dispersion.  Spin-orbit coupling and noncollinear
magnetism are separate, substantially more demanding calculations needed for
the weak canting and Morin-transition anisotropy.

## Sources

- Experimental structure: [American Mineralogist Crystal Structure Database,
  hematite entry 17807](https://larixite.seescience.org/cifs/17807)
- Magnetic configurations and primitive-cell convention: [First-principles
  calculations of hematite by self-consistent DFT+U+V](https://pmc.ncbi.nlm.nih.gov/articles/PMC9941207/)
- QE 7.0 phonon guide: [Quantum ESPRESSO PHonon User's Guide](https://www.quantum-espresso.org/wp-content/uploads/2022/03/ph_user_guide.pdf)
- Pseudopotentials: [Quantum ESPRESSO pseudopotential repository](https://pseudopotentials.quantum-espresso.org/)
