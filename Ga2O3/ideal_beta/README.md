# Ideal beta-Ga2O3 calculation

This directory is a clean replacement for the distorted Ga2O3 calculation
in the parent directory. It uses the experimental monoclinic `C2/m`
beta-Ga2O3 structure with a 10-atom primitive cell.

## Structure source

The conventional-cell parameters and refined independent atomic sites are
from Table 3.1 of Playford's beta-Ga2O3 refinement:

<https://wrap.warwick.ac.uk/id/eprint/55164/1/WRAP_THESIS_Playford_2012.pdf>

They are consistent with the classic structure determined by Geller:

<https://doi.org/10.1063/1.1731237>

The input builder expands the five independent sites with `C2/m` symmetry
and converts the 20-atom conventional cell to the standard 10-atom
primitive cell. Do not manually perturb or round its fractional positions.

## Pseudopotentials

These exact official Quantum ESPRESSO PSLibrary files are included in
`Ga2O3/pseudo/`:

```text
Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
O.pbe-n-kjpaw_psl.1.0.0.UPF
```

They are the same PAW/PBE pseudopotentials used by the earlier calculation:

<https://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library/ga>

<https://pseudopotentials.quantum-espresso.org/legacy_tables/ps-library/o>

## Run order

First perform the SCF calculation on the ideal experimental structure:

```bash
sbatch run_ideal_scf.sbatch
```

Check that `beta_Ga2O3.scf_ideal.out` ends with `JOB DONE`, reports the
expected monoclinic symmetry, and has no electronic-convergence errors.

Then submit the initial `2x2x2` phonon calculation:

```bash
sbatch run_ideal_phonons.sbatch
```

This runs `ph.x`, `q2r.x`, the five disconnected monoclinic path segments,
the phonon DOS, and the plotting script. The path is:

```text
Gamma-Y-F-L-I | I1-Z-F1 | Y-X1 | X-Gamma-N | M-Gamma
```

It is the same path topology as the supplied stable reference plot. The
special-point coordinates were calculated from this exact primitive cell;
they were not copied from Al2O3.

After confirming that the ideal-cell test has no meaningful imaginary
optical modes, replace the phonon input in the Slurm script with
`beta_Ga2O3.ph_grid_4x4x4.in` for the converged production calculation.

`beta_Ga2O3.relax.in` is also provided as an optional fixed-cell ionic
relaxation. Do not use its result for phonons without first copying its
final fractional coordinates into a new SCF input and rerunning `pw.x`.
