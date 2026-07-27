# Fe2O3 V100 finite-displacement phonons

This is a separate Phonopy + Quantum ESPRESSO workflow. It does not modify
the existing DFPT calculation, EOS outputs, or their save directories.

## Model

- source structure: final coordinates from
  `../eos_birch_murnaghan/alpha_Fe2O3v100.eos.out`
- cell: V100 10-atom rhombohedral cell, 100.625 A^3
- magnetic order: collinear AFM, Fe1 up and Fe2 down
- electronic model: PBE+U, U(Fe 3d) = 4 eV, `ortho-atomic` projectors
- supercell: 2x2x2, 80 atoms
- finite displacement: central plus/minus pairs, 0.02 bohr (0.01058 angstrom)
- independent force calculations: 6
- supercell k mesh: 4x4x4

V100 is the structure matching the attempted phonon calculation, but it is
not the zero-pressure EOS minimum. Its relaxation output reports about
99 kbar; the fitted minimum is closer to V106. Use this V100 result as the
compressed-volume point and prepare V106 separately for the near-equilibrium
phonons.

Phonopy identified type-3 AFM magnetic space group UNI 1336. The magnetic
moments used during displacement generation were `0` on O, `+4` on Fe1, and
`-4` on Fe2. The exact displacement metadata is stored in
`phonopy_disp.yaml`.

This avoids the QE 7.0 `ph.x` limitation for DFPT+U with
`U_projection_type='ortho-atomic'`. Each displacement is instead evaluated
with `pw.x`, which supports the projector already used in the EOS.

## Included files

- `unitcell.in`: relaxed V100 cell supplied to Phonopy
- `phonopy_disp.yaml`: complete Phonopy displacement and magnetic metadata
- `supercell.in`: undisplaced 80-atom supercell
- `supercell-001.in` through `supercell-006.in`: raw displaced structures
- `displacements/disp-###/alpha_Fe2O3.fd.scf.in`: ready-to-run QE inputs
- `run_displacements.sbatch`: six-task SLURM array using the existing Jakar
  QE launch configuration
- `collect_forces.sh`: validates all jobs, creates `FORCE_SETS`, and builds
  force constants
- `band.conf`: same rhombohedral path used by the earlier Fe2O3 workflow
- `mesh.conf`: 20x20x20 phonon DOS and 0--1000 K harmonic thermodynamics

## Run on Jakar

Copy the complete `Fe2O3` directory to scratch so the relative pseudopotential
path remains valid. Enter this directory and submit:

```bash
sbatch run_displacements.sbatch
```

The array runs two 80-atom calculations concurrently. To run only one at a
time, change `%2` to `%1`. Check progress with:

```bash
grep -l 'JOB DONE' displacements/disp-*/alpha_Fe2O3.fd.scf.out
grep -H 'total magnetization' displacements/disp-*/alpha_Fe2O3.fd.scf.out
```

There must be six `JOB DONE` lines. The total magnetization should remain near
zero and the AFM state must not collapse. Do not relax any displaced
structure: these are force-only SCF calculations.

The batch script temporarily disables Bash `nounset` while sourcing Intel
oneAPI. This is required because `setvars.sh` may reference optional unset
variables such as `OCL_ICD_FILENAMES`; strict checking is restored immediately
after oneAPI initialization.

The QE launch uses two k-point pools and `-ndiag 1`. The latter deliberately
avoids MKL ScaLAPACK distributed subspace diagonalization, which segfaults for
this 80-atom calculation with the Jakar QE 7.0/OpenMPI/MKL build. MKL and
OpenMP are restricted to one thread per MPI rank to prevent oversubscription.

## Collect forces

Activate the Python environment containing Phonopy, then run:

```bash
bash collect_forces.sh
```

The script refuses to proceed if any output is absent or incomplete. It passes
the six outputs to Phonopy in the same order as the six displaced supercells.

## Dispersion, DOS, and thermal properties

After `collect_forces.sh` succeeds:

```bash
phonopy-load phonopy_disp.yaml --config band.conf --save
phonopy-load phonopy_disp.yaml --config mesh.conf
```

The first command creates the band data and PDF. The second creates the DOS
and harmonic thermal-property files. Inspect the acoustic modes around Gamma
before interpreting small negative values.

This workflow does not include the non-analytical correction, so the LO-TO
splitting at Gamma is absent. Born charges and the dielectric tensor can be
added later from a compatible calculation if needed.

## Recreating the inputs

The displaced supercells were generated with Phonopy 2.31.2:

```bash
phonopy --qe -d --dim 2 2 2 --amplitude 0.02 --pm \
  --magmom 0 0 0 0 0 0 4 4 -4 -4 -c unitcell.in
python prepare_qe_inputs.py
```

Rerunning these commands overwrites generated structures, so it is unnecessary
unless the unit cell, supercell, or displacement amplitude is intentionally
changed.
