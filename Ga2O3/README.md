# Ga2O3 zero-K phonons

This directory contains a completed Quantum ESPRESSO harmonic-phonon
post-processing chain on a `2 x 2 x 2` q grid:

```text
Ga2O3.ph.in -> Ga2O3.dyn1...8 -> Ga2O3.fc -> Ga2O3.freq(.gp)
       ph.x             q2r.x                 matdyn.x
```

`Ga2O3.ph.out`, `Ga2O3.q2r.out`, and `Ga2O3.matdyn.out` all end with
`JOB DONE`.

## Reproduce the dispersion analysis

From this directory:

```bash
/usr/bin/python3 analyze_ga2o3_phonons.py
```

The script reads `Ga2O3.freq.gp` and writes:

- `Ga2O3_phonon_dispersion.png`
- `Ga2O3_Gamma_phonons_sticks.png`
- `Ga2O3_Gamma_freq.dat`
- `Ga2O3_phonon_summary.txt`

The reciprocal-coordinate path is copied from the supplied
`Ga2O3.matdyn.in`. Because the supplied cell has only inversion symmetry,
the coordinate endpoints are shown explicitly rather than assigning
hexagonal high-symmetry labels that may not apply to this structure.

## Validation result

The supplied force constants give a minimum interpolated frequency of
`-76.0629 cm-1`. At Gamma, applying `asr='crystal'` leaves an unstable
optical mode at `-39.4666 cm-1`; the other near-zero values are acoustic
modes. Negative frequencies are retained in the plots.

This does **not** mean that beta-Ga2O3 is physically unstable. It means
that the supplied force constants do not reproduce the expected stable
beta-Ga2O3 spectrum. The DFPT output detects only the `C_i` point group,
whereas an appropriately symmetry-preserving beta-Ga2O3 structure should
retain monoclinic `C_2h` symmetry. The supplied `matdyn` path was also
copied from the Al2O3 example instead of using the monoclinic beta-Ga2O3
path (`Gamma-Y-F-L-I/I1-Z-F1 | Y-X1 | X-Gamma-N | M-Gamma`).

The structure/cell and electronic ground state should therefore be
corrected or verified before recomputing DFPT. A tightly converged
symmetry-preserving relaxation and converged k and q meshes are required.

There are also interpolation-quality warnings worth converging: `q2r.x`
reports an imaginary-term sum of `4.37e-6`, and `matdyn.x` reports a
maximum dynamical-matrix non-Hermiticity of about `0.003964` (roughly
`16.6%` relative for the last sampled point). These do not change the
`JOB DONE` status, but reinforce that the imaginary branches should be
treated as calculation artifacts until the intended beta-phase symmetry
and convergence have been restored.

## Phonon DOS

The DOS input mirrors the Al2O3 example:

```bash
sbatch run_phdos.sbatch
```

The job runs `matdyn.x` using `Ga2O3.phdos.in` and then writes
`Ga2O3_phdos.png`. The local macOS environment used for post-processing
does not contain `matdyn.x`, so the DOS calculation must be run in the
cluster Quantum ESPRESSO environment.

## Reproducibility limitation

The added directory does not contain the Ga2O3 SCF input, pseudopotentials,
or `tmp/Ga2O3-Bands.save` data. The existing dynamical matrices and force
constants are sufficient for q2r/matdyn post-processing, but a new DFPT
calculation from the electronic ground state requires those missing SCF
artifacts.
