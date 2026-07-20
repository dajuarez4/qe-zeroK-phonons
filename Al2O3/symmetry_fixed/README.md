# Symmetry-fixed alpha-Al2O3 phonons

This is a separate improvement workflow. Nothing in the parent `Al2O3/`
calculation has been removed or overwritten.

## Inspection result

The existing Al2O3 spectrum is physically reasonable:

- the structure is alpha-Al2O3 with space group `R-3c` (No. 167)
- every optical branch is positive
- the minimum value, about `-9.64 cm-1`, is confined to acoustic modes
  near Gamma and is a normal numerical artifact at this scale
- `pw.x`, `ph.x`, `q2r.x`, and `matdyn.x` completed successfully

The main numerical issue is symmetry handling. The existing SCF uses
`ibrav=0`; QE explicitly warns against that representation and detects only
identity plus inversion (2 operations), although the structure has 12
`R-3c` operations. The resulting `matdyn.x` calculations repeatedly report
a maximum dynamical-matrix non-Hermiticity of `0.017658`.

This folder retains the same geometry, PBE functional, PAW
pseudopotentials, `80/640 Ry` cutoffs, and `12x12x12` k grid, but uses
QE `ibrav=5` with the exact rhombohedral metric. It also uses:

- symmetry-averaged coordinates (changes below `1e-6` fractionally)
- tighter SCF convergence (`1d-10`)
- dielectric tensor and Born effective charges
- the correct disconnected rhombohedral path
- a `2x2x2` test q grid and a `4x4x4` production input

## Run order

From this directory:

```bash
sbatch run_scf.sbatch
```

Only after `alpha_Al2O3.scf.out` ends with `JOB DONE`, submit:

```bash
sbatch run_phonons.sbatch
```

The first phonon job uses the `2x2x2` test grid. If its outputs are clean,
change the phonon input named in `run_phonons.sbatch` to
`alpha_Al2O3.ph_grid_4x4x4.in` and rerun for production convergence.

The path is:

```text
Gamma-L-B1 | B-Z-Gamma-X | Q-F-P1-Z | L-P
```
