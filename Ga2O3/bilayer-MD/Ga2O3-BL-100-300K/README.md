# Ga2O3 (100) bilayer: 300 K QE molecular dynamics

This folder translates `../Ga2O3-GGA-BL-100.fdf` to a fixed-cell
Quantum ESPRESSO Born-Oppenheimer MD run.

- Ensemble: NVT using stochastic velocity rescaling (`svr`)
- Target temperature: 300 K
- Time step: 20 QE time units = approximately 0.968 fs
- Length: 10,000 steps = approximately 9.68 ps
- Cell and positions: copied exactly from the source FDF
- Electronic setup: PBE, 80/640 Ry cutoffs, `6 x 4 x 1` test k-point mesh
- Slab electrostatics: QE two-dimensional isolation correction

Submit from this directory:

```bash
sbatch run_md.sbatch
```

The pseudopotentials are read from `../../pseudo`. The input assumes the
supplied bilayer geometry is already suitable for MD; relax it first if
the FDF coordinates are only an unrelaxed starting structure.

The source FDF uses a much denser `24 x 14 x 1` mesh. The present
`6 x 4 x 1` mesh is intended for short test runs. Converge the k-point
mesh before using energies, forces, or trajectories for production
analysis.

For a continuation after a cleanly stopped run, preserve `tmp/` and the
output, change `restart_mode` to `'restart'`, and submit again. Do not
run two jobs against the same `tmp/` directory.
