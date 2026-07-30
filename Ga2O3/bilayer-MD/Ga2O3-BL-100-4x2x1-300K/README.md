# Ga2O3 (100) bilayer 4x2x1 supercell at 300 K

This is a `4 x 2 x 1` replication of the 10-atom `(100)` bilayer:

```text
10 atoms x 4 x 2 x 1 = 80 atoms
```

The in-plane lattice vectors are replicated four times along `a` and
twice along `b`; the 30 Angstrom vacuum direction is not replicated.

- Composition: 32 Ga and 48 O atoms
- Ensemble: fixed-cell NVT with stochastic velocity rescaling
- Target temperature: 300 K
- Time step: approximately 0.968 fs
- Requested length: 10,000 steps, approximately 9.68 ps
- Test k-point mesh: `2 x 2 x 1`
- Plane-wave cutoffs: 80/640 Ry

Submit from this directory:

```bash
sbatch run_md.sbatch
```

`build_supercell.py` reproducibly regenerates the QE input from the
10-atom `(100)` input in the neighboring folder.

The source geometry should be relaxed before a production TDEP
trajectory. For a quick execution test, reduce `nstep`; for a TDEP fit,
retain a long equilibrated portion containing hundreds of
configurations.
