# Ga2O3 (001) bilayer 4x2x1 supercell at 300 K

This is a `4 x 2 x 1` replication of the 10-atom `(001)` bilayer:

```text
10 atoms x 4 x 2 x 1 = 80 atoms
```

The in-plane lattice vectors are replicated four times along `a` and
twice along `b`; the 30 Angstrom vacuum direction is not replicated.

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
10-atom `(001)` input in the neighboring folder.

The source geometry has large residual forces and should be relaxed
before a production TDEP trajectory. For a quick execution test, reduce
`nstep`; for a TDEP fit, retain a long equilibrated portion containing
hundreds of configurations.

## Finite-temperature phonon diagnostics

The currently available output has 136 complete MD frames. The
following local workflows process every complete position/force pair
and exclude the unfinished final frame:

```bash
./run_tdep.sh
/opt/anaconda3/bin/python run_held_phonons.py
```

Results are in `TDEP_300K/` and `HELD_300K/`. The sampled 131.594 fs
trajectory remains far above the target (mean 532.732 K), so these are
not converged 300 K spectra. Read the README
in each result folder before interpretation.

The per-step experimental HELD evolution is available as
`HELD_300K/Ga2O3-BL-001-4x2x1-HELD-step-dashboard.gif`.

## No-Vito displacement animation

`Ga2O3-BL-001-4x2x1-30MD-NoVito.gif` was rendered with the local
No-Vito `qe_npz_to_gif.py` utility. It contains exactly the 30 complete
MD frames and displays the step, time, temperature, and QE pressure.
The unfinished final frame is deliberately excluded.

`Ga2O3-BL-001-4x2x1-30MD-NoVito-ZOOM.gif` is the slab-focused version;
it removes most of the vacuum from the camera view without changing
the atomic coordinates.

`Ga2O3-BL-001-80atom-results-dashboard.png` summarizes the temperature
trajectory and the current TDEP and experimental HELD diagnostics.

Regenerate the No-Vito archive with:

```bash
/opt/anaconda3/bin/python qe_to_novito_npz.py \
  Ga2O3-BL-001-4x2x1.md.in Ga2O3-BL-001-4x2x1.md.out \
  -o Ga2O3-BL-001-4x2x1-30MD-NoVito.npz \
  --summary Ga2O3-BL-001-4x2x1-30MD-NoVito-summary.json \
  --md-only
```

The current trajectory spans only 29.028 fs and heats well above the
300 K target. It does not demonstrate thermal or dynamical stability.
