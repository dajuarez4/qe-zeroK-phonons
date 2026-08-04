# Ga2O3 (001): relaxation followed by 80-atom MD

This is the same simple DFT layout as the previous 80-atom folder, with a
10-atom fixed-cell relaxation added before MD. There is no TDEP or HELD step.

The initial surface-to-surface separation between the two five-atom layers is
exactly 4.000 Angstrom. The lower layer was shifted by -0.5 Angstrom and the
upper layer by +0.5 Angstrom relative to the original 3 Angstrom structure, so
the complete slab remains centered in the 30 Angstrom cell.

`Ga2O3-BL-001-4x2x1.md.in` is included so that all MD parameters and all 80
atoms are visible before submission. Its initial coordinates are a 4 x 2 x 1
replication of the starting unit cell. During the job, `build_md_from_relax.py`
overwrites this file using the final relaxed coordinates, and only then does
QE start MD.

The relaxation allows 150 ionic steps and uses a smaller BFGS trust radius.
The batch job explicitly rejects QE's `maximum number of steps has been
reached` termination; that message is not convergence even though QE also
prints `End of BFGS Geometry Optimization` afterward.

Submit only one job:

```bash
sbatch run_relax_md.sbatch
```

The job performs, in order:

1. Relax the 10-atom Ga2O3 (001) bilayer.
2. Confirm that BFGS and QE finished successfully.
3. Build the 4 x 2 x 1 supercell from the relaxed coordinates.
4. Run the 80-atom, 300 K MD calculation.

The generated MD input and output will be:

```text
Ga2O3-BL-001-4x2x1.md.in
Ga2O3-BL-001-4x2x1.md.out
```

MD settings are retained from the original folder: 10,000 steps, `dt=20 au`,
SVR thermostat at 300 K with `nraise=200`, 80/640 Ry cutoffs, 30 Angstrom
vacuum, `assume_isolated='2D'`, and a `2 x 2 x 1` k-point mesh.
