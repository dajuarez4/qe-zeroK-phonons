# Ga2O3 (001): 4 Angstrom bilayer relaxation

This folder runs only the 10-atom fixed-cell relaxation. It deliberately does
not start MD from a perfect replicated structure.

The initial surface-to-surface separation between the two five-atom layers is
exactly 4.000 Angstrom. The lower layer was shifted by -0.5 Angstrom and the
upper layer by +0.5 Angstrom relative to the original 3 Angstrom structure, so
the complete slab remains centered in the 30 Angstrom cell.

The relaxation allows 150 ionic steps and uses a smaller BFGS trust radius.
The batch job explicitly rejects QE's `maximum number of steps has been
reached` termination; that message is not convergence even though QE also
prints `End of BFGS Geometry Optimization` afterward.

Submit the relaxation:

```bash
sbatch run_relax.sbatch
```

When it finishes, place this file back in this folder:

```text
Ga2O3-BL-001.relax.out
```

The next preparation stage will use the relaxed structure for a harmonic
q-grid calculation. A 4 x 2 x 1 q-grid corresponds to an 80-atom supercell;
a 4 x 4 x 1 grid corresponds to 160 atoms. Once stable force constants are
available, TDEP `canonical_configuration` can generate a 300 K displaced
structure. The final QE MD input will also contain explicit Ga/O
`ATOMIC_VELOCITIES { a.u }`, with center-of-mass drift removed and kinetic
temperature rescaled to exactly 300 K, following the Fe/IronCoreMD workflow.
