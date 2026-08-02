# Fe2O3 V106 finite-displacement phonons

This folder is one volume point in the five-volume Phonopy/QHA series
`V100, V102, V104, V106, V108`.

- EOS source: `../eos_birch_murnaghan/alpha_Fe2O3.eosv106.out`
- primitive-cell volume: 106.662445 A^3
- structure: final fixed-volume relaxed coordinates from the EOS output
- model: collinear AFM PBE+U, U(Fe 3d) = 4 eV, ortho-atomic projectors
- supercell: 2x2x2 (80 atoms)
- displacements: three central plus/minus pairs, 0.02 bohr
- independent QE force calculations: 6
- supercell k mesh: 4x4x4

Submit the force calculations from this directory:

```bash
sbatch run_displacements.sbatch
```

After all six outputs contain `JOB DONE`, activate Phonopy and run:

```bash
bash collect_forces.sh
phonopy-load phonopy_disp.yaml --config band.conf --save
phonopy-load phonopy_disp.yaml --config mesh.conf
```

Keep all electronic, magnetic, displacement, and sampling settings identical
between volumes. Do not relax a displaced supercell.
