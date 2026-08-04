# Fe2O3 publication analysis from completed calculations

This directory contains analyses that require no additional DFT calculations.

## Inputs

- `../eos_birch_murnaghan/dftu_eos_diagnostics.csv`
- `../eos_birch_murnaghan/all_completed_birch_murnaghan_summary.txt`
- `../phonopy_V100_2x2x2/projected_dos.dat`
- `../phonopy_V100_2x2x2/band.hdf5`

The V100 projected DOS was evaluated on the existing 20x20x20 force-constant
mesh. Oxygen atoms are indices 1--6 and iron atoms are indices 7--10.

## Reproduce

Generate the atom-projected DOS with the Python environment containing
Phonopy:

```bash
cd ../phonopy_V100_2x2x2
phonopy-load phonopy_disp.yaml --mesh 20 20 20 --gamma-center \
  --pdos '1 2 3 4 5 6, 7 8 9 10' --nowritemesh
```

Then run:

```bash
/opt/anaconda3/envs/siesta_env/bin/python analyze_existing_results.py
```

## Interpretation notes

- Γ-mode Fe/O weights are squared normalized mass-weighted eigenvector
  components.
- Tiny negative acoustic frequencies are numerical translational residuals.
- Raman/IR activities and LO--TO splitting require dielectric-response data
  that are not present in the completed calculations.
- Linear pressure slopes in the summary are descriptive across the sampled
  pressure interval, not substitutes for a physical pressure-dependent model.
