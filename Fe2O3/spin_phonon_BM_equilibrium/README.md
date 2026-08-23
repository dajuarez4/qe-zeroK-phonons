# AFM-FM spin-phonon comparison at the BM equilibrium

This workflow compares collinear AFM (`++--`) and constrained FM (`++++`)
phonons at exactly the same structure.

- BM3 equilibrium volume: 106.261887 A^3 per 10-atom cell
- BM3 hexagonal parameters: a = 5.127825 A, c = 13.999151 A
- internal coordinates: periodic linear interpolation between the independently
  relaxed V104 and V106 structures
- model: PBE+U, U(Fe 3d) = 4 eV, ortho-atomic projectors
- supercells: 2x2x2 (80 atoms), common six central displacements
- supercell k mesh: 4x4x4

The BM point belongs to the existing fixed-shape EOS path. It is the energy-fit
minimum, not a new fully anisotropic zero-stress `vc-relax` structure.

## 1. Validate the two magnetic states

```bash
cd validation
sbatch run_validation.sbatch
python3 check_validation.py
```

Proceed only if both SCFs converged, AFM has total magnetization near 0, FM
has total magnetization near 20 muB per primitive cell, and the reported local
Fe moments preserve their intended signs.

## 2. Submit matched force calculations

```bash
cd ../AFM
sbatch run_displacements.sbatch
cd ../FM
sbatch run_displacements.sbatch
```

Each state has six calculations, limited to two simultaneous jobs. The FM
calculation is constrained to 160 muB in the 80-atom supercell.

## 3. Post-process after all jobs finish

```bash
cd AFM && bash collect_forces.sh
cd ../FM && bash collect_forces.sh
cd .. && python3 plot_afm_fm_comparison.py
```

The same displacement set and the AFM magnetic subgroup symmetry are used for
both fits. This deliberately avoids attributing different symmetry reduction
or different structures to spin-phonon coupling.
