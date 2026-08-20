# Ga2O3 bilayer MD animations

These animations were rendered from the 50 K, 40x40x1 LAMMPS trajectory
(`trajectory_50K.lammpstrj`, 16,000 atoms). A central region is shown so that
the motion and the two layers remain visible.

- `Ga2O3_bilayer_MD_50K_actual_motion.gif`: zoomed-out view with physical
  displacements (x1). This is the recommended animation for interpretation.
- `Ga2O3_bilayer_MD_50K_motion_zoomout.gif`: zoomed-out view with displacements
  amplified x3 for visibility.
- `Ga2O3_bilayer_MD_50K_motion.gif`: closer view with displacements amplified
  x5 for visibility.
- `make_bilayer_md_gif.py`: renderer used to generate the animations.

The x3 and x5 files are visualization aids only. Their apparent atomic
separations are not physical bond-breaking amplitudes.

Example:

```bash
python3 make_bilayer_md_gif.py trajectory_50K.lammpstrj output.gif \
  --amplify 1 --window 34 --fps 6
```

