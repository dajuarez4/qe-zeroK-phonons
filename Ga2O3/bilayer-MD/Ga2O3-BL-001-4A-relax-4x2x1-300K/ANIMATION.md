# Quantum ESPRESSO MD animation at 300 K

`Ga2O3_Bilayer_QE_MD_300K_actual_motion.gif` visualizes the 80-atom bilayer
trajectory stored in `Ga2O3-BL-001-4x2x1.md.out`.

The animation uses 60 uniformly sampled configurations from the 830 complete
`ATOMIC_POSITIONS (crystal)` blocks. Displacements are shown at their physical
magnitude (x1). Periodic-boundary jumps are unwrapped, center-of-mass drift is
removed, and Ga-O bonds are recalculated for every rendered frame.

Reproduce it with:

```bash
python3 make_qe_md_gif.py Ga2O3-BL-001-4x2x1.md.out \
  Ga2O3_Bilayer_QE_MD_300K_actual_motion.gif --max-frames 60 --fps 8
```
