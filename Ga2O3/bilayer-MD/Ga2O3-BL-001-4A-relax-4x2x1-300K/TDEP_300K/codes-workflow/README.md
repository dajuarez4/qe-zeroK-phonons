# TDEP 300 K phonon workflow

Run this workflow from any directory:

```bash
./TDEP_300K/codes-workflow/update_phonons.sh
```

It first rebuilds the `TDEP_300K` inputs from the current
`Ga2O3-BL-001-4x2x1.md.out`, so stale prepared frames are not reused. It then
refits the second-order force constants, copies them to `infile.forceconstant`,
updates the phonon dispersion, and regenerates the diagnostic plot. Before
fitting, it reports both the raw QE MD-step count and the number of complete
position--force configurations prepared for TDEP. It also verifies that the
metadata, positions, forces, and statistics contain matching step counts.

Defaults are a 5.5 Angstrom cutoff and 300 K. Override them if needed:

```bash
TDEP_CUTOFF=6.0 TDEP_TEMPERATURE=300 \
  ./TDEP_300K/codes-workflow/update_phonons.sh
```

If TDEP is installed elsewhere, set its build directory:

```bash
TDEP_BIN_DIR=/path/to/tdep/build/src \
  ./TDEP_300K/codes-workflow/update_phonons.sh
```
