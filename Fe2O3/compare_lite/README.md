# Literature comparison for Fe2O3 V100, V102, and V104 phonons

This folder compares the present harmonic PBE+U finite-displacement results
with representative experimental and first-principles results for bulk
alpha-Fe2O3 (hematite).

## Included comparisons

- `raman_mode_comparison.*`: seven Raman-active Gamma-point modes
  (`2 A1g + 5 Eg`) from experiment, a published bulk DFT calculation, V100,
  V102, and V104.
- `key_characteristics_comparison.*`: maximum phonon energy, the DOS gap near
  50 meV, Fe/O character in low/high-energy modes, and room-temperature
  entropy.
- `raman_modes.csv` and `key_characteristics.csv`: numeric values plotted.
- `plot_literature_comparison.py`: reproducible plotting script.

Run the plotting script from any directory with:

```bash
MPLCONFIGDIR=/tmp/qe-zeroK-mpl python Fe2O3/compare_lite/plot_literature_comparison.py
```

## Literature sources

1. J. L. Shelton and K. E. Knowles, *Polaronic Optical Transitions in
   Hematite (alpha-Fe2O3) Revealed by First-Principles Electron-Phonon
   Coupling*, J. Chem. Phys. 157, 174703 (2022),
   https://doi.org/10.1063/5.0116233 . The paper reports a DFT+U+J phonon
   bandwidth near 81 meV, a phonon-DOS gap just below 50 meV, Fe-dominated
   modes below 30 meV, and O-dominated modes above 50 meV.

2. C. P. Marshall, W. J. B. Dufresne, and C. J. Rufledt, *Polarized Raman
   spectra of hematite and assignment of external modes*, J. Raman Spectrosc.
   51, 1522-1529 (2020), https://doi.org/10.1002/jrs.5824 . This study resolves
   the external Eg modes near 245 and 294 cm-1.

3. M. Hanesch, *Raman spectroscopy of iron oxides and (oxy)hydroxides at low
   laser power and possible applications in environmental magnetic studies*,
   Geophys. J. Int. 177, 941-948 (2009),
   https://doi.org/10.1111/j.1365-246X.2009.04122.x . Representative hematite
   bands occur near 225, 245, 291, 411, 500, and 611 cm-1; the closely spaced
   modes near 293-300 cm-1 may not be resolved in every spectrum.

4. C. Bacaksiz, M. Yagmurcukardes, F. M. Peeters, and M. V. Milosevic,
   *Hematite at its thinnest limit*, 2D Materials 7, 025029 (2020),
   https://doi.org/10.1088/2053-1583/ab6d79 . The bulk DFT
   Raman frequencies used here are 220.0 and 484.3 cm-1 (A1g) and 244.0,
   280.0, 292.1, 410.4, and 607.9 cm-1 (Eg).

5. C. L. Snow, Q. Shi, J. Boerio-Goates, and B. F. Woodfield, *Heat capacity,
   third-law entropy, and low-temperature physical behavior of bulk hematite
   (alpha-Fe2O3)*, J. Chem. Thermodynamics 42, 1136-1141 (2010),
   https://doi.org/10.1016/j.jct.2010.04.010 . The reported
   standard molar entropy at 298.15 K is about 87.32 +/- 2 J mol-1 K-1.

## Interpretation and limitations

- The experimental Raman column uses representative bulk peak positions;
  temperature, sample strain, crystallinity, and peak fitting can shift them
  by a few cm-1.
- The literature DOS-gap location is stated approximately in the paper, so it
  is shown as a shaded 48-50 meV guide rather than as an exact interval.
- V104 no longer has a strictly zero-DOS interval near 50 meV. Its finite-DOS
  minimum is plotted at 46.287 meV and is about 2.02% of its maximum DOS.
- The calculated entropies are harmonic constant-volume values at 300 K. The
  calorimetric value is an experimental standard molar entropy at 298.15 K.
- The present calculations do not include non-analytical corrections. This
  matters for polar IR-active LO-TO splitting near Gamma, but not directly for
  the Raman-active gerade modes compared here.
