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
- `gruneisen_parameters.*`, `gruneisen_temperature.csv`,
  `gruneisen_summary.csv`, and `raman_mode_gruneisen.csv`: mode and
  heat-capacity-weighted Gruneisen parameters from V100--V104.
- `volume_dependent_thermodynamics.csv`: harmonic ZPE, free energy, entropy,
  and heat capacity at each completed phonon volume.
- `preliminary_qha_thermodynamics.*` and `preliminary_qha_pressure_grid.*`:
  equilibrium volume, thermal expansion, bulk modulus, Cp/Cv, and a small
  pressure-temperature grid obtained by combining the nine-point static EOS
  with quadratic phonon free-energy fits.
- `mean_square_displacements.*` and `displacement_tensors_300K.csv`: Fe/O
  mean-square displacements and atom-resolved 300 K displacement tensors.
- `acoustic_velocities.*`: near-Gamma phase velocities and group-velocity
  projections along the three symmetry-equivalent primitive reciprocal axes.
- `species_projected_thermodynamics.*`: Fe- and O-projected ZPE, Cv, and
  entropy derived from the projected DOS.
- `extended_literature_comparison.*`: entropy, 300 K Cp, and bulk-modulus
  benchmarks.
- `compute_extended_properties.py`: reproducible extended-analysis script.

Run the plotting script from any directory with:

```bash
MPLCONFIGDIR=/tmp/qe-zeroK-mpl python Fe2O3/compare_lite/plot_literature_comparison.py
MPLCONFIGDIR=/tmp/qe-zeroK-mpl /opt/anaconda3/bin/python Fe2O3/compare_lite/compute_extended_properties.py
```

## Main extended results

- The static third-order Birch-Murnaghan fit to all nine completed EOS points
  gives V0 = 106.262 A3 per 10-atom cell, B0 = 187.57 GPa, and B0' = 4.234.
- The median mode Gruneisen parameter is 1.328.  The Cv-weighted value is
  1.298 at 300 K and rises to 1.351 at 1000 K.  Raman-mode values range from
  0.708 to 1.964.
- At V102 and 300 K, harmonic S = 90.406 J mol-1 K-1 and Cv = 97.752
  J mol-1 K-1 per Fe2O3.  The projected-DOS partition assigns about 51.0
  J mol-1 K-1 of the entropy to Fe and 38.9 J mol-1 K-1 to O.
- The preliminary QHA result at 300 K is V = 107.567 A3 per 10-atom cell,
  linear alpha = 9.18 x 10-6 K-1, BT = 175.76 GPa, and Cp = 102.11
  J mol-1 K-1.  These four values are extrapolative, as explained below.
- At V102 and 300 K, mean total MSD is 0.01343 A2 for Fe and 0.01586 A2 for
  O.  Near-Gamma acoustic speeds along the equivalent reciprocal primitive
  axes are about 3.28, 4.70, and 8.38 km/s at V102.  These are sorted acoustic
  branches, not polarization-resolved elastic constants.

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

6. NIST Chemistry WebBook, SRD 69, *Hematite (Fe2O3)*,
   https://webbook.nist.gov/cgi/cbook.cgi?ID=C1317608&Mask=2&Plot=on&Type=JANAFS .
   The Chase (1998) Shomate coefficients give Cp = 104.155 J mol-1 K-1 at
   300 K.

7. Y. Zou, P. Wang, Y. Li, H. Chen, C. Zhou, and T. Irifune, *Unveiling
   pressure-induced anomalous shear behavior and thermoelasticity of
   alpha-Fe2O3 hematite at high pressure*, iScience 28, 111905 (2025),
   https://doi.org/10.1016/j.isci.2025.111905 . This article summarizes an
   ambient ultrasonic aggregate value K0 = 206.6 GPa (with G0 = 91.0 GPa)
   from earlier Liebermann measurements; the comparison is not strictly
   like-for-like with a 0 K static isothermal EOS fit.

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
- Only V100, V102, and V104 have completed phonons, whereas the static EOS
  minimum is near V106.26.  Consequently, every zero-pressure QHA equilibrium
  volume lies outside the sampled phonon range.  All equilibrium QHA outputs
  are explicitly flagged `phonon_volume_extrapolated=True`; they are useful
  trends, but V106 and preferably V108 phonons are needed for a defensible
  zero-pressure QHA fit.
- Three volumes determine a quadratic phonon free-energy curve exactly and
  provide no redundancy for testing curvature.  Mode matching uses sorted
  branch indices on a common 12x12x12 mesh, so isolated branch crossings can
  affect individual mode Gruneisen values; distribution averages are more
  robust.
- Acoustic velocities were converted from the QE/phonopy physical length unit
  (bohr) to km/s.  They should be checked against an acoustic-sum-rule and
  supercell-size convergence study before use as precision elastic data.
- Harmonic second-order force constants cannot determine phonon lifetimes or
  lattice thermal conductivity.  Those require third-order force constants
  (or an anharmonic MD/TDEP workflow), which were not inferred here.
