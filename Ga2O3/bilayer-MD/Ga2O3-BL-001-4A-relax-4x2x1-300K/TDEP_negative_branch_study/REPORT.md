# TDEP negative-branch study

This is an independent snapshot study of the 1017-frame TDEP dataset. The negative mode is robust across the tested sampling choices.

## Key results

- Baseline minimum: -1.050 THz with σA = 0.360.
- Late-half minimum: -1.749 THz (509 frames).
- Stride-20 minimum: -1.075 THz (51 frames).
- Cutoff minima: 4.5 Å → -1.082 THz, 5.0 Å → -0.880 THz, 5.5 Å → -1.050 THz.
- Temperature integrated autocorrelation estimate: 45.6 fs; approximate effective temperature samples: 10.8.

## Interpretation

Persistence in the late window and decorrelated fits argues against the branch being caused only by duplicated adjacent MD frames. Strong cutoff sensitivity would instead implicate the finite supercell/FC2 truncation. Persistence across both tests makes a real soft interlayer/flexural mode or a structural-reference/sum-rule issue more plausible. Absolute-valued frequencies must not be used for this diagnosis because they hide the sign.

## Unstable-mode eigenvector

The baseline HDF5 eigenvectors identify the minimum as band 1 at fractional
q = (0, 0.3924, 0), on the Y–Γ segment. The mode is highly localized rather
than a uniform bilayer translation:

- O6 contributes 92.42% of the total site projection.
- Cartesian polarization weights are 65.89% x, 10.50% y, and 23.62% z.
- The visualization is `04_soft_mode_eigenvector.png`; numerical complex
  components are in `soft_mode_eigenvector.csv`.

This localization makes a simple rigid interlayer-sliding interpretation
unlikely. It instead points toward O6's local bonding/reference environment,
an inadequately represented local anharmonic motion, or a force-constant/
constraint artifact centered on that site. The next structural check should
compare O6 bond lengths and force statistics against the other oxygen sites
through the trajectory.
