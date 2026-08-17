# Ga2O3 bilayer: 40x40x1 MD and TDEP at 10 K

This directory contains the 16,000-atom tabGAP molecular-dynamics and
second-order TDEP phonon tests performed at 10 K. The TDEP real-space cutoff
is 5.5 Angstrom and the dispersion path is Gamma-X-S-Y-Gamma.

## Results

| Calculation | MD sampling | TDEP configurations | Fit R2 | Minimum frequency (THz) | Y minimum (THz) |
|---|---:|---:|---:|---:|---:|
| Original | 10 ps | 5 | 0.810 | -0.336 | -0.336 |
| Independent extended window | 25 ps | 5, selected at 0/6/12/19/25 ps | 0.760 | -0.195 | +0.154 |

The extended trajectory contains 26 snapshots separated by 1 ps. Fitting all
26 simultaneously exceeded the available WSL memory because TDEP attempted a
9.88 GB allocation. The reported extended fit therefore uses five evenly
spaced snapshots covering the full 25 ps window.

The Y-point frequency changes sign between the independent samples. The
remaining extended-fit minimum lies between Y and Gamma and its eigenvector is
95.2% out of plane; the Gamma soft mode is 98.2% out of plane. These results
indicate a sampling-sensitive flexural mode rather than a robust Y-point
structural instability.

## Contents

- `original_10ps_5frames/`: original trajectory, inputs, logs, fitted force
  constants, phonon dispersion, and DOS.
- `extended_25ps_5frame_fit/`: 25 ps trajectory and independent full-window
  TDEP fit. The HDF5 dispersion contains eigenvectors used for polarization
  analysis.
- `phonon_sampling_convergence_10K_40x40.png`: comparison with the 0 K finite-
  displacement result and both 10 K TDEP fits.
- `phonon_comparison_0K_FD_vs_10K_TDEP_40x40.png`: original 0 K/10 K comparison.

The tabGAP potential files are maintained elsewhere in the project workspace
and are referenced by relative paths in the LAMMPS inputs; they are not
duplicated here.
