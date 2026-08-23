# Decorrelated TDEP stride and cutoff comparison

This comparison uses the 618-frame, 597.9768 fs Quantum ESPRESSO trajectory.
Samples are anchored at the newest frame so that stride 20 is nested within
stride 10 and both cover nearly the full time interval.

## Sampling

- Full baseline: 618 frames, 0.9676 fs spacing, 5.5 Angstrom cutoff
- Stride 10: frames 8, 18, ..., 618; 62 frames at 9.676 fs spacing
- Stride 20: frames 18, 38, ..., 618; 31 frames at 19.352 fs spacing
- Tested second-order cutoffs: 5.0, 5.5, and 6.0 Angstrom

## Results

| Sample | Cutoff (Angstrom) | Frames | Minimum (THz) | Maximum (THz) | Force RMS (eV/Angstrom) | R2 | sigma_A |
|---|---:|---:|---:|---:|---:|---:|---:|
| Full | 5.5 | 618 | -0.710 | 24.631 | 0.366 | 0.882 | 0.343 |
| Stride 10 | 5.0 | 62 | -0.761 | 24.532 | 0.376 | 0.876 | 0.352 |
| Stride 10 | 5.5 | 62 | -0.694 | 24.618 | 0.367 | 0.882 | 0.343 |
| Stride 10 | 6.0 | 62 | -0.694 | 24.618 | 0.367 | 0.882 | 0.343 |
| Stride 20 | 5.0 | 31 | -0.619 | 24.633 | 0.379 | 0.879 | 0.348 |
| Stride 20 | 5.5 | 31 | -0.695 | 24.694 | 0.369 | 0.885 | 0.339 |
| Stride 20 | 6.0 | 31 | -0.695 | 24.694 | 0.369 | 0.885 | 0.339 |

Stride 10 reproduces the full 5.5 Angstrom result closely while using about one
tenth as many configurations. This confirms that adjacent MD frames carry
substantial redundant information. Stride 20 is still overdetermined but has
only 6.4 force equations per second-order parameter at the 5.5 Angstrom
cutoff, and its frequencies show greater sampling variation.

The 5.0 Angstrom cutoff now makes the lowest mode more negative and generally
worsens the force residual, R2, and anharmonicity measure. The earlier apparent
improvement at 538 frames was therefore not robust to additional sampling.
The 5.5 and 6.0 Angstrom fits remain numerically
identical and contain the same 1,125 force-constant parameters, showing that
there is no additional interaction shell between those cutoffs for this cell.

None of the tested combinations removes the shallow imaginary branch. The
most efficient representative choice is stride 10 with a 5.5 Angstrom cutoff.
The next independent diagnostic should exclude the early thermal transient
and compare an equilibrated late-time window.

Machine-readable results are in `tdep_stride_cutoff_comparison.csv` and
`tdep_stride_cutoff_comparison.json`. The dispersion panels are in
`tdep_stride_cutoff_dispersion_comparison.png`, and the sensitivity metrics are
in `tdep_stride_cutoff_metrics.png`.
