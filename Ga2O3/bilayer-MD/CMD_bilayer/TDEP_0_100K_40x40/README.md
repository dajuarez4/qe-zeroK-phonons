# Ga2O3 bilayer: MD/TDEP phonons from 0 to 100 K

## Protocol

- 40x40x1 bilayer supercell (16,000 atoms).
- tabGAP Ga-O potential; 1 fs timestep.
- 5 ps NVT equilibration and 25 ps NVT production at each finite temperature.
- 26 stored configurations per trajectory; five frames at 0, 6, 12, 19, and 25 ps used per TDEP fit because of the WSL memory ceiling.
- Second-order TDEP cutoff 5.5 Angstrom; Gamma-X-S-Y-Gamma path.

## Results

| T (K) | measured T (K) | Gamma | X | S | Y | global min | R2 | z weight | mean in-plane P (GPa) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.000 | -0.022 | +2.247 | +2.326 | +0.591 | -0.036 | nan | nan% | +nan |
| 10 | 10.025 | -0.107 | +2.005 | +2.047 | +0.154 | -0.195 | 0.760 | 95.2% | -0.170 |
| 20 | 19.995 | -0.098 | +2.126 | +2.122 | +0.092 | -0.223 | 0.749 | 94.2% | -0.212 |
| 30 | 30.084 | -0.135 | +1.541 | +1.560 | -0.127 | -0.227 | 0.703 | 93.6% | -0.277 |
| 40 | 40.031 | -0.086 | +1.055 | +1.042 | -0.159 | -0.228 | 0.684 | 93.5% | -0.327 |
| 50 | 49.916 | -0.140 | +0.521 | +0.472 | -0.260 | -0.273 | 0.538 | 90.5% | -0.350 |
| 60 | 60.101 | -0.351 | +0.348 | +0.409 | -0.342 | -0.351 | 0.461 | 76.0% | -0.362 |
| 70 | 69.921 | -0.267 | +0.280 | +0.255 | -0.253 | -0.310 | 0.350 | 50.6% | -0.357 |
| 80 | 80.016 | -0.291 | +0.199 | +0.307 | -0.348 | -0.372 | 0.325 | 59.0% | -0.358 |
| 90 | 90.138 | -0.243 | +0.348 | +0.340 | -0.411 | -0.411 | 0.302 | 84.5% | -0.367 |
| 100 | 100.080 | -0.286 | +0.345 | +0.257 | -0.393 | -0.445 | 0.297 | 75.6% | -0.365 |

Frequencies are in THz. A negative TDEP frequency denotes an imaginary harmonic mode of the fitted effective force constants; it does not by itself establish a static structural instability.

## Static-mode verification

For the 50 K Gamma eigenvector, 17 tabGAP single-point structures spanning Q = -0.20 to +0.20 Angstrom were evaluated. The local quadratic coefficient is +6.371 eV/Angstrom^2 per 10-atom cell and the fitted minimum is Q = -0.0027 Angstrom. The positive curvature and absence of a double well demonstrate static stability along the Gamma coordinate despite the negative TDEP frequency.

## Interpretation

The low branches are predominantly out of plane and therefore flexural. Their fitted frequencies are affected by finite sampling, the decreasing quality of a purely second-order model, residual anisotropic tensile stress from the fixed in-plane cell, and the known limited phonon accuracy of the general-purpose tabGAP near the lowest acoustic branch. These frequencies should not be reported as proof of a phase instability without direct DFT validation or a finite-temperature cell/stress treatment.
