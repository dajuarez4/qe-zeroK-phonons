# QE Zero-K Phonons

A general Quantum ESPRESSO workflow for calculating harmonic phonons and phonon dispersions at 0 K for crystalline materials.

This repository is designed as a reusable template for different materials. The workflow starts from a converged ground-state calculation and produces phonon frequencies, force constants, and phonon dispersion plots.

---

## Workflow Overview

The general workflow is:

```text
SCF calculation → phonon q-grid → force constants → phonon dispersion → plot
```

Using Quantum ESPRESSO executables:

```text
pw.x      → ground-state SCF calculation
ph.x      → dynamical matrices on a q-point grid
q2r.x     → real-space interatomic force constants
matdyn.x  → phonon frequencies along a selected q-path
```

---

## Repository Structure

```text
qe-zerok-phonons/
├── README.md
├── .gitignore
├── examples/
│   └── material_example/
│       ├── scf.in
│       ├── ph_grid.in
│       ├── q2r.in
│       ├── matdyn.in
│       └── run_all.sh
├── scripts/
│   ├── run_qe_phonons.sh
│   └── plot_phonons.py
├── templates/
│   ├── scf.template.in
│   ├── ph_grid.template.in
│   ├── q2r.template.in
│   └── matdyn.template.in
├── pseudos/
│   └── README.md
└── results/
    └── README.md
```

---

## 1. Ground-State Calculation

The first step is a self-consistent field calculation using `pw.x`.

Example:

```bash
pw.x < scf.in > scf.out
```

The SCF calculation must be well converged before calculating phonons. Important convergence parameters include:

* plane-wave cutoff energy
* charge-density cutoff energy
* k-point mesh
* total energy convergence threshold
* pseudopotential choice
* crystal structure and lattice parameters

For phonons at 0 K, the structure should be relaxed or taken from a reliable experimental/theoretical reference.

---

## 2. Phonon Calculation on a q-Grid

After the SCF calculation, phonons are calculated using `ph.x`.

Example input:

```fortran
Phonon q-grid calculation
&inputph
  prefix  = 'material',
  outdir  = './tmp',
  fildyn  = 'material.dyn',
  tr2_ph  = 1.0d-14,

  ldisp = .true.,
  nq1 = 2,
  nq2 = 2,
  nq3 = 2,
/
```

Run:

```bash
ph.x < ph_grid.in > ph_grid.out
```

This produces dynamical matrices on a uniform q-point grid.

Example output files:

```text
material.dyn0
material.dyn1
material.dyn2
...
```

A small grid such as `2×2×2` is useful for testing. Final calculations should use a converged q-grid, such as `4×4×4`, `6×6×6`, or larger depending on the material.

---

## 3. Convert Dynamical Matrices to Force Constants

The dynamical matrices are converted to real-space interatomic force constants using `q2r.x`.

Example input:

```fortran
&input
  fildyn = 'material.dyn',
  flfrc  = 'material.fc',
  zasr   = 'crystal',
/
```

Run:

```bash
q2r.x < q2r.in > q2r.out
```

Main output:

```text
material.fc
```

This file contains the real-space force constants.

---

## 4. Calculate Phonon Frequencies Along a q-Path

The phonon dispersion is obtained by calculating phonon frequencies along a chosen high-symmetry path using `matdyn.x`.

Example input:

```fortran
&input
  asr = 'crystal',
  flfrc = 'material.fc',
  flfrq = 'material.freq',
  q_in_band_form = .true.,
  q_in_cryst_coord = .true.,
/
8
  0.000 0.000 0.000 40  ! Gamma
  0.500 0.000 0.000 40  ! M
  0.333 0.333 0.000 40  ! K
  0.000 0.000 0.000 40  ! Gamma
  0.000 0.000 0.500 40  ! A
  0.500 0.000 0.500 40  ! L
  0.333 0.333 0.500 40  ! H
  0.000 0.000 0.500 1   ! A
```

Run:

```bash
matdyn.x < matdyn.in > matdyn.out
```

Main outputs:

```text
material.freq
material.freq.gp
```

The `.freq.gp` file can be used directly for plotting.

---

## 5. Plot the Phonon Dispersion

A Python script can be used to plot the phonon dispersion from the `material.freq.gp` file.

Example:

```bash
python ../../scripts/plot_phonons.py material.freq.gp material_phonon_dispersion.png
```

The final result is:

```text
material_phonon_dispersion.png
```

---

## Example Full Run

Inside one material folder:

```bash
bash run_all.sh
```

Example `run_all.sh`:

```bash
#!/bin/bash
set -e

PREFIX="material"

echo "Running SCF calculation..."
pw.x < scf.in > scf.out

echo "Running phonon q-grid calculation..."
ph.x < ph_grid.in > ph_grid.out

echo "Converting dynamical matrices to force constants..."
q2r.x < q2r.in > q2r.out

echo "Calculating phonon dispersion..."
matdyn.x < matdyn.in > matdyn.out

echo "Plotting phonon dispersion..."
python ../../scripts/plot_phonons.py "${PREFIX}.freq.gp" "${PREFIX}_phonon_dispersion.png"

echo "Done."
```

---

## Python Plot Script

Example `plot_phonons.py`:

```python
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

if len(sys.argv) < 3:
    print("Usage: python plot_phonons.py input.freq.gp output.png")
    sys.exit(1)

input_file = Path(sys.argv[1])
output_file = Path(sys.argv[2])

if not input_file.exists():
    print(f"Error: {input_file} not found.")
    sys.exit(1)

data = np.loadtxt(input_file)

x = data[:, 0]
freqs = data[:, 1:]

plt.figure(figsize=(8, 6))

for i in range(freqs.shape[1]):
    plt.plot(x, freqs[:, i], linewidth=1.0)

plt.axhline(0, linestyle="--", linewidth=0.8)

plt.xlabel("Wave vector")
plt.ylabel("Frequency (cm$^{-1}$)")
plt.title("Phonon dispersion")
plt.tight_layout()
plt.savefig(output_file, dpi=300)

print(f"Saved plot to {output_file}")
```

---

## Important Notes

### q-Grid Convergence

The phonon q-grid must be tested for convergence. A small grid is acceptable for testing, but final results require a denser grid.

Example test grids:

```text
2×2×2
4×4×4
6×6×6
```

### Imaginary Frequencies

Negative frequencies in the phonon dispersion usually indicate imaginary modes. These can mean:

* the structure is dynamically unstable
* the structure was not properly relaxed
* the calculation is not converged
* the q-grid is too coarse
* the pseudopotential or calculation settings need checking

### High-Symmetry Path

The q-path depends on the crystal structure. For every new material, the `matdyn.in` path should be changed according to the Brillouin zone of that material.

### Polar Materials

For polar materials, non-analytical corrections may be needed to properly describe LO-TO splitting near Gamma. This requires additional quantities such as Born effective charges and the dielectric tensor.

### Metals

For metallic systems, use appropriate smearing in the SCF calculation and carefully converge the k-point mesh.

---

## Main Files Produced

```text
scf.out              SCF output
ph_grid.out          phonon q-grid output
material.dyn*        dynamical matrices
q2r.out              q2r output
material.fc          real-space force constants
matdyn.out           matdyn output
material.freq        phonon frequencies
material.freq.gp     plottable phonon dispersion data
*.png                phonon dispersion figure
```

---

## Goal of This Repository

The goal of this repository is to provide a clean and reusable 0 K phonon workflow using Quantum ESPRESSO. It can be adapted to different materials by changing the structure, pseudopotentials, convergence settings, q-grid, and high-symmetry path.
