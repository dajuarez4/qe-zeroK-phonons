#!/usr/bin/env python3
"""Generate primitive-cell structures displaced along the 50 K Gamma soft mode."""

from pathlib import Path
import h5py
import numpy as np
from ase.io import read, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
atoms = read(ROOT / "tdep_50K_40x40/infile.ucposcar", format="vasp")

with h5py.File(ROOT / "tdep_50K_40x40/outfile.dispersion_relations.hdf5", "r") as h5:
    frequencies = h5["frequencies"][:]
    mode_index = int(np.argmin(frequencies[0]))
    frequency = float(frequencies[0, mode_index])
    vector = h5["eigenvectors_re"][0, mode_index] + 1j * h5["eigenvectors_im"][0, mode_index]

# Gamma eigenvectors are real up to an arbitrary complex phase. Rotate to the
# phase with the largest real norm, remove rigid translation, and normalize so
# Q is the largest displacement of any atom in Angstrom.
phase = np.exp(-1j * np.angle(vector[np.argmax(np.abs(vector))]))
mode = np.real(vector * phase).reshape(len(atoms), 3)
mode -= mode.mean(axis=0)
mode /= np.linalg.norm(mode, axis=1).max()

amplitudes = np.linspace(-0.20, 0.20, 17)
(HERE / "structures").mkdir(exist_ok=True)
reference = atoms.get_positions().copy()
for index, amplitude in enumerate(amplitudes):
    displaced = atoms.copy()
    displaced.set_positions(reference + amplitude * mode)
    write(HERE / "structures" / f"mode_{index:02d}.data", displaced,
          format="lammps-data", atom_style="atomic", specorder=["Ga", "O"])

np.savetxt(HERE / "mode_vector.dat", mode, header="dx dy dz; max atomic norm = 1")
np.savetxt(HERE / "amplitudes.dat", amplitudes, header="Q_Angstrom")
print(f"Gamma mode {mode_index}, TDEP frequency {frequency:.12f} THz")
print(f"Generated {len(amplitudes)} structures")
