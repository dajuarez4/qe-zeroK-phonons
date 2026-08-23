#!/opt/anaconda3/bin/python
"""Build a diagnostic, underdetermined phonon model from disp-001..011.

This is not a substitute for the production 240-displacement calculation.
Only five complete central-difference directions are available. Unsampled
force-constant blocks are completed with the minimum-norm (zero) solution,
then translation symmetry, Hessian symmetry, and the acoustic sum rule are
imposed.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from phonopy import load
from phonopy.file_IO import write_FORCE_CONSTANTS
from phonopy.interface.qe import parse_set_of_forces
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from phonopy.units import Bohr, Rydberg


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "partial_11_phonons"
OUT.mkdir(exist_ok=True)

phonon = load(str(ROOT / "phonopy_disp.yaml"))
dataset = phonon.dataset["first_atoms"]
natom = len(phonon.supercell)
ndof = natom * 3

files = [
    ROOT / "displacements" / f"disp-{i:03d}" / "Al69_Ga2O3.fd.scf.out"
    for i in range(1, 12)
]
forces = np.asarray(
    parse_set_of_forces(natom, [str(path) for path in files], verbose=False)
)
if forces.shape != (11, natom, 3):
    raise RuntimeError(f"Expected 11 force sets of shape ({natom}, 3), got {forces.shape}")

# The QE parser returns Ry/bohr. Convert to eV/angstrom to match Phonopy's
# displacement and force-constant conventions.
forces *= Rydberg / Bohr

# Complete central-difference pairs. disp-011 has no negative partner yet and
# is deliberately excluded instead of introducing a biased one-sided column.
groups = {
    0: [(0, 1), (2, 3), (4, 5)],
    2: [(6, 7), (8, 9)],
}

raw = np.zeros((ndof, ndof), dtype="double")
sampled_atoms = []
direction_ranks = {}
for displaced_atom, pairs in groups.items():
    delta_u = np.asarray(
        [
            np.asarray(dataset[ip]["displacement"])
            - np.asarray(dataset[im]["displacement"])
            for ip, im in pairs
        ]
    )
    delta_f = np.asarray([forces[ip] - forces[im] for ip, im in pairs])
    # delta_f = -Phi delta_u. The pseudoinverse is exact for atom 1 (rank 3)
    # and gives the minimum-norm completion for atom 3 (rank 2).
    phi = (np.linalg.pinv(delta_u) @ (-delta_f.reshape(len(pairs), ndof))).T
    raw[:, 3 * displaced_atom : 3 * displaced_atom + 3] = phi
    sampled_atoms.append(displaced_atom)
    direction_ranks[displaced_atom + 1] = int(np.linalg.matrix_rank(delta_u))

# The 1x2x1 supercell has one nontrivial pure translation. Populate the
# translational image of each sampled atom and its corresponding response.
translation = phonon.symmetry.atomic_permutations[1]
known_atoms = set(sampled_atoms)
for displaced_atom in sampled_atoms:
    translated_atom = int(translation[displaced_atom])
    block = raw[:, 3 * displaced_atom : 3 * displaced_atom + 3].reshape(natom, 3, 3)
    translated_block = np.zeros_like(block)
    translated_block[translation] = block
    raw[:, 3 * translated_atom : 3 * translated_atom + 3] = translated_block.reshape(ndof, 3)
    known_atoms.add(translated_atom)

# Minimum-norm symmetric completion. If both transposed entries were directly
# inferred, average them; if one was inferred, copy it; otherwise leave zero.
known_columns = np.zeros(ndof, dtype=bool)
for atom in known_atoms:
    known_columns[3 * atom : 3 * atom + 3] = True

hessian = np.zeros_like(raw)
for a in range(ndof):
    for b in range(a, ndof):
        values = []
        if known_columns[b]:
            values.append(raw[a, b])
        if known_columns[a]:
            values.append(raw[b, a])
        if values:
            hessian[a, b] = hessian[b, a] = np.mean(values)

# Orthogonal projection removes the three uniform translations while
# retaining Hessian symmetry: H <- P H P.
translations = np.zeros((ndof, 3))
for atom in range(natom):
    translations[3 * atom : 3 * atom + 3, :] = np.eye(3)
projector = np.eye(ndof) - translations @ np.linalg.inv(
    translations.T @ translations
) @ translations.T
hessian = projector @ hessian @ projector
hessian = (hessian + hessian.T) / 2
force_constants = hessian.reshape(natom, 3, natom, 3).transpose(0, 2, 1, 3)

phonon.force_constants = force_constants
write_FORCE_CONSTANTS(force_constants, filename=str(OUT / "FORCE_CONSTANTS_PARTIAL_11"))

special_points = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [0.5, 0.5, 0.5],
    [0.0, 0.5, 0.5],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.0],
]
labels = ["Gamma", "A", "Z", "M", "L", "V", "Gamma"]
qpaths, connections = get_band_qpoints_and_path_connections(
    [special_points], npoints=81, rec_lattice=np.linalg.inv(phonon.primitive.cell)
)
phonon.run_band_structure(qpaths, path_connections=connections, labels=labels)
band = phonon.get_band_structure_dict()

distances = np.concatenate(band["distances"])
frequencies = np.concatenate(band["frequencies"])
np.savetxt(
    OUT / "partial_11_band_frequencies.csv",
    np.column_stack((distances, frequencies)),
    delimiter=",",
    header="distance," + ",".join(f"mode_{i + 1}_THz" for i in range(frequencies.shape[1])),
    comments="",
)

ticks = [band["distances"][0][0]]
for segment in band["distances"]:
    ticks.append(segment[-1])

fig, ax = plt.subplots(figsize=(8.2, 5.2))
ax.plot(distances, frequencies, color="tab:blue", lw=0.55, alpha=0.75)
ax.axhline(0, color="black", lw=0.7)
for tick in ticks:
    ax.axvline(tick, color="0.82", lw=0.6)
ax.set_xticks(ticks, [r"$\Gamma$", "A", "Z", "M", "L", "V", r"$\Gamma$"])
ax.set_xlim(distances[0], distances[-1])
ax.set_ylabel("Frequency (THz)")
ax.set_title("Partial 11-SCF diagnostic (underdetermined; not production phonons)")
fig.tight_layout()
fig.savefig(OUT / "partial_11_phonon_bands.png", dpi=220)
fig.savefig(OUT / "partial_11_phonon_bands.pdf")
plt.close(fig)

phonon.run_qpoints([[0, 0, 0]])
gamma = np.asarray(phonon.get_qpoints_dict()["frequencies"])[0]
matrix_rank = int(np.linalg.matrix_rank(hessian, tol=1e-8))
near_zero_gamma = int(np.count_nonzero(np.abs(gamma) < 1e-4))
summary = f"""PARTIAL 11-DISPLACEMENT PHONON DIAGNOSTIC

WARNING: This is an underdetermined minimum-norm completion, not a physical
or converged phonon calculation. Do not use it for stability conclusions,
thermodynamics, publication, or comparison with experiment.

Available QE outputs: 11 / 240
Complete central-difference pairs used: 5
Unpaired disp-011: excluded
Sampled primitive atoms (1-based): 1, 3
Independent direction rank for atom 1: {direction_ranks[1]} / 3
Independent direction rank for atom 3: {direction_ranks[3]} / 3
Hessian rank after constraints: {matrix_rank} / {ndof}
Gamma modes with |frequency| < 1e-4 THz: {near_zero_gamma} / {len(gamma)}
Frequency range on requested path: {frequencies.min():.6f} to {frequencies.max():.6f} THz
Gamma frequency range: {gamma.min():.6f} to {gamma.max():.6f} THz

Constraints imposed: 1x2x1 translation symmetry, Hessian symmetry, acoustic
sum rule. All otherwise unsampled force-constant information was assigned the
minimum-norm zero completion.
"""
(OUT / "README.txt").write_text(summary)
print(summary)
