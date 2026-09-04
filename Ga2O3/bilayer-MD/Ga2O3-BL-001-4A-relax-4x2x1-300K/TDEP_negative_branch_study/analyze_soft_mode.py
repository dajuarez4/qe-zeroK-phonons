#!/usr/bin/env python3
"""Extract and visualize the eigenvector of the most negative TDEP mode."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


STUDY = Path(__file__).resolve().parent
CASE = STUDY / "cases/baseline_all_5p5A"
HDF5 = CASE / "outfile.dispersion_relations.hdf5"
POSCAR = CASE / "infile.ucposcar"


def read_poscar(path: Path):
    lines = path.read_text().splitlines()
    scale = float(lines[1])
    lattice = np.array([[float(v) for v in lines[i].split()] for i in range(2, 5)]) * scale
    species = lines[5].split()
    counts = [int(v) for v in lines[6].split()]
    symbols = [symbol for symbol, count in zip(species, counts) for _ in range(count)]
    direct = np.array([[float(v) for v in lines[i].split()[:3]] for i in range(8, 8 + sum(counts))])
    return lattice, symbols, direct, direct @ lattice


def main() -> None:
    lattice, symbols, direct, cartesian = read_poscar(POSCAR)
    with h5py.File(HDF5, "r") as handle:
        frequencies = handle["frequencies"][:]
        q_vectors = handle["q_vector"][:]
        q_axis = handle["q_values"][:]
        minimum = np.unravel_index(np.argmin(frequencies), frequencies.shape)
        q_index, band_index = map(int, minimum)
        q_cart = q_vectors[q_index]
        # TDEP stores q in reciprocal Cartesian coordinates without 2π.
        q_fractional = q_cart @ lattice.T
        eigenvector = (
            handle["eigenvectors_re"][q_index, band_index]
            + 1j * handle["eigenvectors_im"][q_index, band_index]
        ).reshape(len(symbols), 3)
        site_projection = handle["site_projection_per_mode"][q_index, band_index]

    # A finite-q eigenvector has an arbitrary global complex phase. Rotate it
    # so the largest component is real and positive for a reproducible snapshot.
    dominant_flat = int(np.argmax(np.abs(eigenvector)))
    phase_rotation = np.exp(-1j * np.angle(eigenvector.ravel()[dominant_flat]))
    eigenvector *= phase_rotation
    dominant_atom = dominant_flat // 3
    real_snapshot = eigenvector.real
    amplitudes = np.linalg.norm(eigenvector, axis=1)
    direction_weights = np.sum(np.abs(eigenvector) ** 2, axis=0)
    direction_weights /= direction_weights.sum()

    rows = []
    for atom, (symbol, position, vector, projection) in enumerate(
        zip(symbols, cartesian, eigenvector, site_projection), start=1
    ):
        rows.append({
            "atom": atom, "species": symbol,
            "x_A": position[0], "y_A": position[1], "z_A": position[2],
            "eigenvector_x_real": vector[0].real, "eigenvector_x_imag": vector[0].imag,
            "eigenvector_y_real": vector[1].real, "eigenvector_y_imag": vector[1].imag,
            "eigenvector_z_real": vector[2].real, "eigenvector_z_imag": vector[2].imag,
            "amplitude": np.linalg.norm(vector), "site_projection": projection,
        })
    with (STUDY / "soft_mode_eigenvector.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    colors = {"Ga": "#3b78b4", "O": "#d84a3a"}
    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax_top = fig.add_subplot(grid[0, :])
    ax_side = fig.add_subplot(grid[1, 0])
    ax_part = fig.add_subplot(grid[1, 1])

    # Real-space snapshot over repeated unit cells. The Bloch phase shows the
    # modulation associated with q, rather than pretending this is a Γ mode.
    arrow_scale = 1.2 / max(np.linalg.norm(real_snapshot, axis=1).max(), 1e-12)
    for ia in range(4):
        for ib in range(5):
            translation_frac = np.array([ia, ib, 0.0])
            translation_cart = translation_frac @ lattice
            bloch = np.exp(2j * np.pi * np.dot(q_fractional, translation_frac))
            displacement = np.real(eigenvector * bloch) * arrow_scale
            positions = cartesian + translation_cart
            for index, symbol in enumerate(symbols):
                ax_top.scatter(positions[index, 0], positions[index, 1], s=30,
                               color=colors[symbol], edgecolor="white", linewidth=.4, zorder=2)
                ax_top.arrow(positions[index, 0], positions[index, 1],
                             displacement[index, 0], displacement[index, 1],
                             color="#222222", width=.012, head_width=.12,
                             length_includes_head=True, alpha=.85, zorder=3)
    ax_top.set_aspect("equal"); ax_top.set_xlabel("x (Å)"); ax_top.set_ylabel("y (Å)")
    ax_top.set_title("Top view: real-space Bloch displacement snapshot (4×5 repeated cells)")
    ax_top.grid(alpha=.12)

    for index, symbol in enumerate(symbols):
        ax_side.scatter(cartesian[index, 0], cartesian[index, 2], s=90,
                        color=colors[symbol], edgecolor="white", zorder=2)
        ax_side.arrow(cartesian[index, 0], cartesian[index, 2],
                      real_snapshot[index, 0] * arrow_scale,
                      real_snapshot[index, 2] * arrow_scale,
                      color="#222222", width=.018, head_width=.16,
                      length_includes_head=True, zorder=3)
        ax_side.text(cartesian[index, 0], cartesian[index, 2] + .22,
                     f"{symbol}{index+1}", fontsize=8, ha="center")
    ax_side.set_aspect("equal"); ax_side.set_xlabel("x (Å)"); ax_side.set_ylabel("z (Å)")
    ax_side.set_title("Unit-cell side view (phase-selected real component)"); ax_side.grid(alpha=.15)

    labels = [f"{symbol}{index}" for index, symbol in enumerate(symbols, 1)]
    bar_colors = [colors[symbol] for symbol in symbols]
    ax_part.bar(labels, site_projection * 100, color=bar_colors)
    ax_part.set_ylabel("Mode participation (%)"); ax_part.set_title("Atomic participation in unstable mode")
    ax_part.tick_params(axis="x", rotation=45); ax_part.grid(axis="y", alpha=.2)

    frequency = float(frequencies[minimum])
    fig.suptitle(
        f"Most unstable TDEP mode: {frequency:.3f} THz · band {band_index+1} · "
        f"q = ({q_fractional[0]:.3f}, {q_fractional[1]:.3f}, {q_fractional[2]:.3f})",
        fontsize=17, weight="bold",
    )
    figure_path = STUDY / "04_soft_mode_eigenvector.png"
    fig.savefig(figure_path, dpi=190); plt.close(fig)

    summary = {
        "frequency_THz": frequency,
        "q_index_zero_based": q_index,
        "band_one_based": band_index + 1,
        "q_path_coordinate": float(q_axis[q_index]),
        "q_cartesian_reciprocal_Ainv_without_2pi": q_cart.tolist(),
        "q_fractional": q_fractional.tolist(),
        "dominant_atom_one_based": dominant_atom + 1,
        "dominant_species": symbols[dominant_atom],
        "dominant_atom_site_projection": float(site_projection[dominant_atom]),
        "polarization_weight_x": float(direction_weights[0]),
        "polarization_weight_y": float(direction_weights[1]),
        "polarization_weight_z": float(direction_weights[2]),
        "interpretation": (
            "The mode is strongly localized on the dominant atom if its site projection "
            "is near one. Such localization is not a uniform bilayer sliding mode and "
            "should prompt inspection of that atom's local environment and reference forces."
        ),
    }
    (STUDY / "soft_mode_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Saved {figure_path}")


if __name__ == "__main__":
    main()
