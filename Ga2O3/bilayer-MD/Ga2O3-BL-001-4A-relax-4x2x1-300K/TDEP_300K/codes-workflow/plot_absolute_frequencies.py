#!/usr/bin/env python3
"""Plot TDEP phonons with negative frequencies reflected using |frequency|."""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


SCRIPT_DIR = Path(__file__).resolve().parent
TDEP_DIR = SCRIPT_DIR.parent
INPUT = TDEP_DIR / "outfile.dispersion_relations"
OUTPUT = TDEP_DIR / "Ga2O3-BL-001-4A-TDEP-absolute-frequencies.png"


def main() -> None:
    data = np.loadtxt(INPUT)
    x_values = data[:, 0]
    signed_frequencies = data[:, 1:]
    absolute_frequencies = np.abs(signed_frequencies)

    reference_path = TDEP_DIR / "reference_summary.json"
    frames = None
    if reference_path.exists():
        frames = json.loads(reference_path.read_text()).get("n_frames")

    points_per_segment = 80
    tick_indices = [0, 79, 159, 239, 319]
    if len(x_values) <= tick_indices[-1]:
        raise ValueError(f"Expected at least 320 q-points, found {len(x_values)}")
    tick_positions = x_values[tick_indices]

    fig, ax = plt.subplots(figsize=(10.5, 6.6), constrained_layout=True)
    for position in tick_positions:
        ax.axvline(position, color="#aaaaaa", linewidth=0.7, zorder=0)

    # Plot all absolute-valued bands in blue. Overlay only the portions that
    # originated below zero in red, preserving that physical information.
    for band in range(absolute_frequencies.shape[1]):
        y_abs = absolute_frequencies[:, band]
        ax.plot(x_values, y_abs, color="#245b93", linewidth=1.0)
        reflected = np.ma.masked_where(signed_frequencies[:, band] >= 0.0, y_abs)
        ax.plot(x_values, reflected, color="#c43c39", linewidth=1.7)

    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xlim(x_values[0], x_values[-1])
    ax.set_ylim(bottom=0.0)
    ax.set_xticks(tick_positions, [r"$\Gamma$", "X", "S", "Y", r"$\Gamma$"])
    ax.set_ylabel(r"Absolute frequency $|\nu|$ (THz)")
    frame_text = f" · {frames} MD frames" if frames is not None else ""
    ax.set_title(f"Ga₂O₃ 4 Å bilayer · TDEP absolute frequencies{frame_text}")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(
        handles=[
            Line2D([0], [0], color="#245b93", lw=1.2, label="Originally non-negative"),
            Line2D([0], [0], color="#c43c39", lw=1.7, label="Reflected from negative"),
        ],
        loc="upper right",
        frameon=False,
    )
    fig.savefig(OUTPUT, dpi=180)
    plt.close(fig)

    negative_values = signed_frequencies[signed_frequencies < 0.0]
    print(f"Saved: {OUTPUT}")
    print(f"Original minimum frequency: {signed_frequencies.min():.6f} THz")
    print(f"Negative values reflected: {negative_values.size}")


if __name__ == "__main__":
    main()
