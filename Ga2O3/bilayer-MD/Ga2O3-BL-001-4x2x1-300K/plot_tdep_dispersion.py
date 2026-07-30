#!/usr/bin/env python3
"""Plot the current TDEP dispersion and write compact numerical diagnostics."""

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


def main() -> None:
    root = Path(__file__).resolve().parent
    outdir = root / "TDEP_300K"
    data = np.loadtxt(outdir / "outfile.dispersion_relations")
    x_values = data[:, 0]
    frequencies = data[:, 1:]

    # The TDEP path contains four segments with 80 points per segment.
    tick_indices = [0, 79, 159, 239, 319]
    tick_positions = x_values[tick_indices]
    tick_labels = [r"$\Gamma$", "X", "S", "Y", r"$\Gamma$"]

    fig, ax = plt.subplots(figsize=(10.5, 6.6), constrained_layout=True)
    for position in tick_positions:
        ax.axvline(position, color="#aaaaaa", linewidth=0.7, zorder=0)
    for band in range(frequencies.shape[1]):
        ax.plot(x_values, frequencies[:, band], color="#245b93", linewidth=1.0)
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xlim(x_values[0], x_values[-1])
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_ylabel("Frequency (THz)")
    ax.set_title("Ga₂O₃ (001) 4×2×1 TDEP diagnostic (30 MD frames)")
    ax.grid(axis="y", alpha=0.2)
    fig.savefig(
        outdir / "Ga2O3-BL-001-4x2x1-TDEP-diagnostic.png",
        dpi=180,
    )
    plt.close(fig)

    summary = {
        "n_frames": 30,
        "trajectory_time_fs": 29.028,
        "mean_temperature_K": 518.3075248823,
        "minimum_frequency_THz": float(np.min(frequencies)),
        "maximum_frequency_THz": float(np.max(frequencies)),
        "negative_frequency_values_below_minus_1e-6_THz": int(
            np.sum(frequencies < -1.0e-6)
        ),
        "negative_frequency_values_below_minus_0p1_THz": int(
            np.sum(frequencies < -0.1)
        ),
        "total_frequency_values": int(frequencies.size),
    }
    (outdir / "tdep_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
