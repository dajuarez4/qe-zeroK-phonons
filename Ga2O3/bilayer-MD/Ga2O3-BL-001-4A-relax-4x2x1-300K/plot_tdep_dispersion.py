#!/usr/bin/env python3
"""Plot the nine-frame 4 Angstrom TDEP diagnostic."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "TDEP_300K"


def main():
    data = np.loadtxt(OUTDIR / "outfile.dispersion_relations")
    x_values, frequencies = data[:, 0], data[:, 1:]
    reference = json.loads((OUTDIR / "reference_summary.json").read_text())
    fit_log = (OUTDIR / "extract_forceconstants.log").read_text(errors="replace")
    first_order = re.search(r"RMSE total:\s+([-+0-9.Ee]+)", fit_log)
    anharmonicity = re.search(
        r"^\s*second order:\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"<-- anharmonicity measure",
        fit_log,
        re.MULTILINE,
    )
    values = [float(value) for value in anharmonicity.groups()] if anharmonicity else [None] * 5
    tick_indices = [0, 79, 159, 239, 319]
    tick_positions = x_values[tick_indices]

    fig, ax = plt.subplots(figsize=(10.5, 6.6), constrained_layout=True)
    for position in tick_positions:
        ax.axvline(position, color="#aaaaaa", linewidth=0.7, zorder=0)
    for band in range(frequencies.shape[1]):
        ax.plot(x_values, frequencies[:, band], color="#245b93", linewidth=1.0)
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xlim(x_values[0], x_values[-1])
    ax.set_xticks(tick_positions, [r"$\Gamma$", "X", "S", "Y", r"$\Gamma$"])
    ax.set_ylabel("Frequency (THz)")
    ax.set_title(
        f"Ga₂O₃ 4 Å bilayer · TDEP diagnostic ({reference['n_frames']} MD frames)"
    )
    ax.grid(axis="y", alpha=0.2)
    fig.savefig(OUTDIR / "Ga2O3-BL-001-4A-TDEP-diagnostic.png", dpi=180)
    plt.close(fig)

    summary = {
        **reference,
        "fit_temperature_K": 300.0,
        "cutoff_angstrom": 5.5,
        "first_order_reference_force_rmse_eV_per_A": (
            float(first_order.group(1)) if first_order else None
        ),
        "harmonic_predicted_force_rms_eV_per_A": values[0],
        "harmonic_force_fit_residual_rms_eV_per_A": values[1],
        "harmonic_force_fit_residual_std_eV_per_A": values[2],
        "harmonic_force_fit_residual_R2": values[3],
        "anharmonicity_sigma_A": values[4],
        "minimum_frequency_THz": float(np.min(frequencies)),
        "maximum_frequency_THz": float(np.max(frequencies)),
        "negative_frequency_values_below_minus_1e-6_THz": int(np.sum(frequencies < -1.0e-6)),
        "total_frequency_values": int(frequencies.size),
    }
    (OUTDIR / "tdep_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
