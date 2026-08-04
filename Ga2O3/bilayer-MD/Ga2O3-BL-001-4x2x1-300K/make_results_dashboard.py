#!/usr/bin/env python3
"""Create a compact dashboard for the current 80-atom MD phonon diagnostics."""

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

from qe_supercell_to_tdep import parse_qe_input, parse_qe_output


ROOT = Path(__file__).resolve().parent
NAVY = "#10243e"
BLUE = "#2b6cb0"
RED = "#c23b2a"
GOLD = "#d99b28"
LIGHT = "#f4f7fb"
GRAY = "#586575"


def add_card(fig, box, title: str, value: str, subtitle: str, color: str) -> None:
    ax = fig.add_axes(box)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color("#d9e1eb")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.06, 0.76, title, color=GRAY, fontsize=10, weight="bold")
    ax.text(0.06, 0.40, value, color=color, fontsize=22, weight="bold")
    ax.text(0.06, 0.12, subtitle, color=GRAY, fontsize=8.5)


def plot_dispersion(ax, path: Path, color: str, title: str) -> tuple[float, float]:
    data = np.loadtxt(path)
    x_values = data[:, 0]
    frequencies = data[:, 1:]
    for band in range(frequencies.shape[1]):
        ax.plot(x_values, frequencies[:, band], color=color, linewidth=0.65, alpha=0.9)
    ax.axhline(0.0, color="#202832", linewidth=0.8)
    ax.set_xlim(x_values[0], x_values[-1])
    ax.set_title(title, loc="left", fontsize=11, color=NAVY, weight="bold")
    ax.set_ylabel("THz", color=GRAY)
    ax.tick_params(labelsize=8, colors=GRAY)
    ax.grid(axis="y", alpha=0.16)
    return float(np.min(frequencies)), float(np.max(frequencies))


def main() -> None:
    input_file = ROOT / "Ga2O3-BL-001-4x2x1.md.in"
    output_file = ROOT / "Ga2O3-BL-001-4x2x1.md.out"
    _cell, symbols, _reference = parse_qe_input(input_file)
    (
        _positions,
        _forces,
        _potential,
        _kinetic,
        temperatures,
        timestep_fs,
        excluded,
    ) = parse_qe_output(output_file, len(symbols))
    temperatures = np.asarray(temperatures, dtype=float)
    time_fs = np.arange(1, len(temperatures) + 1) * timestep_fs

    tdep = json.loads((ROOT / "TDEP_300K/tdep_summary.json").read_text())
    held = json.loads((ROOT / "HELD_300K/held_summary.json").read_text())

    fig = plt.figure(figsize=(16, 9), facecolor=LIGHT)
    header = fig.add_axes([0.0, 0.89, 1.0, 0.11])
    header.set_facecolor(NAVY)
    header.set_xticks([])
    header.set_yticks([])
    for spine in header.spines.values():
        spine.set_visible(False)
    header.text(
        0.035,
        0.58,
        "Ga₂O₃ (001) bilayer · 80-atom MD diagnostics",
        color="white",
        fontsize=22,
        weight="bold",
        va="center",
    )
    header.text(
        0.965,
        0.58,
        "TDEP + experimental HELD",
        color="#cdd9e8",
        fontsize=11,
        ha="right",
        va="center",
    )

    add_card(
        fig,
        [0.035, 0.735, 0.21, 0.12],
        "COMPLETE MD FRAMES",
        f"{len(temperatures)}",
        f"{excluded} unfinished frame excluded",
        BLUE,
    )
    add_card(
        fig,
        [0.265, 0.735, 0.21, 0.12],
        "TRAJECTORY LENGTH",
        f"{time_fs[-1]:.2f} fs",
        f"timestep {timestep_fs:.4f} fs",
        BLUE,
    )
    add_card(
        fig,
        [0.495, 0.735, 0.21, 0.12],
        "MEAN TEMPERATURE",
        f"{temperatures.mean():.1f} K",
        "target = 300 K",
        RED,
    )
    add_card(
        fig,
        [0.725, 0.735, 0.24, 0.12],
        "THERMAL STATUS",
        "HEATING",
        f"latest complete frame = {temperatures[-1]:.1f} K",
        RED,
    )

    ax_temp = fig.add_axes([0.055, 0.425, 0.42, 0.255], facecolor="white")
    ax_temp.plot(time_fs, temperatures, color=RED, linewidth=2.2, marker="o", markersize=3)
    ax_temp.axhline(300.0, color=BLUE, linestyle="--", linewidth=1.4, label="300 K target")
    ax_temp.fill_between(time_fs, 270.0, 330.0, color=BLUE, alpha=0.08)
    ax_temp.set_title("Temperature trajectory", loc="left", fontsize=12, color=NAVY, weight="bold")
    ax_temp.set_xlabel("Time (fs)", color=GRAY)
    ax_temp.set_ylabel("Temperature (K)", color=GRAY)
    ax_temp.tick_params(colors=GRAY)
    ax_temp.grid(alpha=0.16)
    ax_temp.legend(frameon=False, fontsize=9)

    ax_tdep = fig.add_axes([0.535, 0.425, 0.43, 0.255], facecolor="white")
    tdep_min, tdep_max = plot_dispersion(
        ax_tdep,
        ROOT / "TDEP_300K/outfile.dispersion_relations",
        BLUE,
        "TDEP phonon diagnostic · Γ–X–S–Y–Γ",
    )
    ax_tdep.set_xticks([])

    ax_held = fig.add_axes([0.055, 0.065, 0.42, 0.225], facecolor="white")
    held_min, held_max = plot_dispersion(
        ax_held,
        ROOT / "HELD_300K/held_dispersion.dat",
        RED,
        "Experimental HELD diagnostic · Γ–X–S–Y–Γ",
    )
    ax_held.set_xticks([])

    ax_note = fig.add_axes([0.535, 0.065, 0.43, 0.225], facecolor="white")
    ax_note.set_xticks([])
    ax_note.set_yticks([])
    for spine in ax_note.spines.values():
        spine.set_color("#d9e1eb")
    ax_note.text(0.04, 0.86, "Interpretation", fontsize=13, color=NAVY, weight="bold")
    ax_note.text(
        0.04,
        0.66,
        f"TDEP: {tdep_min:.2f} to {tdep_max:.2f} THz  ·  "
        f"residual RMS {tdep['harmonic_force_fit_residual_rms_eV_per_A']:.3f} eV/Å  ·  "
        f"σA {tdep['anharmonicity_sigma_A']:.3f}",
        fontsize=10.5,
        color=BLUE,
        weight="bold",
    )
    ax_note.text(
        0.04,
        0.48,
        f"HELD: {held_min:.2f} to {held_max:.2f} THz  ·  "
        f"force RMSE {held['force_component_rmse_eV_per_A']:.3f} eV/Å",
        fontsize=10.5,
        color=RED,
        weight="bold",
    )
    ax_note.text(
        0.04,
        0.35,
        "Not a converged 300 K phonon result.\n"
        f"The trajectory is only {time_fs[-1]:.0f} fs and averages "
        f"{temperatures.mean():.0f} K.\n"
        "Relax → equilibrate at 300 K → collect a multi-ps production run.",
        fontsize=11,
        color=GRAY,
        linespacing=1.5,
        va="top",
    )

    output = ROOT / "Ga2O3-BL-001-80atom-results-dashboard.png"
    fig.savefig(output, dpi=170, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
