#!/usr/bin/env python3
"""Prepare and plot cumulative TDEP anharmonicity fits."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "TDEP_300K"
RESULTS = SOURCE / "anharmonicity_convergence"


def checkpoints(nframes):
    values = list(range(10, nframes + 1, 10))
    if not values or values[-1] != nframes:
        values.append(nframes)
    return values


def prepare():
    natoms, nframes, timestep_fs, _temperature = (
        SOURCE / "infile.meta"
    ).read_text().split()
    natoms, nframes = int(natoms), int(nframes)
    timestep_fs = float(timestep_fs)
    positions = np.loadtxt(SOURCE / "infile.positions").reshape(nframes, natoms, 3)
    forces = np.loadtxt(SOURCE / "infile.forces").reshape(nframes, natoms, 3)
    stat_lines = (SOURCE / "infile.stat").read_text().splitlines()
    temperatures = np.asarray([float(line.split()[5]) for line in stat_lines])
    RESULTS.mkdir(parents=True, exist_ok=True)
    counts = checkpoints(nframes)
    for count in counts:
        folder = RESULTS / f"nframes_{count:04d}"
        folder.mkdir(parents=True, exist_ok=True)
        for name in ("infile.ucposcar", "infile.ssposcar"):
            shutil.copy2(SOURCE / name, folder / name)
        np.savetxt(folder / "infile.positions", positions[:count].reshape(-1, 3), fmt="%.16e")
        np.savetxt(folder / "infile.forces", forces[:count].reshape(-1, 3), fmt="%.16e")
        (folder / "infile.stat").write_text("\n".join(stat_lines[:count]) + "\n")
        (folder / "infile.meta").write_text(
            f"{natoms}\n{count}\n{timestep_fs:.10f}\n{np.mean(temperatures[:count]):.10f}\n"
        )
    (RESULTS / "checkpoints.json").write_text(json.dumps(counts) + "\n")
    print(f"Prepared cumulative fits at {counts}")


def parse_force_line(text):
    match = re.search(
        r"^\s*second order:\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"<-- anharmonicity measure",
        text,
        re.MULTILINE,
    )
    if not match:
        raise ValueError("Could not find TDEP anharmonicity line")
    values = [float(value) for value in match.groups()]
    return {
        "harmonic_force_rms_eV_per_A": values[0],
        "residual_force_rms_eV_per_A": values[1],
        "residual_force_std_eV_per_A": values[2],
        "residual_R2": values[3],
        "sigma_A": values[4],
    }


def analyze():
    counts = json.loads((RESULTS / "checkpoints.json").read_text())
    rows = []
    for count in counts:
        folder = RESULTS / f"nframes_{count:04d}"
        meta = (folder / "infile.meta").read_text().split()
        rows.append({
            "frames": count,
            "time_fs": count * float(meta[2]),
            "mean_temperature_K": float(meta[3]),
            **parse_force_line((folder / "extract_forceconstants.log").read_text()),
        })
    with (RESULTS / "tdep_anharmonicity_convergence.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "definition": "sigma_A = std(F_DFT - F_harmonic) / std(F_DFT)",
        "final": rows[-1],
        "checkpoints": rows,
        "interpretation": (
            "Diagnostic only: the cumulative sample includes the initial cooling "
            "transient and is not an equilibrated 300 K ensemble."
        ),
    }
    (RESULTS / "tdep_anharmonicity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    time_fs = np.asarray([row["time_fs"] for row in rows])
    sigma = np.asarray([row["sigma_A"] for row in rows])
    residual_std = np.asarray([row["residual_force_std_eV_per_A"] for row in rows])
    temperatures = np.asarray([row["mean_temperature_K"] for row in rows])
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), constrained_layout=True)
    fig.suptitle(
        "Ga₂O₃ 4 Å bilayer · cumulative TDEP anharmonicity",
        fontsize=16, weight="bold",
    )
    axes[0].plot(time_fs, sigma, marker="o", color="#c23b2a", linewidth=2)
    axes[0].set_ylabel(r"Anharmonicity $\sigma_A$")
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Normalized force residual")
    axes[1].plot(time_fs, residual_std, marker="o", color="#2b6cb0", linewidth=2)
    axes[1].set_ylabel("Residual-force std. (eV/Å)")
    axes[1].set_title("Unexplained force fluctuations")
    axes[2].plot(time_fs, temperatures, marker="o", color="#d99b28", linewidth=2)
    axes[2].axhline(300, color="#2b6cb0", linestyle="--", linewidth=1.2)
    axes[2].set_ylabel("Cumulative mean temperature (K)")
    axes[2].set_title("Sampled temperature")
    for axis in axes:
        axis.set_xlabel("Trajectory included (fs)")
        axis.grid(alpha=0.2)
    axes[0].annotate(
        f"final = {rows[-1]['sigma_A']:.3f}",
        xy=(time_fs[-1], sigma[-1]), xytext=(-65, -28),
        textcoords="offset points", color="#c23b2a", weight="bold",
        arrowprops={"arrowstyle": "->", "color": "#c23b2a"},
    )
    fig.savefig(
        RESULTS / "Ga2O3-BL-001-4A-TDEP-anharmonicity-convergence.png", dpi=190
    )
    plt.close(fig)
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "analyze"))
    args = parser.parse_args()
    prepare() if args.mode == "prepare" else analyze()


if __name__ == "__main__":
    main()
