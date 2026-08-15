#!/usr/bin/env python3
"""Prepare and analyze a selected inclusive MD-frame window for TDEP."""

from __future__ import annotations

import argparse
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


def output_directory(start: int, end: int) -> Path:
    return ROOT / f"TDEP_300K_frames_{start:04d}_{end:04d}"


def prepare(start: int, end: int) -> Path:
    meta = (SOURCE / "infile.meta").read_text().split()
    natoms, nframes = int(meta[0]), int(meta[1])
    timestep_fs = float(meta[2])
    if start < 1 or end < start or end > nframes:
        raise ValueError(
            f"Requested inclusive frames {start}:{end}, but available range is 1:{nframes}"
        )

    positions = np.loadtxt(SOURCE / "infile.positions").reshape(nframes, natoms, 3)
    forces = np.loadtxt(SOURCE / "infile.forces").reshape(nframes, natoms, 3)
    stat_lines = (SOURCE / "infile.stat").read_text().splitlines()
    temperatures = np.asarray([float(line.split()[5]) for line in stat_lines])
    selection = slice(start - 1, end)
    selected_positions = positions[selection]
    selected_forces = forces[selection]
    selected_stat = stat_lines[selection]
    selected_temperatures = temperatures[selection]

    outdir = output_directory(start, end)
    outdir.mkdir(parents=True, exist_ok=True)
    for name in ("infile.ucposcar", "infile.ssposcar", "infile.qpoints_dispersion"):
        shutil.copy2(SOURCE / name, outdir / name)
    np.savetxt(
        outdir / "infile.positions", selected_positions.reshape(-1, 3), fmt="%.16e"
    )
    np.savetxt(
        outdir / "infile.forces", selected_forces.reshape(-1, 3), fmt="%.16e"
    )
    (outdir / "infile.stat").write_text("\n".join(selected_stat) + "\n")
    (outdir / "infile.meta").write_text(
        f"{natoms}\n{len(selected_positions)}\n{timestep_fs:.10f}\n"
        f"{np.mean(selected_temperatures):.10f}\n"
    )
    preparation = {
        "source_total_frames": nframes,
        "frame_start_inclusive": start,
        "frame_end_inclusive": end,
        "n_selected_frames": len(selected_positions),
        "timestep_fs": timestep_fs,
        "selected_time_span_fs": len(selected_positions) * timestep_fs,
        "mean_temperature_K": float(np.mean(selected_temperatures)),
        "temperature_std_K": float(np.std(selected_temperatures, ddof=1)),
        "first_temperature_K": float(selected_temperatures[0]),
        "last_temperature_K": float(selected_temperatures[-1]),
    }
    (outdir / "window_preparation.json").write_text(
        json.dumps(preparation, indent=2) + "\n"
    )
    print(json.dumps(preparation, indent=2))
    return outdir


def analyze(start: int, end: int, cutoff: float) -> None:
    outdir = output_directory(start, end)
    data = np.loadtxt(outdir / "outfile.dispersion_relations")
    x, frequencies = data[:, 0], data[:, 1:]
    log = (outdir / "extract_forceconstants.log").read_text(errors="replace")
    first_order = re.search(r"RMSE total:\s+([-+0-9.Ee]+)", log)
    anharmonicity = re.search(
        r"^\s*second order:\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"<-- anharmonicity measure",
        log,
        re.MULTILINE,
    )
    values = [float(value) for value in anharmonicity.groups()]
    preparation = json.loads((outdir / "window_preparation.json").read_text())
    branch_minima = frequencies.min(axis=0)
    summary = {
        **preparation,
        "cutoff_angstrom": cutoff,
        "first_order_reference_force_rmse_eV_per_A": (
            float(first_order.group(1)) if first_order else None
        ),
        "harmonic_predicted_force_rms_eV_per_A": values[0],
        "harmonic_force_fit_residual_rms_eV_per_A": values[1],
        "harmonic_force_fit_residual_std_eV_per_A": values[2],
        "harmonic_force_fit_residual_R2": values[3],
        "anharmonicity_sigma_A": values[4],
        "minimum_frequency_THz": float(frequencies.min()),
        "maximum_frequency_THz": float(frequencies.max()),
        "negative_branches_below_minus_1e-6_THz": int(
            np.sum(branch_minima < -1.0e-6)
        ),
        "negative_branches_below_minus_0.1_THz": int(
            np.sum(branch_minima < -0.1)
        ),
        "negative_branches_below_minus_1_THz": int(np.sum(branch_minima < -1.0)),
        "negative_frequency_values_below_minus_1e-6_THz": int(
            np.sum(frequencies < -1.0e-6)
        ),
    }
    (outdir / "tdep_window_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    ticks = x[[0, 79, 159, 239, 319]]
    fig, ax = plt.subplots(figsize=(10.5, 6.6), constrained_layout=True)
    for position in ticks:
        ax.axvline(position, color="#aaaaaa", linewidth=0.7, zorder=0)
    for branch in frequencies.T:
        ax.plot(x, branch, color="#245b93", linewidth=1.0)
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xlim(x[0], x[-1])
    ax.set_xticks(ticks, [r"$\Gamma$", "X", "S", "Y", r"$\Gamma$"])
    ax.set_ylabel("Frequency (THz)")
    ax.set_title(
        f"Ga₂O₃ 4 Å bilayer · TDEP frames {start}–{end} "
        f"({len(range(start, end + 1))} frames, {cutoff:g} Å cutoff)"
    )
    ax.grid(axis="y", alpha=0.2)
    figure = outdir / f"Ga2O3-BL-001-4A-TDEP-frames-{start:04d}-{end:04d}.png"
    fig.savefig(figure, dpi=190)
    plt.close(fig)
    print(json.dumps(summary, indent=2))
    print(f"saved {figure}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "analyze"))
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--cutoff", type=float, default=5.5)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.start, args.end)
    else:
        analyze(args.start, args.end, args.cutoff)


if __name__ == "__main__":
    main()
