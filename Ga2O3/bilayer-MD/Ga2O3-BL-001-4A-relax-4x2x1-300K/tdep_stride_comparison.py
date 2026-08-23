#!/usr/bin/env python3
"""Prepare and summarize decorrelated TDEP stride/cutoff comparisons."""

from __future__ import annotations

import argparse
import csv
import json
import os
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
RESULTS = ROOT / "TDEP_stride_comparison"
STRIDES = (10, 20)
CUTOFFS = (5.0, 5.5, 6.0)


def cutoff_label(cutoff: float) -> str:
    return f"{cutoff:g}".replace(".", "p")


def output_directory(stride: int, cutoff: float) -> Path:
    return ROOT / f"TDEP_300K_stride_{stride:02d}_cutoff_{cutoff_label(cutoff)}A"


def selected_indices(nframes: int, stride: int) -> np.ndarray:
    """Return zero-based, nested samples ending at the newest frame."""
    return np.asarray(list(range(nframes - 1, -1, -stride))[::-1], dtype=int)


def prepare() -> None:
    meta = (SOURCE / "infile.meta").read_text().split()
    natoms, nframes = int(meta[0]), int(meta[1])
    timestep_fs = float(meta[2])
    positions = np.loadtxt(SOURCE / "infile.positions").reshape(nframes, natoms, 3)
    forces = np.loadtxt(SOURCE / "infile.forces").reshape(nframes, natoms, 3)
    stat_lines = (SOURCE / "infile.stat").read_text().splitlines()
    temperatures = np.asarray([float(line.split()[5]) for line in stat_lines])
    source_reference = json.loads((SOURCE / "reference_summary.json").read_text())

    prepared = []
    for stride in STRIDES:
        indices = selected_indices(nframes, stride)
        for cutoff in CUTOFFS:
            outdir = output_directory(stride, cutoff)
            outdir.mkdir(parents=True, exist_ok=True)
            for name in (
                "infile.ucposcar",
                "infile.ssposcar",
                "infile.qpoints_dispersion",
            ):
                shutil.copy2(SOURCE / name, outdir / name)
            np.savetxt(
                outdir / "infile.positions",
                positions[indices].reshape(-1, 3),
                fmt="%.16e",
            )
            np.savetxt(
                outdir / "infile.forces",
                forces[indices].reshape(-1, 3),
                fmt="%.16e",
            )
            (outdir / "infile.stat").write_text(
                "\n".join(stat_lines[index] for index in indices) + "\n"
            )
            mean_temperature = float(np.mean(temperatures[indices]))
            (outdir / "infile.meta").write_text(
                f"{natoms}\n{len(indices)}\n{stride * timestep_fs:.10f}\n"
                f"{mean_temperature:.10f}\n"
            )
            reference = {
                **source_reference,
                "n_frames": int(len(indices)),
                "source_total_frames": nframes,
                "selection_stride_steps": stride,
                "source_timestep_fs": timestep_fs,
                "effective_sample_spacing_fs": stride * timestep_fs,
                "first_selected_frame_1based": int(indices[0] + 1),
                "last_selected_frame_1based": int(indices[-1] + 1),
                "selected_trajectory_span_fs": float(
                    (indices[-1] - indices[0]) * timestep_fs
                ),
                "mean_temperature_K": mean_temperature,
                "latest_complete_temperature_K": float(temperatures[indices[-1]]),
                "cutoff_angstrom": cutoff,
            }
            (outdir / "reference_summary.json").write_text(
                json.dumps(reference, indent=2) + "\n"
            )
            prepared.append(
                {
                    "directory": outdir.name,
                    "stride": stride,
                    "cutoff_angstrom": cutoff,
                    "n_frames": int(len(indices)),
                    "first_frame": int(indices[0] + 1),
                    "last_frame": int(indices[-1] + 1),
                }
            )
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "prepared_fits.json").write_text(
        json.dumps(prepared, indent=2) + "\n"
    )
    print(json.dumps(prepared, indent=2))


def analyze() -> None:
    rows = []
    for stride in STRIDES:
        for cutoff in CUTOFFS:
            outdir = output_directory(stride, cutoff)
            summary = json.loads((outdir / "tdep_summary.json").read_text())
            rows.append(
                {
                    "sample": f"stride_{stride}",
                    "stride_steps": stride,
                    "sample_spacing_fs": summary["effective_sample_spacing_fs"],
                    "n_frames": summary["n_frames"],
                    "cutoff_angstrom": cutoff,
                    "mean_temperature_K": summary["mean_temperature_K"],
                    "residual_rms_eV_per_A": summary[
                        "harmonic_force_fit_residual_rms_eV_per_A"
                    ],
                    "residual_R2": summary["harmonic_force_fit_residual_R2"],
                    "sigma_A": summary["anharmonicity_sigma_A"],
                    "minimum_frequency_THz": summary["minimum_frequency_THz"],
                    "maximum_frequency_THz": summary["maximum_frequency_THz"],
                    "negative_frequency_values": summary[
                        "negative_frequency_values_below_minus_1e-6_THz"
                    ],
                    "directory": outdir.name,
                }
            )

    baseline = json.loads((SOURCE / "tdep_summary.json").read_text())
    baseline_row = {
        "sample": "all_frames",
        "stride_steps": 1,
        "sample_spacing_fs": baseline["timestep_fs"],
        "n_frames": baseline["n_frames"],
        "cutoff_angstrom": baseline["cutoff_angstrom"],
        "mean_temperature_K": baseline["mean_temperature_K"],
        "residual_rms_eV_per_A": baseline[
            "harmonic_force_fit_residual_rms_eV_per_A"
        ],
        "residual_R2": baseline["harmonic_force_fit_residual_R2"],
        "sigma_A": baseline["anharmonicity_sigma_A"],
        "minimum_frequency_THz": baseline["minimum_frequency_THz"],
        "maximum_frequency_THz": baseline["maximum_frequency_THz"],
        "negative_frequency_values": baseline[
            "negative_frequency_values_below_minus_1e-6_THz"
        ],
        "directory": SOURCE.name,
    }
    all_rows = [baseline_row, *rows]
    RESULTS.mkdir(parents=True, exist_ok=True)
    with (RESULTS / "tdep_stride_cutoff_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    (RESULTS / "tdep_stride_cutoff_comparison.json").write_text(
        json.dumps({"baseline": baseline_row, "decorrelated_fits": rows}, indent=2)
        + "\n"
    )

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.0), sharex=True, sharey=True)
    finite_min, finite_max = 0.0, 0.0
    datasets = {}
    for row, stride in enumerate(STRIDES):
        for column, cutoff in enumerate(CUTOFFS):
            outdir = output_directory(stride, cutoff)
            data = np.loadtxt(outdir / "outfile.dispersion_relations")
            x, frequencies = data[:, 0], data[:, 1:]
            datasets[(stride, cutoff)] = (x, frequencies)
            finite_min = min(finite_min, float(frequencies.min()))
            finite_max = max(finite_max, float(frequencies.max()))
            axis = axes[row, column]
            ticks = x[[0, 79, 159, 239, 319]]
            for position in ticks:
                axis.axvline(position, color="#aaaaaa", linewidth=0.6)
            for branch in frequencies.T:
                axis.plot(x, branch, color="#245b93", linewidth=0.75)
            axis.axhline(0, color="black", linewidth=0.8)
            axis.set_title(
                f"stride {stride} · {cutoff:g} Å\n"
                f"{len(selected_indices(baseline['n_frames'], stride))} frames"
            )
            axis.set_xticks(ticks, ["Γ", "X", "S", "Y", "Γ"])
            axis.grid(axis="y", alpha=0.15)
    ypad = 0.04 * (finite_max - finite_min)
    for axis in axes.flat:
        axis.set_ylim(finite_min - ypad, finite_max + ypad)
    for axis in axes[:, 0]:
        axis.set_ylabel("Frequency (THz)")
    fig.suptitle(
        "Ga₂O₃ 4 Å bilayer · decorrelated TDEP stride/cutoff comparison",
        fontsize=16,
        weight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(RESULTS / "tdep_stride_cutoff_dispersion_comparison.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4), constrained_layout=True)
    for stride, color in zip(STRIDES, ("#245b93", "#c23b2a")):
        selected = [item for item in rows if item["stride_steps"] == stride]
        cutoffs = [item["cutoff_angstrom"] for item in selected]
        axes[0].plot(
            cutoffs,
            [item["minimum_frequency_THz"] for item in selected],
            marker="o",
            color=color,
            label=f"stride {stride}",
        )
        axes[1].plot(
            cutoffs,
            [item["residual_rms_eV_per_A"] for item in selected],
            marker="o",
            color=color,
        )
        axes[2].plot(
            cutoffs,
            [item["sigma_A"] for item in selected],
            marker="o",
            color=color,
        )
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_ylabel("Minimum frequency (THz)")
    axes[1].set_ylabel("Force residual RMS (eV/Å)")
    axes[2].set_ylabel(r"Anharmonicity $\sigma_A$")
    for axis in axes:
        axis.set_xlabel("Second-order cutoff (Å)")
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False)
    fig.suptitle("Decorrelated TDEP sensitivity metrics", fontsize=15, weight="bold")
    fig.savefig(RESULTS / "tdep_stride_cutoff_metrics.png", dpi=190)
    plt.close(fig)
    print(json.dumps({"baseline": baseline_row, "decorrelated_fits": rows}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "analyze"))
    args = parser.parse_args()
    prepare() if args.mode == "prepare" else analyze()


if __name__ == "__main__":
    main()
