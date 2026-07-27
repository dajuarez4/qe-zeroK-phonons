#!/usr/bin/env python3
"""Plot the original raw Al2O3 dispersion without clipping negative values."""

from pathlib import Path
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
SEGMENTS = [
    ("Al2O3.mp_seg1.freq.gp", [0, 30, 60], [r"$\Gamma$", "L", r"B$_1$"]),
    ("Al2O3.mp_seg2.freq.gp", [0, 30, 60, 90], ["B", "Z", r"$\Gamma$", "X"]),
    ("Al2O3.mp_seg3.freq.gp", [0, 30, 60, 90], ["Q", "F", r"P$_1$", "Z"]),
    ("Al2O3.mp_seg4.freq.gp", [0, 30], ["L", "P"]),
]


def main():
    datasets = []
    negative_rows = []
    for segment_index, (filename, ticks, labels) in enumerate(SEGMENTS, start=1):
        data = np.loadtxt(ROOT / filename)
        frequencies = data[:, 1:]  # Raw values: deliberately no clipping.
        datasets.append((data, frequencies, ticks, labels))
        for row_index, mode_index in np.argwhere(frequencies < 0.0):
            negative_rows.append(
                {
                    "segment": segment_index,
                    "path_file": filename,
                    "row_index_zero_based": int(row_index),
                    "path_coordinate": float(data[row_index, 0]),
                    "mode": int(mode_index + 1),
                    "frequency_cm-1": float(frequencies[row_index, mode_index]),
                }
            )

    output_csv = ROOT / "Al2O3_raw_negative_frequency_samples.csv"
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=negative_rows[0].keys())
        writer.writeheader()
        writer.writerows(negative_rows)

    widths = [2, 3, 3, 1]
    fig = plt.figure(figsize=(13.2, 8.4), facecolor="white")
    outer = fig.add_gridspec(
        2, 4, height_ratios=[1.25, 1], width_ratios=widths,
        left=0.075, right=0.985, bottom=0.095, top=0.885,
        hspace=0.19, wspace=0.045,
    )

    overview_axes = []
    zoom_axes = []
    for column, (data, frequencies, ticks, labels) in enumerate(datasets):
        overview = fig.add_subplot(outer[0, column], sharey=overview_axes[0] if overview_axes else None)
        zoom = fig.add_subplot(outer[1, column], sharey=zoom_axes[0] if zoom_axes else None)
        overview_axes.append(overview)
        zoom_axes.append(zoom)

        x = np.arange(len(data))
        for mode in range(frequencies.shape[1]):
            branch = frequencies[:, mode]
            overview.plot(x, branch, color="#1559a6", lw=0.82, alpha=0.92)
            zoom.plot(x, branch, color="#1559a6", lw=1.0, alpha=0.95)

        negative_indices = np.argwhere(frequencies < 0.0)
        if len(negative_indices):
            overview.scatter(
                negative_indices[:, 0], frequencies[negative_indices[:, 0], negative_indices[:, 1]],
                color="#d62728", s=18, zorder=5,
            )
            zoom.scatter(
                negative_indices[:, 0], frequencies[negative_indices[:, 0], negative_indices[:, 1]],
                color="#d62728", edgecolor="white", linewidth=0.35, s=35, zorder=6,
            )

        for axis in (overview, zoom):
            axis.axhline(0.0, color="#111111", lw=0.9, ls="--", zorder=1)
            for tick in ticks:
                axis.axvline(tick, color="0.78", lw=0.65, zorder=0)
            axis.set_xlim(0, len(data) - 1)
            axis.set_xticks(ticks, labels)
            axis.tick_params(axis="x", labelsize=10)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)

        overview.set_ylim(-15.0, 880.0)
        zoom.set_ylim(-12.0, 80.0)
        overview.set_xticklabels([])
        if column > 0:
            overview.tick_params(labelleft=False)
            zoom.tick_params(labelleft=False)

    overview_axes[0].set_ylabel(r"Raw frequency (cm$^{-1}$)", fontsize=12)
    zoom_axes[0].set_ylabel(r"Raw frequency (cm$^{-1}$)", fontsize=12)
    overview_axes[0].text(
        0.01, 0.96, "Full dispersion — raw QE values",
        transform=overview_axes[0].transAxes, va="top", fontsize=11, fontweight="bold",
    )
    zoom_axes[0].text(
        0.01, 0.94, "Zoom around zero",
        transform=zoom_axes[0].transAxes, va="top", fontsize=11, fontweight="bold",
    )

    minimum = min(row["frequency_cm-1"] for row in negative_rows)
    fig.suptitle(
        r"Original $alpha$-Al$_2$O$_3$ phonons: unclipped acoustic frequencies near $Gamma$",
        fontsize=17, y=0.965,
    )
    fig.text(
        0.075, 0.915,
        f"Raw minimum = {minimum:.4f} cm⁻¹  •  {len(negative_rows)} negative samples  •  "
        "red markers are values that earlier presentation plots set to zero",
        fontsize=10.5, color="#4a4a4a",
    )
    fig.text(0.53, 0.025, "Wave-vector path", ha="center", fontsize=12)

    png = ROOT / "Al2O3_raw_imaginary_phonons_zoom.png"
    pdf = ROOT / "Al2O3_raw_imaginary_phonons_zoom.pdf"
    fig.savefig(png, dpi=230)
    fig.savefig(pdf)
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")
    print(f"Wrote {output_csv}")
    print(f"Minimum raw frequency: {minimum:.4f} cm^-1")


if __name__ == "__main__":
    main()

