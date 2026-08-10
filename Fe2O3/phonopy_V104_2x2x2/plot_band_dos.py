#!/usr/bin/env python3
"""Plot the V104 phonon dispersion and total DOS."""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
GROUPS = ((0, 1), (2, 3), (4, 5, 6, 7), (8,))


def format_label(label):
    label = label.decode() if isinstance(label, bytes) else str(label)
    return f"${label}$" if "_" in label and "$" not in label else label


with h5py.File(ROOT / "band.hdf5") as band:
    distances = band["distance"][:]
    frequencies = band["frequency"][:]
    labels = band["label"][:]

dos = np.loadtxt(ROOT / "total_dos.dat")
widths = [distances[g[-1], -1] - distances[g[0], 0] for g in GROUPS]
fig, axes = plt.subplots(
    1,
    len(GROUPS) + 1,
    figsize=(6.4, 4.8),
    sharey=True,
    gridspec_kw={"width_ratios": [*widths, max(widths) * 0.68], "wspace": 0.12},
)

for ax, group in zip(axes[:-1], GROUPS):
    left = distances[group[0], 0]
    ticks = [left]
    ticklabels = [format_label(labels[group[0], 0])]
    for segment in group:
        x = distances[segment]
        ax.plot(x, frequencies[segment], color="red", linewidth=1.05)
        ticks.append(x[-1])
        ticklabels.append(format_label(labels[segment, 1]))
    ax.axhline(0, color="blue", linestyle=":", linewidth=0.7)
    ax.set_xlim(left, distances[group[-1], -1])
    ax.set_xticks(ticks, ticklabels)
    ax.tick_params(direction="in", top=True, right=True)

axes[0].set_ylabel("Frequency (THz)")
axes[-1].plot(dos[:, 1], dos[:, 0], linewidth=1.4)
axes[-1].axhline(0, color="blue", linestyle=":", linewidth=0.7)
axes[-1].set_xlabel("DOS")
axes[-1].set_xlim(left=0)
axes[-1].tick_params(direction="in", top=True, right=True)

ymin = min(-0.5, float(np.nanmin(frequencies)) - 0.5)
ymax = max(float(np.nanmax(frequencies)), float(np.nanmax(dos[:, 0]))) + 0.5
axes[0].set_ylim(ymin, ymax)
fig.suptitle(r"Fe$_2$O$_3$ V104 phonons (2$\times$2$\times$2)")
fig.subplots_adjust(top=0.84, bottom=0.16, left=0.12, right=0.98)
fig.savefig(ROOT / "Fe2O3_V104_phonon_band_dos.pdf")
