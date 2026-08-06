#!/usr/bin/env python3
"""Create quantitative literature comparisons for the V100/V102 phonons."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
COLORS = {
    "Experiment": "black",
    "Literature DFT": "#2ca02c",
    "V100": "#d62728",
    "V102": "#1f77b4",
}


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")


raman = read_csv(ROOT / "raman_modes.csv")
x = np.arange(len(raman))
labels = [
    r"$A_{1g}^{(1)}$",
    r"$E_g^{(1)}$",
    r"$E_g^{(2)}$",
    r"$E_g^{(3)}$",
    r"$E_g^{(4)}$",
    r"$A_{1g}^{(2)}$",
    r"$E_g^{(5)}$",
]

fig, ax = plt.subplots(figsize=(8.0, 4.8))
series = (
    ("Experiment", "experiment_cm1", "o", "-"),
    ("Literature DFT", "literature_dft_cm1", "s", "--"),
    ("V100", "V100_cm1", "^", "-"),
    ("V102", "V102_cm1", "D", "-"),
)
for name, field, marker, linestyle in series:
    ax.plot(
        x,
        raman[field],
        marker=marker,
        linestyle=linestyle,
        linewidth=1.4,
        markersize=6,
        color=COLORS[name],
        label=name,
    )
ax.set_xticks(x, labels)
ax.set_ylabel(r"Raman frequency (cm$^{-1}$)")
ax.set_title(r"Bulk $\alpha$-Fe$_2$O$_3$: Raman-active $\Gamma$ modes")
ax.grid(axis="y", alpha=0.25)
ax.legend(frameon=False, ncol=2)
fig.tight_layout()
fig.savefig(ROOT / "raman_mode_comparison.pdf")
fig.savefig(ROOT / "raman_mode_comparison.png", dpi=220)
plt.close(fig)


fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0))

# Overall phonon bandwidth.
ax = axes[0, 0]
names = ["Literature\nDFT+U+J", "V100", "V102"]
values = [81.0, 82.8272, 80.6321]
colors = ["#9467bd", COLORS["V100"], COLORS["V102"]]
bars = ax.bar(names, values, color=colors, width=0.65)
ax.bar_label(bars, fmt="%.1f", padding=3)
ax.set_ylim(0, 90)
ax.set_ylabel("Maximum phonon energy (meV)")
ax.set_title("Phonon bandwidth")

# Reported and computed gap around 50 meV.
ax = axes[0, 1]
ax.axvspan(48, 50, color="#9467bd", alpha=0.18, label="Literature: just below 50 meV")
ax.hlines(1, 48.340, 50.824, color=COLORS["V100"], linewidth=9, label="V100")
ax.hlines(0, 47.545, 48.512, color=COLORS["V102"], linewidth=9, label="V102")
ax.scatter([48.340, 50.824], [1, 1], color=COLORS["V100"], s=20)
ax.scatter([47.545, 48.512], [0, 0], color=COLORS["V102"], s=20)
ax.set_xlim(45.5, 52.5)
ax.set_ylim(-0.7, 1.7)
ax.set_yticks([0, 1], ["V102", "V100"])
ax.set_xlabel("Phonon energy (meV)")
ax.set_title("Low-DOS gap near 50 meV")
ax.legend(frameon=False, fontsize=8, loc="upper right")

# Element-resolved mode character.
ax = axes[1, 0]
regions = ["Below 30 meV", "Above 50 meV"]
v100_dominant = [85.30, 91.94]
v102_dominant = [86.25, 92.84]
xx = np.arange(2)
width = 0.34
b1 = ax.bar(xx - width / 2, v100_dominant, width, color=COLORS["V100"], label="V100")
b2 = ax.bar(xx + width / 2, v102_dominant, width, color=COLORS["V102"], label="V102")
ax.bar_label(b1, fmt="%.1f%%", padding=2, fontsize=8)
ax.bar_label(b2, fmt="%.1f%%", padding=2, fontsize=8)
ax.set_xticks(xx, regions)
ax.set_ylabel("Dominant-species contribution (%)")
ax.set_ylim(0, 105)
ax.set_title("Fe character (low E) / O character (high E)")
ax.text(0, 5, "Fe", ha="center", fontweight="bold")
ax.text(1, 5, "O", ha="center", fontweight="bold")
ax.legend(frameon=False)

# Harmonic entropy compared with calorimetry.
ax = axes[1, 1]
names = ["Experiment\n298.15 K", "V100\n300 K", "V102\n300 K"]
values = [87.32, 88.010, 90.406]
errors = [2.0, 0.0, 0.0]
bars = ax.bar(names, values, yerr=errors, capsize=4, color=["#555555", COLORS["V100"], COLORS["V102"]])
ax.bar_label(bars, fmt="%.1f", padding=4)
ax.set_ylim(75, 96)
ax.set_ylabel(r"Entropy (J mol$^{-1}$ K$^{-1}$)")
ax.set_title("Molar entropy per Fe$_2$O$_3$")

fig.suptitle(r"Key literature benchmarks for $\alpha$-Fe$_2$O$_3$ phonons", fontsize=15)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(ROOT / "key_characteristics_comparison.pdf")
fig.savefig(ROOT / "key_characteristics_comparison.png", dpi=220)
plt.close(fig)
