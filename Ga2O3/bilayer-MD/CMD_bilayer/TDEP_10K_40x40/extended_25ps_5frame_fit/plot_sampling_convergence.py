#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load(path):
    rows = []
    for line in path.read_text().splitlines():
        try:
            rows.append([float(x) for x in line.split()])
        except ValueError:
            pass
    data = np.asarray(rows)
    return np.linspace(0, 4, len(data)), data[:, 1:]


datasets = [
    ("0 K finite displacement", ROOT / "phonopy_0K_tabgap_40x40_amp002" / "dispersion_0K_tabgap.dat", "#2563eb"),
    ("10 K TDEP: original 10 ps window", ROOT / "tdep_10K_40x40" / "outfile.dispersion_relations", "#f59e0b"),
    ("10 K TDEP: independent 25 ps window", HERE / "outfile.dispersion_relations", "#dc2626"),
]

loaded = [(name, *load(path), color) for name, path, color in datasets]
fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                         gridspec_kw={"height_ratios": [2.1, 1]})
for ax in axes:
    ax.axhline(0, color="black", lw=0.9)
    for tick in range(5):
        ax.axvline(tick, color="0.85", lw=0.8, zorder=0)
    for name, x, y, color in loaded:
        for branch in range(y.shape[1]):
            ax.plot(x, y[:, branch], color=color, lw=0.5, alpha=0.58)
    ax.grid(axis="y", color="0.93", lw=0.6)

for name, x, y, color in loaded:
    axes[0].plot([], [], color=color, lw=1.7,
                 label=f"{name} (min {y.min():.3f} THz)")
axes[0].legend(loc="upper right", fontsize=8.5)
axes[0].set_ylabel("Frequency (THz)")
axes[0].set_title("16,000-atom Ga$_2$O$_3$ bilayer: 10 K TDEP sampling convergence\n"
                  "Second-order cutoff 5.5 Angstrom")

axes[1].set_ylim(-0.55, 3.5)
axes[1].fill_between([0, 4], -0.55, 0, color="#fee2e2", alpha=0.45)
axes[1].set_title("Low-frequency zoom")
axes[1].set_ylabel("Frequency (THz)")
axes[1].set_xlabel("Wave-vector path")
axes[1].set_xticks(range(5), [r"$\Gamma$", "X", "S", "Y", r"$\Gamma$"])

fig.tight_layout()
out = ROOT / "phonon_sampling_convergence_10K_40x40.png"
fig.savefig(out, dpi=240, bbox_inches="tight")
print(out)
