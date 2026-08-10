#!/usr/bin/env python3
"""Overlay V100, V102, and V104 phonon bands and total DOS."""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
CASES = {
    "V100": (ROOT / "phonopy_V100_2x2x2", "#d62728"),
    "V102": (ROOT / "phonopy_V102_2x2x2", "#1f77b4"),
    "V104": (ROOT / "phonopy_V104_2x2x2", "#ff7f0e"),
}
GROUPS = ((0, 1), (2, 3), (4, 5, 6, 7), (8,))


def load_case(directory):
    with h5py.File(directory / "band.hdf5") as band:
        result = {
            "distance": band["distance"][:],
            "frequency": band["frequency"][:],
            "label": band["label"][:],
        }
    result["dos"] = np.loadtxt(directory / "total_dos.dat")
    return result


def format_label(label):
    label = label.decode() if isinstance(label, bytes) else str(label)
    return f"${label}$" if "_" in label and "$" not in label else label


data = {name: load_case(directory) for name, (directory, _) in CASES.items()}
reference = data["V100"]
widths = [
    reference["distance"][group[-1], -1] - reference["distance"][group[0], 0]
    for group in GROUPS
]

fig, axes = plt.subplots(
    1,
    len(GROUPS) + 1,
    figsize=(7.6, 5.1),
    sharey=True,
    gridspec_kw={"width_ratios": [*widths, max(widths) * 0.78], "wspace": 0.13},
)

for ax, group in zip(axes[:-1], GROUPS):
    ref_left = reference["distance"][group[0], 0]
    ticks = [ref_left]
    ticklabels = [format_label(reference["label"][group[0], 0])]
    for segment in group:
        ref_x = reference["distance"][segment]
        ticks.append(ref_x[-1])
        ticklabels.append(format_label(reference["label"][segment, 1]))
        for name, (_, color) in CASES.items():
            source_x = data[name]["distance"][segment]
            aligned_x = np.interp(source_x, (source_x[0], source_x[-1]), (ref_x[0], ref_x[-1]))
            ax.plot(aligned_x, data[name]["frequency"][segment], color=color, linewidth=0.9, alpha=0.82)
    ax.axhline(0, color="black", linestyle=":", linewidth=0.7)
    ax.set_xlim(ref_left, reference["distance"][group[-1], -1])
    ax.set_xticks(ticks, ticklabels)
    ax.tick_params(direction="in", top=True, right=True)

axes[0].set_ylabel("Frequency (THz)")
for name, (_, color) in CASES.items():
    dos = data[name]["dos"]
    axes[-1].plot(dos[:, 1], dos[:, 0], color=color, linewidth=1.35, label=name)
axes[-1].axhline(0, color="black", linestyle=":", linewidth=0.7)
axes[-1].set_xlabel("DOS")
axes[-1].set_xlim(left=0)
axes[-1].tick_params(direction="in", top=True, right=True)
axes[-1].legend(frameon=False, fontsize=9, loc="upper right")

all_frequencies = np.concatenate([case["frequency"].ravel() for case in data.values()])
all_dos_frequencies = np.concatenate([case["dos"][:, 0] for case in data.values()])
axes[0].set_ylim(
    min(-0.5, float(all_frequencies.min()) - 0.5),
    max(float(all_frequencies.max()), float(all_dos_frequencies.max())) + 0.5,
)
fig.suptitle(r"Fe$_2$O$_3$ phonons: V100, V102, and V104 (2$\times$2$\times$2)")
fig.subplots_adjust(top=0.85, bottom=0.16, left=0.10, right=0.98)
fig.savefig(ROOT / "Fe2O3_V100_V102_V104_phonon_comparison.pdf")
fig.savefig(ROOT / "Fe2O3_V100_V102_V104_phonon_comparison.png", dpi=220)
