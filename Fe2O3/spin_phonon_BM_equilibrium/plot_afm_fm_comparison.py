#!/usr/bin/env python3
"""Plot matched AFM/FM phonons and tabulate Gamma-point frequency shifts."""

from pathlib import Path
import csv

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parent
CASES = {
    "AFM (++--)": (ROOT / "AFM", "#1f77b4"),
    "FM (++++)": (ROOT / "FM", "#d62728"),
}
GROUPS = ((0, 1), (2, 3), (4, 5, 6, 7), (8,))


def load_case(directory):
    with h5py.File(directory / "band.hdf5") as band:
        result = {
            "distance": band["distance"][:],
            "frequency": band["frequency"][:],
            "label": band["label"][:],
            "eigenvector": band["eigenvector"][:],
        }
    result["dos"] = np.loadtxt(directory / "total_dos.dat")
    return result


def format_label(label):
    label = label.decode() if isinstance(label, bytes) else str(label)
    return f"${label}$" if "_" in label and "$" not in label else label


data = {name: load_case(directory) for name, (directory, _) in CASES.items()}
reference = data["AFM (++--)"]
widths = [
    reference["distance"][group[-1], -1] - reference["distance"][group[0], 0]
    for group in GROUPS
]
fig, axes = plt.subplots(
    1,
    len(GROUPS) + 1,
    figsize=(7.8, 5.1),
    sharey=True,
    gridspec_kw={"width_ratios": [*widths, max(widths) * 0.78], "wspace": 0.13},
)

for ax, group in zip(axes[:-1], GROUPS):
    left = reference["distance"][group[0], 0]
    ticks = [left]
    ticklabels = [format_label(reference["label"][group[0], 0])]
    for segment in group:
        ref_x = reference["distance"][segment]
        ticks.append(ref_x[-1])
        ticklabels.append(format_label(reference["label"][segment, 1]))
        for name, (_, color) in CASES.items():
            source_x = data[name]["distance"][segment]
            aligned_x = np.interp(
                source_x, (source_x[0], source_x[-1]), (ref_x[0], ref_x[-1])
            )
            ax.plot(
                aligned_x,
                data[name]["frequency"][segment],
                color=color,
                linewidth=0.95,
                alpha=0.82,
            )
    ax.axhline(0, color="black", linestyle=":", linewidth=0.7)
    ax.set_xlim(left, reference["distance"][group[-1], -1])
    ax.set_xticks(ticks, ticklabels)
    ax.tick_params(direction="in", top=True, right=True)

axes[0].set_ylabel("Frequency (THz)")
for name, (_, color) in CASES.items():
    dos = data[name]["dos"]
    axes[-1].plot(dos[:, 1], dos[:, 0], color=color, linewidth=1.4, label=name)
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
fig.suptitle(
    r"Fe$_2$O$_3$ spin--phonon comparison at BM equilibrium "
    r"(2$\times$2$\times$2)"
)
fig.subplots_adjust(top=0.85, bottom=0.16, left=0.10, right=0.98)
fig.savefig(ROOT / "Fe2O3_BMeq_AFM_FM_phonon_comparison.pdf")
fig.savefig(ROOT / "Fe2O3_BMeq_AFM_FM_phonon_comparison.png", dpi=220)

afm_gamma = data["AFM (++--)"]["frequency"][0, 0]
fm_gamma = data["FM (++++)"]["frequency"][0, 0]
afm_vectors = data["AFM (++--)"]["eigenvector"][0, 0]
fm_vectors = data["FM (++++)"]["eigenvector"][0, 0]
overlap = np.abs(afm_vectors.conj() @ fm_vectors.T) ** 2
afm_indices, fm_indices = linear_sum_assignment(-overlap)
matches = sorted(zip(afm_indices, fm_indices), key=lambda pair: afm_gamma[pair[0]])
with (ROOT / "Fe2O3_BMeq_AFM_FM_Gamma_shifts.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(
        ["mode_sorted", "AFM_THz", "FM_THz", "delta_FM_minus_AFM_THz", "delta_cm-1", "eigenvector_overlap"]
    )
    for index, (afm_mode, fm_mode) in enumerate(matches, 1):
        afm_frequency = afm_gamma[afm_mode]
        fm_frequency = fm_gamma[fm_mode]
        delta = fm_frequency - afm_frequency
        writer.writerow(
            [index, f"{afm_frequency:.10f}", f"{fm_frequency:.10f}", f"{delta:.10f}", f"{delta * 33.35640952:.8f}", f"{overlap[afm_mode, fm_mode]:.8f}"]
        )

print("Wrote AFM-FM comparison PDF/PNG and eigenvector-matched Gamma shifts CSV.")
