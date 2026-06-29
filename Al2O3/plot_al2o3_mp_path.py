from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SEGMENTS = [
    {
        "file": Path("Al2O3.mp_seg1.freq.gp"),
        "labels": [r"$\Gamma$", "L", r"$B_1$"],
        "tick_indices": [0, 30, 60],
    },
    {
        "file": Path("Al2O3.mp_seg2.freq.gp"),
        "labels": ["B", "Z", r"$\Gamma$", "X"],
        "tick_indices": [0, 30, 60, 90],
    },
    {
        "file": Path("Al2O3.mp_seg3.freq.gp"),
        "labels": ["Q", "F", r"$P_1$", "Z"],
        "tick_indices": [0, 30, 60, 90],
    },
    {
        "file": Path("Al2O3.mp_seg4.freq.gp"),
        "labels": ["L", "P"],
        "tick_indices": [0, 30],
    },
]

GAP = 0.08


def load_segment(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected format in {path}")
    x = data[:, 0].copy()
    freqs = data[:, 1:].copy()
    # Small negative values near Gamma are clipped only for plotting.
    freqs[freqs < 0.0] = 0.0
    return x, freqs


plt.rcParams.update(
    {
        "font.size": 14,
        "axes.labelsize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "axes.linewidth": 1.2,
    }
)

fig, ax = plt.subplots(figsize=(9, 7))

offset = 0.0
all_tick_positions: list[float] = []
all_tick_labels: list[str] = []

for i, seg in enumerate(SEGMENTS):
    x_local, freqs = load_segment(seg["file"])
    x_plot = x_local + offset

    for mode in range(freqs.shape[1]):
        ax.plot(x_plot, freqs[:, mode], color="blue", linewidth=1.2)

    for j, idx in enumerate(seg["tick_indices"]):
        if idx >= len(x_plot):
            continue
        pos = x_plot[idx]
        label = seg["labels"][j]
        if i > 0 and j == 0:
            # Start of a disconnected segment: keep the label, but do not
            # duplicate the vertical line if it lands on the previous end.
            all_tick_positions.append(pos)
            all_tick_labels.append(label)
            ax.axvline(pos, color="black", linewidth=1.0)
            continue
        all_tick_positions.append(pos)
        all_tick_labels.append(label)
        ax.axvline(pos, color="black", linewidth=1.0)

    offset = x_plot[-1] + GAP

ax.set_xticks(all_tick_positions)
ax.set_xticklabels(all_tick_labels)
ax.set_ylabel(r"Frequencies (cm$^{-1}$)")
ax.set_xlabel("Wave Vector")
ax.set_ylim(bottom=0)
ax.set_xlim(all_tick_positions[0], all_tick_positions[-1])
ax.tick_params(direction="in", top=True, right=True, length=6, width=1.1)
ax.spines["top"].set_visible(True)
ax.spines["right"].set_visible(True)

plt.tight_layout()
plt.savefig("Al2O3_mp_style_phonon_dispersion.png", dpi=300, bbox_inches="tight")
plt.show()
