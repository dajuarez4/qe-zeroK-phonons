from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ORIGINAL_FILE = Path("Al2O3.freq.gp")
ORIGINAL_LABELS = [r"$\Gamma$", "X", "K", r"$\Gamma$", "Z", "U", "H", "Z"]
ORIGINAL_TICK_INDICES = [0, 30, 60, 90, 120, 150, 180, 210]

MP_SEGMENTS = [
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

MP_GAP = 0.08


def load_gp(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected format in {path}")
    x = data[:, 0].copy()
    freqs = data[:, 1:].copy()
    freqs[freqs < 0.0] = 0.0
    return x, freqs


def plot_original(ax: plt.Axes) -> None:
    x, freqs = load_gp(ORIGINAL_FILE)

    for mode in range(freqs.shape[1]):
        ax.plot(x, freqs[:, mode], color="black", linewidth=1.2)

    tick_positions = [x[i] for i in ORIGINAL_TICK_INDICES if i < len(x)]
    tick_labels = ORIGINAL_LABELS[: len(tick_positions)]

    for xpos in tick_positions:
        ax.axvline(x=xpos, color="gray", linewidth=1.0)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(bottom=0)
    ax.set_title("QE Path")
    ax.set_xlabel("Wave Vector")
    ax.tick_params(direction="in", top=True, right=True, length=6, width=1.1)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


def plot_mp_style(ax: plt.Axes) -> None:
    offset = 0.0
    all_tick_positions: list[float] = []
    all_tick_labels: list[str] = []

    for seg in MP_SEGMENTS:
        x_local, freqs = load_gp(seg["file"])
        x_plot = x_local + offset

        for mode in range(freqs.shape[1]):
            ax.plot(x_plot, freqs[:, mode], color="blue", linewidth=1.2)

        for idx, label in zip(seg["tick_indices"], seg["labels"]):
            if idx >= len(x_plot):
                continue
            pos = x_plot[idx]
            all_tick_positions.append(pos)
            all_tick_labels.append(label)
            ax.axvline(x=pos, color="black", linewidth=1.0)

        offset = x_plot[-1] + MP_GAP

    ax.set_xticks(all_tick_positions)
    ax.set_xticklabels(all_tick_labels)
    ax.set_xlim(all_tick_positions[0], all_tick_positions[-1])
    ax.set_ylim(bottom=0)
    ax.set_title("MP-Style Path")
    ax.set_xlabel("Wave Vector")
    ax.tick_params(direction="in", top=True, right=True, length=6, width=1.1)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


plt.rcParams.update(
    {
        "font.size": 13,
        "axes.labelsize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.2,
    }
)

fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

plot_original(axes[0])
plot_mp_style(axes[1])

axes[0].set_ylabel(r"Frequencies (cm$^{-1}$)")

plt.tight_layout()
plt.savefig("Al2O3_compare_dispersion_paths.png", dpi=300, bbox_inches="tight")
plt.show()
