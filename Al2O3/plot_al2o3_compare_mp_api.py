import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


LOCAL_SEGMENTS = [
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

MP_RAW_FILE = Path("mp-1143-band.json")
MP_LABEL_SEQUENCES = [
    [r"\Gamma", "L", "B_1"],
    ["B", "Z", r"\Gamma", "X"],
    ["Q", "F", "P_1", "Z"],
    ["L", "P"],
]

GAP = 0.08


def load_gp(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected format in {path}")
    x = data[:, 0].copy()
    freqs = data[:, 1:].copy()
    freqs[freqs < 0.0] = 0.0
    return x, freqs


def plot_local_qe(ax: plt.Axes) -> None:
    offset = 0.0
    tick_positions: list[float] = []
    tick_labels: list[str] = []

    for seg in LOCAL_SEGMENTS:
        x_local, freqs = load_gp(seg["file"])
        x_plot = x_local + offset

        for mode in range(freqs.shape[1]):
            ax.plot(x_plot, freqs[:, mode], color="black", linewidth=1.1)

        for idx, label in zip(seg["tick_indices"], seg["labels"]):
            if idx >= len(x_plot):
                continue
            pos = x_plot[idx]
            tick_positions.append(pos)
            tick_labels.append(label)
            ax.axvline(pos, color="gray", linewidth=1.0)

        offset = x_plot[-1] + GAP

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(tick_positions[0], tick_positions[-1])
    ax.set_ylim(bottom=0)
    ax.set_title("QE On MP-Style Path")
    ax.set_xlabel("Wave Vector")
    ax.tick_params(direction="in", top=True, right=True, length=6, width=1.1)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


def find_label_indices(qpoints, labels_dict, label_sequence, start_index):
    indices = []
    start = start_index
    for label in label_sequence:
        target = labels_dict[label]
        match = None
        for i in range(start, len(qpoints)):
            if all(abs(a - b) < 1.0e-8 for a, b in zip(qpoints[i], target)):
                match = i
                break
        if match is None:
            raise ValueError(f"Could not locate label {label} in MP band data")
        indices.append(match)
        start = match + 1
    return indices


def cumulative_distance(qpoints: np.ndarray) -> np.ndarray:
    diffs = np.diff(qpoints, axis=0)
    steps = np.linalg.norm(diffs, axis=1)
    return np.concatenate(([0.0], np.cumsum(steps)))


def plot_mp_raw(ax: plt.Axes) -> None:
    with MP_RAW_FILE.open() as f:
        data = json.load(f)

    qpoints = np.asarray(data["qpoints"], dtype=float)
    freqs = np.asarray(data["frequencies"], dtype=float).T
    freqs[freqs < 0.0] = 0.0
    labels_dict = data["labels_dict"]

    offset = 0.0
    search_start = 0
    tick_positions: list[float] = []
    tick_labels: list[str] = []

    for label_sequence in MP_LABEL_SEQUENCES:
        idx = find_label_indices(qpoints, labels_dict, label_sequence, search_start)
        seg_q = qpoints[idx[0] : idx[-1] + 1]
        seg_f = freqs[idx[0] : idx[-1] + 1]
        x_local = cumulative_distance(seg_q)
        x_plot = x_local + offset

        for mode in range(seg_f.shape[1]):
            ax.plot(x_plot, seg_f[:, mode], color="blue", linewidth=1.1)

        for label, i_global in zip(label_sequence, idx):
            pos = x_plot[i_global - idx[0]]
            tick_positions.append(pos)
            tick_labels.append(
                r"$\Gamma$" if label == r"\Gamma" else rf"${label}$" if "_" in label else label
            )
            ax.axvline(pos, color="black", linewidth=1.0)

        offset = x_plot[-1] + GAP
        search_start = idx[-1] + 1

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(tick_positions[0], tick_positions[-1])
    ax.set_ylim(bottom=0)
    ax.set_title("Materials Project Raw Band")
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

plot_local_qe(axes[0])
plot_mp_raw(axes[1])

axes[0].set_ylabel(r"Frequencies (cm$^{-1}$)")

plt.tight_layout()
plt.savefig("Al2O3_compare_with_mp_api.png", dpi=300, bbox_inches="tight")
plt.show()
