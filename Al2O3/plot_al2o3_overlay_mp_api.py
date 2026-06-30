import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parent

QE_SEGMENTS = [
    {"file": ROOT / "Al2O3.mp_seg1.freq.gp", "labels": [r"\Gamma", "L", "B_1"], "tick_indices": [0, 30, 60]},
    {
        "file": ROOT / "Al2O3.mp_seg2.freq.gp",
        "labels": ["B", "Z", r"\Gamma", "X"],
        "tick_indices": [0, 30, 60, 90],
    },
    {
        "file": ROOT / "Al2O3.mp_seg3.freq.gp",
        "labels": ["Q", "F", "P_1", "Z"],
        "tick_indices": [0, 30, 60, 90],
    },
    {"file": ROOT / "Al2O3.mp_seg4.freq.gp", "labels": ["L", "P"], "tick_indices": [0, 30]},
]

MP_BAND_FILE = ROOT / "mp-1143-band.json"
MP_LABEL_SEQUENCES = [
    [r"\Gamma", "L", "B_1"],
    ["B", "Z", r"\Gamma", "X"],
    ["Q", "F", "P_1", "Z"],
    ["L", "P"],
]

CM1_PER_THZ = 33.35640952
GAP = 0.10
QE_COLOR = "black"
MP_COLOR = "#1f5aa6"
PROJECT_LABEL = "This work"
REFERENCE_LABEL = "Petretto et al. (2018)"
OUTPUT_PNG = ROOT / "Al2O3_qe_vs_mp_overlay.png"
OUTPUT_PDF = ROOT / "Al2O3_qe_vs_mp_overlay.pdf"


def format_label(label: str) -> str:
    return rf"${label}$"


def load_qe_segment(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected format in {path}")

    x = data[:, 0].astype(float)
    x -= x[0]
    freqs = data[:, 1:].astype(float)
    freqs[freqs < 0.0] = 0.0
    return x, freqs


def cumulative_distance(qpoints: np.ndarray) -> np.ndarray:
    diffs = np.diff(qpoints, axis=0)
    steps = np.linalg.norm(diffs, axis=1)
    return np.concatenate(([0.0], np.cumsum(steps)))


def find_label_indices(
    qpoints: np.ndarray, labels_dict: dict[str, list[float]], label_sequence: list[str], start_index: int
) -> list[int]:
    indices: list[int] = []
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


def load_mp_segments() -> list[dict[str, np.ndarray | list[str]]]:
    with MP_BAND_FILE.open() as f:
        data = json.load(f)

    qpoints = np.asarray(data["qpoints"], dtype=float)
    # The raw MP values are in THz-scale units; convert to cm^-1 so they match
    # the QE plot and the units shown on the MP website figure.
    freqs = np.asarray(data["frequencies"], dtype=float).T * CM1_PER_THZ
    freqs[freqs < 0.0] = 0.0
    labels_dict = data["labels_dict"]

    segments = []
    search_start = 0

    for label_sequence in MP_LABEL_SEQUENCES:
        idx = find_label_indices(qpoints, labels_dict, label_sequence, search_start)
        seg_q = qpoints[idx[0] : idx[-1] + 1]
        seg_f = freqs[idx[0] : idx[-1] + 1]
        seg_x = cumulative_distance(seg_q)
        tick_positions = [seg_x[i - idx[0]] for i in idx]

        segments.append(
            {
                "x": seg_x,
                "freqs": seg_f,
                "labels": label_sequence,
                "ticks": tick_positions,
            }
        )
        search_start = idx[-1] + 1

    return segments


def map_qe_x_to_mp(qe_x: np.ndarray, qe_tick_indices: list[int], mp_tick_positions: np.ndarray) -> np.ndarray:
    qe_tick_positions = qe_x[qe_tick_indices]
    return np.interp(qe_x, qe_tick_positions, mp_tick_positions)


def format_joint_label(left: str, right: str) -> str:
    return rf"${left}\,|\,{right}$"


plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 14,
        "axes.labelsize": 19,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "axes.linewidth": 1.25,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
    }
)

mp_segments = load_mp_segments()

if len(mp_segments) != len(QE_SEGMENTS):
    raise ValueError("Mismatch between QE segments and MP segments")

fig, ax = plt.subplots(figsize=(10.2, 6.6))

offset = 0.0
segment_tick_data: list[dict[str, list[float] | list[str]]] = []

for qe_seg, mp_seg in zip(QE_SEGMENTS, mp_segments):
    qe_x, qe_freqs = load_qe_segment(qe_seg["file"])
    mp_x = np.asarray(mp_seg["x"], dtype=float)
    mp_freqs = np.asarray(mp_seg["freqs"], dtype=float)
    labels = list(mp_seg["labels"])
    ticks = np.asarray(mp_seg["ticks"], dtype=float)

    qe_x_scaled = map_qe_x_to_mp(qe_x, qe_seg["tick_indices"], ticks)
    qe_x_plot = qe_x_scaled + offset
    mp_x_plot = mp_x + offset

    for mode in range(mp_freqs.shape[1]):
        ax.plot(
            mp_x_plot,
            mp_freqs[:, mode],
            color=MP_COLOR,
            linewidth=1.45,
            alpha=0.45,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=1,
        )

    for mode in range(qe_freqs.shape[1]):
        ax.plot(
            qe_x_plot,
            qe_freqs[:, mode],
            color=QE_COLOR,
            linewidth=0.95,
            alpha=0.92,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=2,
        )

    for pos, label in zip(ticks + offset, labels):
        ax.axvline(pos, color="0.70", linewidth=0.9, zorder=0)

    segment_tick_data.append(
        {
            "positions": [float(pos) for pos in (ticks + offset)],
            "labels": labels,
        }
    )

    offset += mp_x[-1] + GAP

tick_positions = [
    segment_tick_data[0]["positions"][0],
    segment_tick_data[0]["positions"][1],
]
tick_labels = [
    format_label(segment_tick_data[0]["labels"][0]),
    format_label(segment_tick_data[0]["labels"][1]),
]

for i in range(len(segment_tick_data) - 1):
    left = segment_tick_data[i]
    right = segment_tick_data[i + 1]
    seam_position = 0.5 * (left["positions"][-1] + right["positions"][0])
    tick_positions.append(seam_position)
    tick_labels.append(format_joint_label(left["labels"][-1], right["labels"][0]))

    for pos, label in zip(right["positions"][1:-1], right["labels"][1:-1]):
        tick_positions.append(pos)
        tick_labels.append(format_label(label))

tick_positions.append(segment_tick_data[-1]["positions"][-1])
tick_labels.append(format_label(segment_tick_data[-1]["labels"][-1]))

legend_handles = [
    Line2D([0], [0], color=QE_COLOR, linewidth=1.6, label=PROJECT_LABEL),
    Line2D([0], [0], color=MP_COLOR, linewidth=2.0, alpha=0.8, label=REFERENCE_LABEL),
]

ax.legend(
    handles=legend_handles,
    loc="upper right",
    bbox_to_anchor=(1.03, 0.99),
    frameon=False,
    handlelength=2.5,
    fontsize=12,
)
ax.set_ylabel(r"Frequency (cm$^{-1}$)")
ax.set_xlabel("Wave Vector")
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels)
ax.set_xlim(tick_positions[0], tick_positions[-1])
ax.set_ylim(0, 900)
ax.tick_params(top=True, right=True)
ax.spines["top"].set_visible(True)
ax.spines["right"].set_visible(True)

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight")
fig.savefig(OUTPUT_PDF, bbox_inches="tight")

print(f"saved {OUTPUT_PNG.name}")
print(f"saved {OUTPUT_PDF.name}")
