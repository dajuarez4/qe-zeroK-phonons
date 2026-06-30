import json
import csv
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
MP_DOS_FILE = ROOT / "mp-1143-dos.json"
QE_DOS_FILE = ROOT / "Al2O3.phdos.dat"
EXPERIMENT_CSV = ROOT / "sapphire_20K_experiment_points.csv"

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
EXPERIMENT_LABEL = "Schober et al. (1993)"

OUTPUT_PNG = ROOT / "Al2O3_qe_vs_mp_band_dos.png"
OUTPUT_PDF = ROOT / "Al2O3_qe_vs_mp_band_dos.pdf"


def format_label(label: str) -> str:
    return rf"${label}$"


def format_joint_label(left: str, right: str) -> str:
    return rf"${left}\,|\,{right}$"


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


def build_dispersion_ticks(
    segment_tick_data: list[dict[str, list[float] | list[str]]],
) -> tuple[list[float], list[str]]:
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
    return tick_positions, tick_labels


def load_qe_dos() -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(QE_DOS_FILE, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected DOS format in {QE_DOS_FILE}")

    freq = data[:, 0].astype(float)
    dos = data[:, 1].astype(float)
    mask = freq >= 0.0
    return freq[mask], dos[mask]


def load_mp_dos() -> tuple[np.ndarray, np.ndarray]:
    with MP_DOS_FILE.open() as f:
        data = json.load(f)

    freq = np.asarray(data["frequencies"], dtype=float) * CM1_PER_THZ
    dos = np.asarray(data["densities"], dtype=float) / CM1_PER_THZ
    mask = freq >= 0.0
    return freq[mask], dos[mask]


def load_experiment_gamma_l() -> tuple[np.ndarray, np.ndarray]:
    xi_values: list[float] = []
    freqs_cm1: list[float] = []

    with EXPERIMENT_CSV.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["segment"] != "gamma_a":
                continue

            xi_values.append(float(row["xi"]))
            freqs_cm1.append(float(row["frequency_thz"]) * CM1_PER_THZ)

    return np.asarray(xi_values, dtype=float), np.asarray(freqs_cm1, dtype=float)


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
qe_dos_freq, qe_dos = load_qe_dos()
mp_dos_freq, mp_dos = load_mp_dos()
exp_gamma_l_xi, exp_gamma_l_freq = load_experiment_gamma_l()

if len(mp_segments) != len(QE_SEGMENTS):
    raise ValueError("Mismatch between QE segments and MP segments")

band_frequency_max = max(float(np.max(np.asarray(seg["freqs"], dtype=float))) for seg in mp_segments)
frequency_max = max(band_frequency_max, float(np.max(qe_dos_freq)), float(np.max(mp_dos_freq)))
frequency_max = 50.0 * np.ceil(frequency_max / 50.0)
dos_max = 1.08 * max(float(np.max(qe_dos)), float(np.max(mp_dos)))

fig, (ax_band, ax_dos) = plt.subplots(
    1,
    2,
    figsize=(12.0, 6.9),
    sharey=True,
    gridspec_kw={"width_ratios": [5.4, 1.45], "wspace": 0.05},
)

offset = 0.0
segment_tick_data: list[dict[str, list[float] | list[str]]] = []

for segment_index, (qe_seg, mp_seg) in enumerate(zip(QE_SEGMENTS, mp_segments)):
    qe_x, qe_freqs = load_qe_segment(qe_seg["file"])
    mp_x = np.asarray(mp_seg["x"], dtype=float)
    mp_freqs = np.asarray(mp_seg["freqs"], dtype=float)
    labels = list(mp_seg["labels"])
    ticks = np.asarray(mp_seg["ticks"], dtype=float)

    qe_x_scaled = map_qe_x_to_mp(qe_x, qe_seg["tick_indices"], ticks)
    qe_x_plot = qe_x_scaled + offset
    mp_x_plot = mp_x + offset

    for mode in range(mp_freqs.shape[1]):
        ax_band.plot(
            mp_x_plot,
            mp_freqs[:, mode],
            color=MP_COLOR,
            linewidth=1.45,
            alpha=0.42,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=1,
        )

    for mode in range(qe_freqs.shape[1]):
        ax_band.plot(
            qe_x_plot,
            qe_freqs[:, mode],
            color=QE_COLOR,
            linewidth=0.95,
            alpha=0.92,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=2,
        )

    if segment_index == 0:
        exp_x_plot = ticks[0] + (exp_gamma_l_xi / 0.5) * (ticks[1] - ticks[0]) + offset
        ax_band.scatter(
            exp_x_plot,
            exp_gamma_l_freq,
            s=18,
            facecolors="black",
            edgecolors="white",
            linewidths=0.22,
            alpha=0.96,
            zorder=3,
        )

    for pos in ticks + offset:
        ax_band.axvline(pos, color="0.72", linewidth=0.9, zorder=0)

    segment_tick_data.append(
        {
            "positions": [float(pos) for pos in (ticks + offset)],
            "labels": labels,
        }
    )

    offset += mp_x[-1] + GAP

tick_positions, tick_labels = build_dispersion_ticks(segment_tick_data)

ax_dos.plot(
    mp_dos,
    mp_dos_freq,
    color=MP_COLOR,
    linewidth=2.0,
    alpha=0.80,
    solid_capstyle="round",
    solid_joinstyle="round",
    zorder=1,
)
ax_dos.plot(
    qe_dos,
    qe_dos_freq,
    color=QE_COLOR,
    linewidth=1.55,
    solid_capstyle="round",
    solid_joinstyle="round",
    zorder=2,
)

legend_handles = [
    Line2D([0], [0], color=QE_COLOR, linewidth=1.7, label=PROJECT_LABEL),
    Line2D([0], [0], color=MP_COLOR, linewidth=2.2, label=REFERENCE_LABEL),
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor="black",
        markeredgecolor="white",
        markeredgewidth=0.35,
        markersize=5.5,
        label=EXPERIMENT_LABEL,
    ),
]

ax_band.set_ylabel(r"Frequency (cm$^{-1}$)")
ax_band.set_xlabel("Wave Vector")
ax_band.set_xticks(tick_positions)
ax_band.set_xticklabels(tick_labels)
ax_band.set_xlim(tick_positions[0], tick_positions[-1])
ax_band.set_ylim(0.0, frequency_max)
ax_band.tick_params(top=True, right=False)
ax_band.spines["top"].set_visible(True)
ax_band.spines["right"].set_visible(False)

ax_dos.set_xlabel("Phonon DOS")
ax_dos.set_xlim(0.0, dos_max)
ax_dos.tick_params(top=True, right=True, left=False, labelleft=False)
ax_dos.spines["top"].set_visible(True)
ax_dos.spines["left"].set_visible(False)
ax_dos.spines["right"].set_visible(True)
ax_dos.legend(
    handles=legend_handles,
    loc="upper right",
    bbox_to_anchor=(1.05, 0.99),
    frameon=False,
    handlelength=2.5,
    fontsize=12,
)

fig.subplots_adjust(left=0.082, right=0.988, bottom=0.145, top=0.985, wspace=0.05)
fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight")
fig.savefig(OUTPUT_PDF, bbox_inches="tight")

print(f"saved {OUTPUT_PNG.name}")
print(f"saved {OUTPUT_PDF.name}")
