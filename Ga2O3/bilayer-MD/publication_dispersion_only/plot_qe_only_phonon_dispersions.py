"""Publication-style QE-only phonon dispersions for beta-Ga2O3 and alpha-Al2O3."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


OUTPUT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = OUTPUT_DIR.parents[2]
GA_ROOT = PROJECT_ROOT / "Ga2O3" / "ideal_beta"
AL_ROOT = PROJECT_ROOT / "Al2O3"
TDEP_ROOT = (
    PROJECT_ROOT
    / "Ga2O3"
    / "bilayer-MD"
    / "Ga2O3-BL-001-4A-relax-4x2x1-300K"
    / "TDEP_300K"
)

CM1_PER_THZ = 33.35640952
GAP_FRACTION = 0.10

GA_SEGMENTS = [
    {"file": GA_ROOT / "beta_Ga2O3.path1.freq.gp", "labels": [r"\Gamma", "Y", "F", "L", "I"], "ticks": [0, 30, 60, 90, 120]},
    {"file": GA_ROOT / "beta_Ga2O3.path2.freq.gp", "labels": ["I_1", "Z", "F_1"], "ticks": [0, 30, 60]},
    {"file": GA_ROOT / "beta_Ga2O3.path3.freq.gp", "labels": ["Y", "X_1"], "ticks": [0, 30]},
    {"file": GA_ROOT / "beta_Ga2O3.path4.freq.gp", "labels": ["X", r"\Gamma", "N"], "ticks": [0, 30, 60]},
    {"file": GA_ROOT / "beta_Ga2O3.path5.freq.gp", "labels": ["M", r"\Gamma"], "ticks": [0, 30]},
]

AL_SEGMENTS = [
    {"file": AL_ROOT / "Al2O3.mp_seg1.freq.gp", "labels": [r"\Gamma", "L", "B_1"], "ticks": [0, 30, 60]},
    {"file": AL_ROOT / "Al2O3.mp_seg2.freq.gp", "labels": ["B", "Z", r"\Gamma", "X"], "ticks": [0, 30, 60, 90]},
    {"file": AL_ROOT / "Al2O3.mp_seg3.freq.gp", "labels": ["Q", "F", "P_1", "Z"], "ticks": [0, 30, 60, 90]},
    {"file": AL_ROOT / "Al2O3.mp_seg4.freq.gp", "labels": ["L", "P"], "ticks": [0, 30]},
]


plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 22,
        "axes.labelsize": 30,
        "xtick.labelsize": 25,
        "ytick.labelsize": 24,
        "axes.linewidth": 1.25,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
    }
)


def format_label(label: str) -> str:
    return rf"${label}$"


def format_joint_label(left: str, right: str) -> str:
    return rf"${left}\,|\,{right}$"


def load_segment(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected phonon data format in {path}")
    x = data[:, 0].astype(float)
    x -= x[0]
    frequencies = data[:, 1:].astype(float) / CM1_PER_THZ
    # Match the reference figures, which display tiny Gamma numerical noise at zero.
    return x, np.maximum(frequencies, 0.0)


def build_ticks(segment_tick_data: list[dict]) -> tuple[list[float], list[str]]:
    positions = list(segment_tick_data[0]["positions"][:-1])
    labels = [format_label(label) for label in segment_tick_data[0]["labels"][:-1]]

    for left, right in zip(segment_tick_data[:-1], segment_tick_data[1:]):
        positions.append(0.5 * (left["positions"][-1] + right["positions"][0]))
        labels.append(format_joint_label(left["labels"][-1], right["labels"][0]))
        for position, label in zip(right["positions"][1:-1], right["labels"][1:-1]):
            positions.append(position)
            labels.append(format_label(label))

    positions.append(segment_tick_data[-1]["positions"][-1])
    labels.append(format_label(segment_tick_data[-1]["labels"][-1]))
    return positions, labels


def plot_dispersion(segments: list[dict], stem: str, frequency_max: float) -> None:
    loaded = [(segment, *load_segment(segment["file"])) for segment in segments]
    characteristic_width = np.median(
        [x[-1] / max(len(segment["labels"]) - 1, 1) for segment, x, _ in loaded]
    )
    gap = GAP_FRACTION * characteristic_width

    fig, ax = plt.subplots(figsize=(13.0, 8.2))
    offset = 0.0
    tick_data: list[dict] = []

    for segment, x, frequencies in loaded:
        x_plot = x + offset
        tick_positions = x[segment["ticks"]] + offset

        for branch in frequencies.T:
            ax.plot(
                x_plot,
                branch,
                color="black",
                linewidth=0.95,
                alpha=0.92,
                solid_capstyle="round",
                solid_joinstyle="round",
                zorder=2,
            )

        for position in tick_positions:
            ax.axvline(position, color="0.72", linewidth=0.9, zorder=0)

        tick_data.append(
            {
                "positions": [float(value) for value in tick_positions],
                "labels": list(segment["labels"]),
            }
        )
        offset += x[-1] + gap

    tick_positions, tick_labels = build_ticks(tick_data)
    ax.axhline(0.0, color="0.35", linewidth=0.8, zorder=1)
    ax.set_ylabel("Frequency (THz)")
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(tick_positions[0], tick_positions[-1])
    ax.set_ylim(0.0, frequency_max)
    ax.tick_params(axis="x", labelsize=27, pad=10, top=True)
    ax.tick_params(axis="y", labelsize=25, right=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)

    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.15, top=0.975)
    png_path = OUTPUT_DIR / f"{stem}.png"
    pdf_path = OUTPUT_DIR / f"{stem}.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {png_path.name}")
    print(f"saved {pdf_path.name}")


def plot_tdep_bilayer() -> None:
    """Plot the latest computed 4-Angstrom bilayer TDEP dispersion."""
    data = np.loadtxt(TDEP_ROOT / "outfile.dispersion_relations")
    x = data[:, 0]
    frequencies = data[:, 1:]
    tick_indices = [0, 79, 159, 239, 319]
    tick_positions = x[tick_indices]
    tick_labels = [r"$\Gamma$", r"$X$", r"$S$", r"$Y$", r"$\Gamma$"]

    fig, ax = plt.subplots(figsize=(13.0, 8.2))
    for position in tick_positions:
        ax.axvline(position, color="0.72", linewidth=0.9, zorder=0)
    for branch in frequencies.T:
        ax.plot(
            x,
            branch,
            color="black",
            linewidth=0.95,
            alpha=0.92,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=2,
        )

    ax.axhline(0.0, color="0.35", linewidth=0.8, zorder=1)
    ax.set_ylabel("Frequency (THz)")
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-5.0, 20.0)
    ax.tick_params(axis="x", labelsize=27, pad=10, top=True)
    ax.tick_params(axis="y", labelsize=25, right=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)

    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.15, top=0.975)
    stem = "Ga2O3_4A_bilayer_300K_TDEP_this_work_phonon_dispersion"
    png_path = OUTPUT_DIR / f"{stem}.png"
    pdf_path = OUTPUT_DIR / f"{stem}.pdf"
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {png_path.name}")
    print(f"saved {pdf_path.name}")


if __name__ == "__main__":
    plot_dispersion(GA_SEGMENTS, "beta_Ga2O3_this_work_phonon_dispersion", 30.0)
    plot_dispersion(AL_SEGMENTS, "alpha_Al2O3_this_work_phonon_dispersion", 30.0)
    plot_tdep_bilayer()
