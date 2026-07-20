"""Validate and plot the Ga2O3 phonons produced by Quantum ESPRESSO."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FREQ_FILE = ROOT / "Ga2O3.freq.gp"
DISPERSION_PLOT = ROOT / "Ga2O3_phonon_dispersion.png"
GAMMA_PLOT = ROOT / "Ga2O3_Gamma_phonons_sticks.png"
GAMMA_TABLE = ROOT / "Ga2O3_Gamma_freq.dat"
SUMMARY_FILE = ROOT / "Ga2O3_phonon_summary.txt"

N_MODES = 30
POINTS_PER_SEGMENT = 40
PATH_LABELS = [
    r"$\Gamma$",
    r"$q_1$",
    r"$q_2$",
    r"$\Gamma$",
    r"$q_3$",
    r"$q_4$",
    r"$q_5$",
    r"$q_3$",
]


def load_frequencies() -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(FREQ_FILE)
    if data.ndim != 2 or data.shape[1] != N_MODES + 1:
        raise ValueError(
            f"Expected distance plus {N_MODES} modes in {FREQ_FILE}, "
            f"found shape {data.shape}"
        )
    return data[:, 0], data[:, 1:]


def path_ticks(q_distance: np.ndarray) -> tuple[list[int], np.ndarray]:
    expected_points = POINTS_PER_SEGMENT * (len(PATH_LABELS) - 1) + 1
    if len(q_distance) != expected_points:
        raise ValueError(
            f"Expected {expected_points} q points for the path in Ga2O3.matdyn.in, "
            f"found {len(q_distance)}"
        )
    indices = [
        min(segment * POINTS_PER_SEGMENT, len(q_distance) - 1)
        for segment in range(len(PATH_LABELS))
    ]
    return indices, q_distance[indices]


def write_gamma_table(gamma: np.ndarray) -> None:
    with GAMMA_TABLE.open("w", encoding="utf-8") as handle:
        handle.write("# mode frequency_cm-1 frequency_THz imaginary\n")
        for mode, value in enumerate(gamma, start=1):
            handle.write(
                f"{mode:2d} {value:12.4f} {value / 33.35640952:12.6f} "
                f"{str(value < 0).lower()}\n"
            )


def write_summary(q_distance: np.ndarray, frequencies: np.ndarray) -> str:
    minimum_index = np.unravel_index(np.argmin(frequencies), frequencies.shape)
    negative = frequencies < 0.0
    gamma = frequencies[0]
    gamma_optical_instabilities = gamma[np.abs(gamma) > 10.0]
    gamma_optical_instabilities = gamma_optical_instabilities[
        gamma_optical_instabilities < 0.0
    ]

    summary = "\n".join(
        [
            "Ga2O3 phonon validation summary",
            "================================",
            f"q points: {len(q_distance)}",
            f"branches: {frequencies.shape[1]}",
            f"minimum frequency: {frequencies[minimum_index]:.4f} cm^-1",
            f"minimum at path distance: {q_distance[minimum_index[0]]:.6f}",
            f"minimum branch: {minimum_index[1] + 1}",
            f"maximum frequency: {frequencies.max():.4f} cm^-1",
            f"negative sampled frequencies: {negative.sum()}",
            f"q points with at least one negative mode: {negative.any(axis=1).sum()}",
            "Gamma frequencies below -10 cm^-1: "
            + (
                ", ".join(f"{value:.4f}" for value in gamma_optical_instabilities)
                if len(gamma_optical_instabilities)
                else "none"
            ),
            "",
            "Interpretation: the supplied force constants contain imaginary modes",
            "and do not reproduce the expected stable beta-Ga2O3 dispersion. This",
            "indicates a structure/symmetry/convergence mismatch in this calculation;",
            "it is not evidence that beta-Ga2O3 itself is unstable.",
            "",
        ]
    )
    SUMMARY_FILE.write_text(summary, encoding="utf-8")
    return summary


def plot_dispersion(q_distance: np.ndarray, frequencies: np.ndarray) -> None:
    _, ticks = path_ticks(q_distance)
    fig, ax = plt.subplots(figsize=(10.5, 6.8))
    for branch in frequencies.T:
        ax.plot(q_distance, branch, color="black", linewidth=0.75)
    for position in ticks:
        ax.axvline(position, color="0.72", linewidth=0.7, zorder=0)
    ax.axhspan(frequencies.min() - 5.0, 0.0, color="#d62728", alpha=0.08)
    ax.axhline(0.0, color="#b22222", linewidth=0.9, linestyle="--")
    ax.set_xticks(ticks, PATH_LABELS)
    ax.set_xlim(q_distance[0], q_distance[-1])
    ax.set_ylim(frequencies.min() - 5.0, frequencies.max() + 15.0)
    ax.set_xlabel("Wave vector along the specified reciprocal-coordinate path")
    ax.set_ylabel(r"Frequency (cm$^{-1}$)")
    ax.set_title(r"Ga$_2$O$_3$ phonon dispersion (2$\times$2$\times$2 q grid)")
    ax.tick_params(direction="in")
    ax.text(
        0.5,
        -0.16,
        r"$q_1=(1/2,0,0)$; $q_2=(1/3,1/3,0)$; "
        r"$q_3=(0,0,1/2)$; $q_4=(1/2,0,1/2)$; "
        r"$q_5=(1/3,1/3,1/2)$",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=9,
    )
    fig.subplots_adjust(bottom=0.23)
    fig.savefig(DISPERSION_PLOT, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_gamma(gamma: np.ndarray) -> None:
    modes = np.arange(1, len(gamma) + 1)
    colors = np.where(gamma < 0.0, "#c62828", "black")
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.vlines(modes, 0.0, gamma, color=colors, linewidth=1.5)
    ax.scatter(modes, gamma, c=colors, s=18, zorder=3)
    ax.axhline(0.0, color="0.35", linewidth=0.8)
    ax.set_xlim(0.3, len(gamma) + 0.7)
    ax.set_xlabel("Mode index")
    ax.set_ylabel(r"Frequency (cm$^{-1}$)")
    ax.set_title(r"Ga$_2$O$_3$ frequencies at $\Gamma$ after crystal ASR")
    ax.tick_params(direction="in")
    fig.tight_layout()
    fig.savefig(GAMMA_PLOT, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    q_distance, frequencies = load_frequencies()
    write_gamma_table(frequencies[0])
    plot_dispersion(q_distance, frequencies)
    plot_gamma(frequencies[0])
    print(write_summary(q_distance, frequencies), end="")


if __name__ == "__main__":
    main()
