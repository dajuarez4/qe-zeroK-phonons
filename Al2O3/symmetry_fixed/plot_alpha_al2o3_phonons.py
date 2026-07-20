"""Plot the disconnected rhombohedral alpha-Al2O3 path and phonon DOS."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SEGMENTS = [
    (["Γ", "L", "B₁"], [0, 30, 60]),
    (["B", "Z", "Γ", "X"], [0, 30, 60, 90]),
    (["Q", "F", "P₁", "Z"], [0, 30, 60, 90]),
    (["L", "P"], [0, 30]),
]


def main() -> None:
    fig, ax = plt.subplots(figsize=(9.0, 6.5))
    offset = 0.0
    ticks: list[float] = []
    labels: list[str] = []

    for index, (segment_labels, tick_indices) in enumerate(SEGMENTS, start=1):
        data = np.loadtxt(ROOT / f"alpha_Al2O3.path{index}.freq.gp")
        x = data[:, 0] - data[0, 0] + offset
        for branch in data[:, 1:].T:
            ax.plot(x, branch, color="black", linewidth=0.9)
        for position, label in zip(x[tick_indices], segment_labels):
            if ticks and np.isclose(position, ticks[-1]):
                labels[-1] += f"|{label}"
            else:
                ticks.append(float(position))
                labels.append(label)
        offset = x[-1]

    for position in ticks:
        ax.axvline(position, color="0.35", linewidth=0.7)
    ax.axhline(0.0, color="black", linewidth=0.7)
    ax.set_xticks(ticks, labels)
    ax.set_xlim(ticks[0], ticks[-1])
    ax.set_xlabel("Wave vector")
    ax.set_ylabel(r"Frequency (cm$^{-1}$)")
    ax.set_title(r"Symmetry-fixed $\alpha$-Al$_2$O$_3$ phonon dispersion")
    ax.tick_params(direction="in")
    fig.tight_layout()
    fig.savefig(ROOT / "alpha_Al2O3_phonon_dispersion.png", dpi=300)
    plt.close(fig)

    dos_file = ROOT / "alpha_Al2O3.phdos.dat"
    if dos_file.exists():
        dos = np.loadtxt(dos_file)
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        ax.plot(dos[:, 0], dos[:, 1], color="black", linewidth=1.4)
        ax.axvline(0.0, color="black", linewidth=0.7)
        ax.set_xlabel(r"Frequency (cm$^{-1}$)")
        ax.set_ylabel("Density of states")
        ax.set_ylim(bottom=0.0)
        ax.tick_params(direction="in")
        fig.tight_layout()
        fig.savefig(ROOT / "alpha_Al2O3_phdos.png", dpi=300)
        plt.close(fig)


if __name__ == "__main__":
    main()
