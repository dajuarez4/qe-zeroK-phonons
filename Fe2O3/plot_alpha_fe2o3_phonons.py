#!/usr/bin/env python3
"""Plot the four disconnected hematite phonon paths and the phonon DOS."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


SEGMENTS = [
    ("alpha_Fe2O3.path1.freq.gp", [0, 30, 60], [r"$\Gamma$", "L", r"B$_1$"]),
    ("alpha_Fe2O3.path2.freq.gp", [0, 30, 60, 90], ["B", "Z", r"$\Gamma$", "X"]),
    ("alpha_Fe2O3.path3.freq.gp", [0, 30, 60, 90], ["Q", "F", r"P$_1$", "Z"]),
    ("alpha_Fe2O3.path4.freq.gp", [0, 30], ["L", "P"]),
]


def main() -> None:
    missing = [name for name, _, _ in SEGMENTS if not Path(name).exists()]
    if missing:
        raise SystemExit("Missing matdyn files: " + ", ".join(missing))

    fig, axes = plt.subplots(
        1, 4, figsize=(11, 6), sharey=True,
        gridspec_kw={"width_ratios": [2, 3, 3, 1], "wspace": 0.04},
    )
    for ax, (name, ticks, labels) in zip(axes, SEGMENTS):
        data = np.loadtxt(name)
        x = np.arange(len(data))
        for branch in data[:, 1:].T:
            ax.plot(x, branch, color="firebrick", lw=0.9)
        ax.axhline(0, color="0.35", lw=0.7)
        for tick in ticks:
            ax.axvline(tick, color="0.75", lw=0.6)
        ax.set_xlim(0, len(data) - 1)
        ax.set_xticks(ticks, labels)
    axes[0].set_ylabel(r"Frequency (cm$^{-1}$)")
    fig.supxlabel("Wave vector")
    fig.suptitle(r"AFM $\alpha$-Fe$_2$O$_3$ phonon dispersion (PBE+U)")
    fig.subplots_adjust(bottom=0.12, top=0.91, left=0.09, right=0.98)
    fig.savefig("alpha_Fe2O3_phonon_dispersion.png", dpi=220)
    fig.savefig("alpha_Fe2O3_phonon_dispersion.pdf")
    plt.close(fig)

    dos_path = Path("alpha_Fe2O3.phdos.dat")
    if dos_path.exists():
        dos = np.loadtxt(dos_path)
        fig, ax = plt.subplots(figsize=(5.2, 6))
        ax.plot(dos[:, 1], dos[:, 0], color="firebrick", lw=1.2)
        ax.axhline(0, color="0.35", lw=0.7)
        ax.set_xlabel("Phonon DOS (states / cm$^{-1}$)")
        ax.set_ylabel(r"Frequency (cm$^{-1}$)")
        ax.set_title(r"AFM $\alpha$-Fe$_2$O$_3$ phonon DOS")
        fig.tight_layout()
        fig.savefig("alpha_Fe2O3_phdos.png", dpi=220)
        fig.savefig("alpha_Fe2O3_phdos.pdf")
        plt.close(fig)


if __name__ == "__main__":
    main()

