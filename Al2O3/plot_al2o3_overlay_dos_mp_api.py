import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parent

QE_DOS_FILE = ROOT / "Al2O3.phdos.dat"
MP_DOS_FILE = ROOT / "mp-1143-dos.json"

CM1_PER_THZ = 33.35640952
QE_COLOR = "black"
MP_COLOR = "#1f5aa6"
PROJECT_LABEL = "This work"
REFERENCE_LABEL = "Petretto et al. (2018)"

OUTPUT_PNG = ROOT / "Al2O3_qe_vs_mp_dos_overlay.png"
OUTPUT_PDF = ROOT / "Al2O3_qe_vs_mp_dos_overlay.pdf"


def load_qe_dos(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected DOS format in {path}")

    freq = data[:, 0].astype(float)
    dos = data[:, 1].astype(float)

    mask = freq >= 0.0
    return freq[mask], dos[mask]


def load_mp_dos(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open() as f:
        data = json.load(f)

    freq = np.asarray(data["frequencies"], dtype=float) * CM1_PER_THZ
    # The MP DOS is stored per THz. After converting the x-axis to cm^-1,
    # divide the DOS by the same factor to preserve the total number of modes.
    dos = np.asarray(data["densities"], dtype=float) / CM1_PER_THZ

    mask = freq >= 0.0
    return freq[mask], dos[mask]


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

qe_freq, qe_dos = load_qe_dos(QE_DOS_FILE)
mp_freq, mp_dos = load_mp_dos(MP_DOS_FILE)

fig, ax = plt.subplots(figsize=(8.4, 5.8))

ax.plot(
    mp_freq,
    mp_dos,
    color=MP_COLOR,
    linewidth=2.0,
    alpha=0.80,
    solid_capstyle="round",
    solid_joinstyle="round",
    zorder=1,
)
ax.plot(
    qe_freq,
    qe_dos,
    color=QE_COLOR,
    linewidth=1.55,
    solid_capstyle="round",
    solid_joinstyle="round",
    zorder=2,
)

legend_handles = [
    Line2D([0], [0], color=QE_COLOR, linewidth=1.7, label=PROJECT_LABEL),
    Line2D([0], [0], color=MP_COLOR, linewidth=2.3, label=REFERENCE_LABEL),
]

ax.legend(
    handles=legend_handles,
    loc="upper right",
    bbox_to_anchor=(1.03, 0.99),
    frameon=False,
    handlelength=2.6,
    fontsize=12,
)
ax.set_xlabel(r"Frequency (cm$^{-1}$)")
ax.set_ylabel(r"Phonon DOS")
ax.set_xlim(0, max(qe_freq.max(), mp_freq.max()))
ax.set_ylim(bottom=0.0)
ax.tick_params(top=True, right=True)
ax.spines["top"].set_visible(True)
ax.spines["right"].set_visible(True)

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight")
fig.savefig(OUTPUT_PDF, bbox_inches="tight")

print(f"saved {OUTPUT_PNG.name}")
print(f"saved {OUTPUT_PDF.name}")
