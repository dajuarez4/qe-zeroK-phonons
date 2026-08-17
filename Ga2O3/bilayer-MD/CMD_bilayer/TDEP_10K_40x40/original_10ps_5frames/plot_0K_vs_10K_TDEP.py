#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FD = HERE.parent / "phonopy_0K_tabgap_40x40_amp002" / "dispersion_0K_tabgap.dat"
TDEP = HERE / "outfile.dispersion_relations"
OUT = HERE.parent / "phonon_comparison_0K_FD_vs_10K_TDEP_40x40.png"


def read_table(path):
    rows = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        try:
            rows.append([float(x) for x in s.split()])
        except ValueError:
            pass
    a = np.asarray(rows)
    if a.ndim != 2 or a.shape[1] < 2:
        raise RuntimeError(f"No numeric dispersion table found in {path}")
    return a


def segmented_x(n):
    # Both calculations use four equal path sections: Gamma-X-S-Y-Gamma.
    return np.linspace(0.0, 4.0, n)

fd = read_table(FD)
td = read_table(TDEP)
xf, xt = segmented_x(len(fd)), segmented_x(len(td))
yf, yt = fd[:, 1:], td[:, 1:]

fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                         gridspec_kw={"height_ratios": [2.1, 1]})

for ax in axes:
    ax.axhline(0, color="black", lw=0.9, zorder=1)
    for x in range(5):
        ax.axvline(x, color="0.82", lw=0.8, zorder=0)
    for j in range(yf.shape[1]):
        ax.plot(xf, yf[:, j], color="#2563eb", lw=0.55, alpha=0.70)
    for j in range(yt.shape[1]):
        ax.plot(xt, yt[:, j], color="#dc2626", lw=0.55, alpha=0.70)
    ax.grid(axis="y", color="0.92", lw=0.6)

axes[0].set_ylabel("Frequency (THz)")
axes[0].set_title("16,000-atom Ga$_2$O$_3$ bilayer: harmonic phonons\n"
                  "0 K finite displacement vs 10 K TDEP (cutoff 5.5 Angstrom)")
axes[0].plot([], [], color="#2563eb", lw=1.5,
             label=f"0 K finite displacement (min {yf.min():.3f} THz)")
axes[0].plot([], [], color="#dc2626", lw=1.5,
             label=f"10 K TDEP (min {yt.min():.3f} THz; fit $R^2$=0.810)")
axes[0].legend(loc="upper right", frameon=True, fontsize=9)

axes[1].set_ylim(-0.55, 3.5)
axes[1].set_ylabel("Frequency (THz)")
axes[1].set_xlabel("Wave-vector path")
axes[1].set_title("Low-frequency zoom")
axes[1].fill_between([0, 4], -0.55, 0, color="#fee2e2", alpha=0.45, zorder=-1)
axes[1].set_xticks(range(5), [r"$\Gamma$", "X", "S", "Y", r"$\Gamma$"])

fig.tight_layout()
fig.savefig(OUT, dpi=240, bbox_inches="tight")
print(OUT)
