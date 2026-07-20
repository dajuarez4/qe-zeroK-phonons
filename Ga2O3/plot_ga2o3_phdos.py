"""Plot the total Ga2O3 phonon DOS written by matdyn.x."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
dos_file = ROOT / "Ga2O3.phdos.dat"
data = np.loadtxt(dos_file)
if data.ndim != 2 or data.shape[1] < 2:
    raise ValueError(f"Unexpected DOS format in {dos_file}: {data.shape}")

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(data[:, 0], data[:, 1], color="black", linewidth=1.5)
ax.axvline(0.0, color="#b22222", linewidth=0.9, linestyle="--")
ax.set_xlabel(r"Frequency (cm$^{-1}$)")
ax.set_ylabel("Density of states")
ax.set_ylim(bottom=0.0)
ax.tick_params(direction="in")
fig.tight_layout()
fig.savefig(ROOT / "Ga2O3_phdos.png", dpi=300, bbox_inches="tight")
