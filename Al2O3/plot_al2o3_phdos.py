from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


dos_file = Path("Al2O3.phdos.dat")

data = np.loadtxt(dos_file)
if data.ndim != 2 or data.shape[1] < 2:
    raise ValueError(f"Unexpected DOS format in {dos_file}")

freq = data[:, 0]
dos = data[:, 1]

plt.rcParams.update(
    {
        "font.size": 14,
        "axes.labelsize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "axes.linewidth": 1.2,
    }
)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(freq, dos, color="black", linewidth=1.5)
ax.set_xlabel(r"Frequency (cm$^{-1}$)")
ax.set_ylabel(r"Density of States")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.tick_params(direction="in", length=6, width=1.2)
ax.spines["top"].set_visible(True)
ax.spines["right"].set_visible(True)

plt.tight_layout()
plt.savefig("Al2O3_phdos.png", dpi=300, bbox_inches="tight")
plt.show()
