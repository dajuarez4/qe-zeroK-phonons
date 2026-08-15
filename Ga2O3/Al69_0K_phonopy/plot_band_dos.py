#!/usr/bin/env python3
"""Plot the 69% Al alloy phonon band structure and total DOS."""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent

with h5py.File(ROOT / "band.hdf5") as band:
    distance = band["distance"][:]
    frequency = band["frequency"][:]
    labels = band["label"][:]

dos = np.loadtxt(ROOT / "total_dos.dat")
fig, (ax, dax) = plt.subplots(
    1, 2, figsize=(8.2, 5.0), sharey=True,
    gridspec_kw={"width_ratios": [4.5, 1.2], "wspace": 0.08},
)

for segment in range(distance.shape[0]):
    ax.plot(distance[segment], frequency[segment], color="#185fa5", lw=0.75)

ticks = [distance[0, 0], *[x[-1] for x in distance]]
ticklabels = [labels[0, 0], *labels[:, 1]]
ticklabels = [x.decode() if isinstance(x, bytes) else str(x) for x in ticklabels]
ax.set_xticks(ticks, ticklabels)
for tick in ticks:
    ax.axvline(tick, color="0.82", lw=0.6)
ax.axhline(0, color="black", ls=":", lw=0.8)
ax.set_xlim(ticks[0], ticks[-1])
ax.set_xlabel("Wave vector")
ax.set_ylabel("Frequency (THz)")

dax.plot(dos[:, 1], dos[:, 0], color="#d1495b", lw=1.2)
dax.axhline(0, color="black", ls=":", lw=0.8)
dax.set_xlim(left=0)
dax.set_xlabel("DOS")

ymin = min(-0.5, float(np.nanmin(frequency)) - 0.2)
ymax = max(float(np.nanmax(frequency)), float(np.nanmax(dos[:, 0]))) + 0.3
ax.set_ylim(ymin, ymax)
fig.suptitle(r"$(\mathrm{Al}_{0.6875}\mathrm{Ga}_{0.3125})_2\mathrm{O}_3$ harmonic phonons at 0 K")
fig.subplots_adjust(top=0.90, bottom=0.14, left=0.10, right=0.98)
fig.savefig(ROOT / "Al69_Ga2O3_phonon_band_dos.pdf")
fig.savefig(ROOT / "Al69_Ga2O3_phonon_band_dos.png", dpi=220)
plt.close(fig)
