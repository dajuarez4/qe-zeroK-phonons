#!/usr/bin/env python3
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
SOURCES = {
    0: ROOT / "phonopy_0K_tabgap_40x40_amp002" / "dispersion_0K_tabgap.dat",
    10: ROOT / "tdep_10K_40x40_extended" / "outfile.dispersion_relations",
    20: ROOT / "tdep_20K_40x40" / "outfile.dispersion_relations",
    30: ROOT / "tdep_30K_40x40" / "outfile.dispersion_relations",
    40: ROOT / "tdep_40K_40x40" / "outfile.dispersion_relations",
    50: ROOT / "tdep_50K_40x40" / "outfile.dispersion_relations",
    60: ROOT / "tdep_60K_40x40" / "outfile.dispersion_relations",
    70: ROOT / "tdep_70K_40x40" / "outfile.dispersion_relations",
    80: ROOT / "tdep_80K_40x40" / "outfile.dispersion_relations",
    90: ROOT / "tdep_90K_40x40" / "outfile.dispersion_relations",
    100: ROOT / "tdep_100K_40x40" / "outfile.dispersion_relations",
}
# The TDEP paths contain 100 rows per segment. The 0 K path differs slightly
# in sampling, so fractional path positions are used rather than fixed indices.
POINTS = {"Gamma": 0.0, "X": 0.25, "S": 0.50, "Y": 0.75}


def frequencies(path):
    rows = []
    for line in path.read_text().splitlines():
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            pass
    table = np.asarray(rows)
    return table[:, 1:]


temperatures = np.array(sorted(SOURCES))
values = {point: [] for point in POINTS}
values["Global minimum"] = []

for temperature in temperatures:
    modes = frequencies(SOURCES[int(temperature)])
    last = len(modes) - 1
    for point, fraction in POINTS.items():
        index = int(round(fraction * last))
        values[point].append(float(modes[index].min()))
    values["Global minimum"].append(float(modes.min()))

csv_path = ROOT / "phonon_stability_vs_temperature.csv"
with csv_path.open("w", newline="") as output:
    writer = csv.writer(output)
    writer.writerow(["temperature_K", "Gamma_THz", "X_THz", "S_THz", "Y_THz", "global_min_THz"])
    for i, temperature in enumerate(temperatures):
        writer.writerow([temperature] + [values[key][i] for key in ("Gamma", "X", "S", "Y", "Global minimum")])

colors = {"Gamma": "#7c3aed", "X": "#2563eb", "S": "#16a34a", "Y": "#dc2626"}
markers = {"Gamma": "o", "X": "s", "S": "^", "Y": "D"}
fig, (ax, zoom) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                               gridspec_kw={"height_ratios": [1.35, 1]})

for axis in (ax, zoom):
    axis.axhspan(-0.5, 0, color="#fee2e2", alpha=0.65, label="Imaginary / unstable" if axis is ax else None)
    axis.axhline(0, color="black", lw=1.1)
    for point in ("Gamma", "X", "S", "Y"):
        axis.plot(temperatures, values[point], color=colors[point], marker=markers[point],
                  ms=6.5, lw=2, label=point if axis is ax else None)
    axis.plot(temperatures, values["Global minimum"], color="#111827", marker="x",
              ms=7, lw=1.6, ls="--", label="Global minimum" if axis is ax else None)
    axis.grid(True, color="0.90", lw=0.7)

ax.set_ylabel("Lowest frequency (THz)")
ax.set_title("Ga$_2$O$_3$ bilayer phonon stability versus temperature\n"
             "40x40x1 supercell; 10-100 K TDEP cutoff 5.5 Angstrom")
ax.legend(ncol=3, fontsize=9, loc="best")
ax.set_ylim(-0.45, 2.55)

zoom.set_ylim(-0.30, 0.75)
zoom.set_ylabel("Lowest frequency (THz)")
zoom.set_xlabel("Temperature (K)")
zoom.set_title("Soft-mode zoom: frequency below zero indicates an imaginary mode")
zoom.set_xticks(temperatures)

for point in ("Gamma", "Y", "Global minimum"):
    offset = 0.025 if point != "Global minimum" else -0.045
    for temperature, value in zip(temperatures, values[point]):
        zoom.annotate(f"{value:+.3f}", (temperature, value), xytext=(0, 5 if offset > 0 else -12),
                      textcoords="offset points", ha="center", fontsize=7.5,
                      color=colors.get(point, "#111827"))

fig.tight_layout()
plot_path = ROOT / "phonon_stability_vs_temperature.png"
fig.savefig(plot_path, dpi=260, bbox_inches="tight")
print(plot_path)
print(csv_path)
