#!/usr/bin/env python3
"""Plot and preliminarily fit only the currently completed Fe2O3 EOS points."""

from pathlib import Path
import csv
import re

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


POINTS = ["V092", "V094", "V096", "V098", "V100", "V102", "V104", "V106", "V108"]
RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
EV_A3_TO_GPA = 160.21766208


def birch_murnaghan(volume, e0, v0, b0, bp):
    x = (v0 / volume) ** (2.0 / 3.0)
    return e0 + 9.0 * v0 * b0 / 16.0 * (
        bp * (x - 1.0) ** 3 + (x - 1.0) ** 2 * (6.0 - 4.0 * x)
    )


def extract(path):
    if not path.exists():
        return None
    text = path.read_text(errors="replace")
    if "JOB DONE" not in text:
        return None
    energies = re.findall(r"!\s+total energy\s+=\s+([-+0-9.Ee]+)\s+Ry", text)
    volumes = re.findall(r"unit-cell volume\s+=\s+([-+0-9.Ee]+)", text, re.I)
    if not energies or not volumes:
        return None
    energy = float(energies[-1]) * RY_TO_EV
    volume = float(volumes[-1]) * BOHR_TO_ANG**3
    return volume, energy


rows = []
for point in POINTS:
    candidates = [Path(point) / "alpha_Fe2O3.eos.out"]
    if point == "V100":
        candidates.append(Path("..") / "alpha_Fe2O3.relax.out")

    result = None
    source = None
    for candidate in candidates:
        result = extract(candidate)
        if result is not None:
            source = candidate
            break

    if result is None:
        print(f"Skipping {point}: output missing or incomplete")
        continue

    volume, energy = result
    rows.append((point, str(source), volume, energy))
    print(f"Using {point}: {source}")

if len(rows) < 2:
    raise SystemExit("At least two completed points are needed to make a plot.")

rows.sort(key=lambda row: row[2])
labels = [row[0] for row in rows]
volume = np.array([row[2] for row in rows])
energy = np.array([row[3] for row in rows])
relative_energy = (energy - energy.min()) / 2.0  # primitive cell has 2 Fe2O3

with Path("completed_eos_points.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["point", "source_file", "volume_A3", "energy_eV", "relative_eV_per_Fe2O3"])
    for row, relative in zip(rows, relative_energy):
        writer.writerow([row[0], row[1], row[2], row[3], relative])

fig, ax = plt.subplots(figsize=(6.4, 4.9))
ax.scatter(volume, relative_energy, color="firebrick", s=42, zorder=3, label="completed QE points")
for label, xvalue, yvalue in zip(labels, volume, relative_energy):
    ax.annotate(label, (xvalue, yvalue), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8)

fit_message = "No BM fit: fewer than five completed points."
if len(rows) >= 5:
    guess = [float(energy.min()), float(volume[np.argmin(energy)]), 1.2, 4.0]
    lower = [energy.min() - 10.0, volume.min() * 0.85, 1.0e-6, 0.5]
    upper = [energy.min() + 10.0, volume.max() * 1.20, 10.0, 12.0]
    try:
        parameters, _ = curve_fit(
            birch_murnaghan, volume, energy, p0=guess,
            bounds=(lower, upper), maxfev=200000,
        )
        e0, v0, b0, bp = parameters
        grid = np.linspace(min(volume.min(), v0) * 0.995, max(volume.max(), v0) * 1.005, 500)
        fitted_relative = (birch_murnaghan(grid, *parameters) - energy.min()) / 2.0
        ax.plot(grid, fitted_relative, color="black", lw=1.2, label="preliminary BM fit")
        ax.axvline(v0, color="0.5", ls="--", lw=0.8)
        bracketed = volume.min() < v0 < volume.max()
        fit_message = (
            "PRELIMINARY third-order Birch-Murnaghan fit\n"
            f"Completed points: {', '.join(labels)}\n"
            f"V0 = {v0:.6f} A^3 (primitive cell)\n"
            f"B0 = {b0 * EV_A3_TO_GPA:.4f} GPa\n"
            f"B0' = {bp:.6f}\n"
            f"V0 bracketed by completed data: {bracketed}\n"
            "Do not report this as final until the remaining points finish.\n"
        )
    except Exception as exc:
        fit_message = f"Preliminary BM fit failed: {exc}"

ax.set_xlabel(r"Primitive-cell volume ($\mathrm{\AA}^3$)")
ax.set_ylabel(r"Relative energy (eV / Fe$_2$O$_3$)")
ax.set_title(r"Completed AFM $\alpha$-Fe$_2$O$_3$ EOS points")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("alpha_Fe2O3_completed_eos_points.png", dpi=220)
fig.savefig("alpha_Fe2O3_completed_eos_points.pdf")

Path("completed_points_fit.txt").write_text(fit_message + "\n")
print("\n" + fit_message)
print("Wrote completed_eos_points.csv")
print("Wrote alpha_Fe2O3_completed_eos_points.png")

