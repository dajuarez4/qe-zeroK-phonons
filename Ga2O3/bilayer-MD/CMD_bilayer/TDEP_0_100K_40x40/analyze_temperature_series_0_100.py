#!/usr/bin/env python3
"""Summarize and plot the 0--100 K Ga2O3 bilayer phonon series."""

from pathlib import Path
import csv
import re
import h5py
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
TEMPERATURES = list(range(0, 101, 10))


def numeric_table(path):
    rows = []
    for line in path.read_text().splitlines():
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            pass
    return np.asarray(rows)


def production_average(log_path):
    blocks, current, capture = [], [], False
    for line in log_path.read_text(errors="ignore").splitlines():
        if re.match(r"^\s*Step\s+Time\s+Temp", line):
            if current:
                blocks.append(current)
            current, capture = [], True
            continue
        if capture:
            fields = line.split()
            if len(fields) >= 10:
                try:
                    current.append([float(value) for value in fields[:10]])
                except ValueError:
                    pass
    if current:
        blocks.append(current)
    block = np.asarray(max(blocks, key=len))
    return block[:, 2].mean(), block[:, 7:10].mean(axis=0) * 1e-4


def fit_r2(path):
    text = path.read_text(errors="ignore")
    match = re.search(r"second order:\s+[\d.Ee+-]+\s+[\d.Ee+-]+\s+[\d.Ee+-]+\s+([\d.Ee+-]+)", text)
    return float(match.group(1)) if match else np.nan


records = []
dispersions = {}
for nominal in TEMPERATURES:
    if nominal == 0:
        table = numeric_table(ROOT / "phonopy_0K_tabgap_40x40_amp002/dispersion_0K_tabgap.dat")
        measured, stress, r2, z_fraction = 0.0, np.full(3, np.nan), np.nan, np.nan
    else:
        folder = ROOT / ("tdep_10K_40x40_extended" if nominal == 10 else f"tdep_{nominal}K_40x40")
        table = numeric_table(folder / "outfile.dispersion_relations")
        if nominal == 10:
            log = ROOT / "tdep_10K_40x40/log_extended.lammps"
            fitlog = folder / "extract_forceconstants_5frames_25ps.log"
        else:
            log = folder / "log.lammps"
            fitlog = folder / "extract_forceconstants.log"
        measured, stress = production_average(log)
        r2 = fit_r2(fitlog)
        with h5py.File(folder / "outfile.dispersion_relations.hdf5", "r") as h5:
            frequencies = h5["frequencies"][:]
            iq, im = np.unravel_index(np.argmin(frequencies), frequencies.shape)
            vector = h5["eigenvectors_re"][iq, im] + 1j * h5["eigenvectors_im"][iq, im]
            weights = np.abs(vector.reshape(10, 3)) ** 2
            z_fraction = float(weights[:, 2].sum() / weights.sum())
    modes = table[:, 1:]
    last = len(modes) - 1
    gamma = float(modes[0].min())
    xpoint = float(modes[round(0.25 * last)].min())
    spoint = float(modes[round(0.50 * last)].min())
    ypoint = float(modes[round(0.75 * last)].min())
    global_min = float(modes.min())
    negative_qrows = int(np.any(modes < 0, axis=1).sum())
    dispersions[nominal] = modes
    records.append(dict(nominal=nominal, measured=measured, gamma=gamma, X=xpoint,
                        S=spoint, Y=ypoint, global_min=global_min,
                        negative_qrows=negative_qrows, r2=r2, z_fraction=z_fraction,
                        pxx=stress[0], pyy=stress[1], pzz=stress[2]))

csv_path = ROOT / "phonon_temperature_report_0_100K.csv"
fields = list(records[0])
with csv_path.open("w", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=fields)
    writer.writeheader()
    writer.writerows(records)

# High-symmetry stability plot.
fig, (full, zoom) = plt.subplots(2, 1, figsize=(9.5, 8.5), sharex=True,
                                gridspec_kw={"height_ratios": [1.25, 1]})
styles = {"gamma": ("Gamma", "#7c3aed", "o"), "X": ("X", "#2563eb", "s"),
          "S": ("S", "#16a34a", "^"), "Y": ("Y", "#dc2626", "D")}
temps = np.array([row["nominal"] for row in records])
for axis in (full, zoom):
    axis.axhspan(-1, 0, color="#fee2e2", alpha=0.60)
    axis.axhline(0, color="black", lw=1)
    for key, (label, color, marker) in styles.items():
        axis.plot(temps, [row[key] for row in records], marker=marker, color=color,
                  lw=1.8, ms=5.5, label=label if axis is full else None)
    axis.plot(temps, [row["global_min"] for row in records], "x--", color="#111827",
              lw=1.5, ms=6, label="Global minimum" if axis is full else None)
    axis.grid(color="0.90")
full.set_ylim(-0.9, 2.6)
full.set_ylabel("Lowest frequency (THz)")
full.set_title("Ga$_2$O$_3$ bilayer stability at high-symmetry points, 0--100 K\n"
               "40x40x1 cell; finite displacement at 0 K and TDEP cutoff 5.5 Angstrom")
full.legend(ncol=3, fontsize=9)
zoom.set_ylim(-0.8, 0.8)
zoom.set_ylabel("Lowest frequency (THz)")
zoom.set_xlabel("Temperature (K)")
zoom.set_title("Soft-mode zoom; negative values are imaginary TDEP frequencies")
zoom.set_xticks(temps)
fig.tight_layout()
fig.savefig(ROOT / "phonon_stability_vs_temperature_0_100K.png", dpi=250, bbox_inches="tight")

# TDEP quality and residual-stress diagnostics.
finite = records[1:]
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
t = np.array([row["nominal"] for row in finite])
axes[0].plot(t, [row["r2"] for row in finite], "o-", color="#2563eb")
axes[0].set_ylabel("TDEP force-fit R-squared")
axes[0].set_ylim(0, 1)
axes[0].set_title("Harmonic fit quality")
axes[1].plot(t, [row["pxx"] for row in finite], "o-", label="Pxx")
axes[1].plot(t, [row["pyy"] for row in finite], "s-", label="Pyy")
axes[1].plot(t, [0.5 * (row["pxx"] + row["pyy"]) for row in finite], "D--", color="#111827", label="in-plane mean")
axes[1].axhline(0, color="black", lw=0.8)
axes[1].set_ylabel("Mean pressure (GPa)")
axes[1].set_title("Fixed-cell residual stress")
axes[1].legend(fontsize=8)
axes[2].plot(t, 100 * np.array([row["z_fraction"] for row in finite]), "o-", color="#7c3aed")
axes[2].set_ylabel("Out-of-plane weight (%)")
axes[2].set_ylim(0, 105)
axes[2].set_title("Global soft-mode polarization")
for axis in axes:
    axis.set_xlabel("Temperature (K)")
    axis.grid(color="0.91")
fig.tight_layout()
fig.savefig(ROOT / "tdep_quality_stress_polarization_10_100K.png", dpi=250, bbox_inches="tight")

# Full dispersion overlay. Low alpha keeps the 11-data-set plot readable.
cmap = plt.get_cmap("turbo")
fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                         gridspec_kw={"height_ratios": [2.1, 1]})
for axis in axes:
    axis.axhline(0, color="black", lw=0.9)
    for tick in range(5):
        axis.axvline(tick, color="0.87", lw=0.7)
    for index, nominal in enumerate(TEMPERATURES):
        modes = dispersions[nominal]
        xaxis = np.linspace(0, 4, len(modes))
        color = "#111827" if nominal == 0 else cmap(nominal / 100)
        for branch in range(modes.shape[1]):
            axis.plot(xaxis, modes[:, branch], color=color, lw=0.35, alpha=0.35)
        axis.plot([], [], color=color, lw=1.5, label=f"{nominal} K")
    axis.grid(axis="y", color="0.94")
axes[0].set_ylabel("Frequency (THz)")
axes[0].set_title("Ga$_2$O$_3$ bilayer phonons, 0--100 K")
axes[0].legend(ncol=6, fontsize=7.5)
axes[1].set_ylim(-0.9, 3.5)
axes[1].fill_between([0, 4], -0.9, 0, color="#fee2e2", alpha=0.45)
axes[1].set_ylabel("Frequency (THz)")
axes[1].set_xlabel("Wave-vector path")
axes[1].set_title("Low-frequency zoom")
axes[1].set_xticks(range(5), [r"$\Gamma$", "X", "S", "Y", r"$\Gamma$"])
fig.tight_layout()
fig.savefig(ROOT / "phonon_temperature_comparison_0K_to_100K_40x40.png", dpi=250, bbox_inches="tight")

# Markdown report.
report = [
    "# Ga2O3 bilayer: MD/TDEP phonons from 0 to 100 K",
    "",
    "## Protocol",
    "",
    "- 40x40x1 bilayer supercell (16,000 atoms).",
    "- tabGAP Ga-O potential; 1 fs timestep.",
    "- 5 ps NVT equilibration and 25 ps NVT production at each finite temperature.",
    "- 26 stored configurations per trajectory; five frames at 0, 6, 12, 19, and 25 ps used per TDEP fit because of the WSL memory ceiling.",
    "- Second-order TDEP cutoff 5.5 Angstrom; Gamma-X-S-Y-Gamma path.",
    "",
    "## Results",
    "",
    "| T (K) | measured T (K) | Gamma | X | S | Y | global min | R2 | z weight | mean in-plane P (GPa) |",
    "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
]
for row in records:
    pin = np.nan if row["nominal"] == 0 else 0.5 * (row["pxx"] + row["pyy"])
    report.append(f'| {row["nominal"]} | {row["measured"]:.3f} | {row["gamma"]:+.3f} | {row["X"]:+.3f} | {row["S"]:+.3f} | {row["Y"]:+.3f} | {row["global_min"]:+.3f} | {row["r2"]:.3f} | {100*row["z_fraction"]:.1f}% | {pin:+.3f} |')
report += [
    "",
    "Frequencies are in THz. A negative TDEP frequency denotes an imaginary harmonic mode of the fitted effective force constants; it does not by itself establish a static structural instability.",
    "",
    "## Static-mode verification",
    "",
    "For the 50 K Gamma eigenvector, 17 tabGAP single-point structures spanning Q = -0.20 to +0.20 Angstrom were evaluated. The local quadratic coefficient is +6.371 eV/Angstrom^2 per 10-atom cell and the fitted minimum is Q = -0.0027 Angstrom. The positive curvature and absence of a double well demonstrate static stability along the Gamma coordinate despite the negative TDEP frequency.",
    "",
    "## Interpretation",
    "",
    "The low branches are predominantly out of plane and therefore flexural. Their fitted frequencies are affected by finite sampling, the decreasing quality of a purely second-order model, residual anisotropic tensile stress from the fixed in-plane cell, and the known limited phonon accuracy of the general-purpose tabGAP near the lowest acoustic branch. These frequencies should not be reported as proof of a phase instability without direct DFT validation or a finite-temperature cell/stress treatment.",
]
(ROOT / "REPORT_MD_TDEP_0_100K.md").write_text("\n".join(report) + "\n")

print(csv_path)
print(ROOT / "REPORT_MD_TDEP_0_100K.md")
