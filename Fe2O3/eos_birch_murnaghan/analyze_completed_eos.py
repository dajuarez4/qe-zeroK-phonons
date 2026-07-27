#!/usr/bin/env python3
"""Analyze all unique completed QE EOS outputs, including every relaxation energy."""

from pathlib import Path
import csv
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
EV_A3_TO_GPA = 160.21766208
VREF = 100.624948
AREF = 5.0355
CREF = 13.7471


def bm3(volume, e0, v0, b0, bp):
    x = (v0 / volume) ** (2.0 / 3.0)
    return e0 + 9.0 * v0 * b0 / 16.0 * (
        bp * (x - 1.0) ** 3 + (x - 1.0) ** 2 * (6.0 - 4.0 * x)
    )


def floats(pattern, text):
    return [float(value) for value in re.findall(pattern, text, flags=re.I)]


def fit(volume, energy):
    guess = [float(energy.min()), 1.055 * VREF, 1.2, 4.0]
    lower = [energy.min() - 10.0, volume.min() * 0.95, 0.05, 0.5]
    upper = [energy.min() + 10.0, 1.20 * VREF, 5.0, 12.0]
    parameters, _ = curve_fit(
        bm3, volume, energy, p0=guess, bounds=(lower, upper), maxfev=200000
    )
    fitted = bm3(volume, *parameters)
    rmse = 1000.0 * np.sqrt(np.mean((energy - fitted) ** 2))
    return parameters, rmse


raw_records = []
for path in sorted(Path(".").glob("*.out")):
    text = path.read_text(errors="replace")
    if "JOB DONE" not in text:
        continue
    volumes = floats(r"unit-cell volume\s+=\s+([-+0-9.Ee]+)", text)
    energies = floats(r"!\s+total energy\s+=\s+([-+0-9.Ee]+)\s+Ry", text)
    forces = floats(r"Total force\s+=\s+([-+0-9.Ee]+)", text)
    total_mag = floats(r"total magnetization\s+=\s+([-+0-9.Ee]+)", text)
    absolute_mag = floats(r"absolute magnetization\s+=\s+([-+0-9.Ee]+)", text)
    if not volumes or not energies:
        continue
    volume = volumes[-1] * BOHR_TO_ANG**3
    ratio = volume / VREF
    label = f"V{int(round(100.0 * ratio)):03d}"
    raw_records.append(
        {
            "label": label,
            "path": path,
            "volume": volume,
            "energies": np.array(energies) * RY_TO_EV,
            "force": forces[-1] if forces else np.nan,
            "total_mag": total_mag[-1] if total_mag else np.nan,
            "absolute_mag": absolute_mag[-1] if absolute_mag else np.nan,
        }
    )

raw_records.sort(key=lambda record: record["volume"])
unique = {}
duplicates = []
for record in raw_records:
    key = round(record["volume"], 5)
    previous = unique.get(key)
    if previous is None:
        unique[key] = record
    elif record["force"] < previous["force"]:
        duplicates.append(previous["path"])
        unique[key] = record
    else:
        duplicates.append(record["path"])
records = sorted(unique.values(), key=lambda record: record["volume"])
if len(records) < 5:
    raise SystemExit(f"Need at least five unique completed volumes; found {len(records)}")

final_volume = np.array([record["volume"] for record in records])
final_energy = np.array([record["energies"][-1] for record in records])
all_volume = np.concatenate(
    [np.full(len(record["energies"]), record["volume"]) for record in records]
)
all_energy = np.concatenate([record["energies"] for record in records])

final_parameters, final_rmse = fit(final_volume, final_energy)
all_parameters, all_rmse = fit(all_volume, all_energy)


def derived(parameters):
    e0, v0, b0, bp = parameters
    scale = (v0 / VREF) ** (1.0 / 3.0)
    return {
        "E0": e0,
        "V0": v0,
        "B0_GPa": b0 * EV_A3_TO_GPA,
        "B0p": bp,
        "a_hex": AREF * scale,
        "c_hex": CREF * scale,
    }


final_result = derived(final_parameters)
all_result = derived(all_parameters)

with Path("all_completed_relaxation_steps.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(
        ["volume_label", "source_file", "ionic_step", "volume_A3", "energy_Ry", "energy_eV", "is_final"]
    )
    for record in records:
        for step, energy in enumerate(record["energies"]):
            writer.writerow(
                [
                    record["label"], record["path"], step, record["volume"],
                    energy / RY_TO_EV, energy, step == len(record["energies"]) - 1,
                ]
            )

lines = [
    f"{len(records)} completed AFM alpha-Fe2O3 EOS calculations",
    "",
    "Completed-point checks:",
]
for record in records:
    lines.append(
        f"{record['label']}: V={record['volume']:.6f} A^3, "
        f"steps={len(record['energies'])}, Efinal={record['energies'][-1] / RY_TO_EV:.8f} Ry, "
        f"force={record['force']:.6f} Ry/bohr, Mtot={record['total_mag']:.2f}, "
        f"|M|={record['absolute_mag']:.2f}"
    )

lines.extend(
    [
        "",
        f"Excluded duplicate-volume outputs: {', '.join(str(path) for path in duplicates) if duplicates else 'none'}",
        "",
        f"REQUESTED FIT: every reported relaxation energy ({len(all_energy)} observations)",
        f"V0       = {all_result['V0']:.6f} A^3",
        f"a_hex    = {all_result['a_hex']:.6f} A",
        f"c_hex    = {all_result['c_hex']:.6f} A",
        f"B0       = {all_result['B0_GPa']:.6f} GPa",
        f"B0'      = {all_result['B0p']:.6f}",
        f"RMSE     = {all_rmse:.6f} meV/primitive cell",
        "Warning: intermediate ionic configurations are not equilibrium EOS points.",
        "Repeated volumes are weighted according to their number of reported relaxation energies.",
        "",
        f"STANDARD FIT: final relaxed energy at each volume ({len(records)} observations)",
        f"V0       = {final_result['V0']:.6f} A^3",
        f"a_hex    = {final_result['a_hex']:.6f} A",
        f"c_hex    = {final_result['c_hex']:.6f} A",
        f"B0       = {final_result['B0_GPa']:.6f} GPa",
        f"B0'      = {final_result['B0p']:.6f}",
        f"RMSE     = {final_rmse:.6f} meV/primitive cell",
        f"Calculated energy minimum bracketed: {0 < int(np.argmin(final_energy)) < len(final_energy)-1}",
        f"Fitted V0 inside sampled range: {final_volume.min() < final_result['V0'] < final_volume.max()}",
    ]
)
summary = "\n".join(lines) + "\n"
Path("all_completed_birch_murnaghan_summary.txt").write_text(summary)
print(summary)

energy_reference = final_energy.min()
grid = np.linspace(
    min(final_volume.min(), all_result["V0"], final_result["V0"]) * 0.995,
    max(final_volume.max(), all_result["V0"], final_result["V0"]) * 1.005,
    600,
)

fig, (ax_eos, ax_relax) = plt.subplots(2, 1, figsize=(7.2, 8.2), gridspec_kw={"height_ratios": [1.35, 1]})

ax_eos.scatter(
    all_volume, (all_energy - energy_reference) / 2.0,
    s=22, facecolors="none", edgecolors="tab:blue", alpha=0.65,
    label=f"all {len(all_energy)} relaxation energies",
)
ax_eos.scatter(
    final_volume, (final_energy - energy_reference) / 2.0,
    s=48, color="firebrick", zorder=4, label=f"{len(records)} final relaxed energies",
)
ax_eos.plot(
    grid, (bm3(grid, *all_parameters) - energy_reference) / 2.0,
    color="tab:blue", ls="--", lw=1.2, label="fit using all relaxation energies",
)
ax_eos.plot(
    grid, (bm3(grid, *final_parameters) - energy_reference) / 2.0,
    color="black", lw=1.3, label="fit using final energies",
)
for record, xvalue, yvalue in zip(records, final_volume, (final_energy - energy_reference) / 2.0):
    ax_eos.annotate(record["label"], (xvalue, yvalue), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8)
if max(all_result["V0"], final_result["V0"]) > final_volume.max():
    ax_eos.axvspan(final_volume.max(), grid.max(), color="0.92", zorder=-5, label="extrapolated region")
elif min(all_result["V0"], final_result["V0"]) < final_volume.min():
    ax_eos.axvspan(grid.min(), final_volume.min(), color="0.92", zorder=-5, label="extrapolated region")
ax_eos.set_ylabel(r"Relative energy (eV / Fe$_2$O$_3$)")
ax_eos.set_title(rf"{len(records)} completed AFM $\alpha$-Fe$_2$O$_3$ EOS points")
ax_eos.legend(frameon=False, fontsize=8)

for record in records:
    delta = (record["energies"] - record["energies"][-1]) * 1000.0 / 2.0
    ax_relax.plot(np.arange(len(delta)), delta, marker="o", ms=3, lw=1.0, label=record["label"])
ax_relax.set_xlabel("Reported relaxation energy index")
ax_relax.set_ylabel(r"$E-E_{final}$ (meV / Fe$_2$O$_3$)")
ax_relax.set_yscale("symlog", linthresh=0.01)
ax_relax.set_title("Internal-coordinate relaxation histories")
ax_relax.legend(frameon=False, ncol=5, fontsize=8)

fig.tight_layout()
fig.savefig("alpha_Fe2O3_all_completed_points_bm.png", dpi=230)
fig.savefig("alpha_Fe2O3_all_completed_points_bm.pdf")
