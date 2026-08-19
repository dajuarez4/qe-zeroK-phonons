#!/usr/bin/env python3
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
q = np.loadtxt(HERE / "amplitudes.dat")
energy = np.array([float(line.split()[1]) for line in (HERE / "energies.raw").read_text().splitlines()])
relative_mev = (energy - energy.min()) * 1000.0

# Fit the local |Q| <= 0.1 A region to avoid higher-order terms dominating.
mask = np.abs(q) <= 0.1000001
quadratic = np.polyfit(q[mask], energy[mask], 2)
vertex = -quadratic[1] / (2.0 * quadratic[0])


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
            elif current and "Loop time" in line:
                blocks.append(current)
                current, capture = [], False
    if current:
        blocks.append(current)
    block = np.asarray(max(blocks, key=len))
    # LAMMPS metal pressure is bar; 1 bar = 1e-4 GPa.
    return block[:, 2].mean(), block[:, 7:10].mean(axis=0) * 1e-4


temperatures, stresses = [], []
for nominal in (10, 20, 30, 40, 50):
    folder = ROOT / ("tdep_10K_40x40" if nominal == 10 else f"tdep_{nominal}K_40x40")
    log = folder / ("log_extended.lammps" if nominal == 10 else "log.lammps")
    measured, stress = production_average(log)
    temperatures.append(measured)
    stresses.append(stress)
stresses = np.asarray(stresses)

np.savetxt(HERE / "gamma_energy_scan.csv", np.column_stack((q, energy, relative_mev)),
           delimiter=",", header="Q_A,total_energy_eV,relative_energy_meV", comments="")
np.savetxt(HERE / "residual_stress_vs_temperature.csv",
           np.column_stack((temperatures, stresses)), delimiter=",",
           header="temperature_K,Pxx_GPa,Pyy_GPa,Pzz_GPa", comments="")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].plot(q, relative_mev, "o", color="#2563eb", label="tabGAP single points")
fine = np.linspace(q.min(), q.max(), 400)
fit_relative = (np.polyval(quadratic, fine) - energy.min()) * 1000.0
axes[0].plot(fine, fit_relative, color="#dc2626", lw=1.8,
             label=f"local quadratic fit; vertex {vertex:+.4f} A")
axes[0].axvline(0, color="black", lw=0.8)
axes[0].set_xlabel("Mode amplitude Q (A; maximum atomic displacement)")
axes[0].set_ylabel("Relative energy (meV per 10-atom cell)")
axes[0].set_title("50 K TDEP Gamma soft eigenvector\nevaluated on the static tabGAP surface")
axes[0].legend(fontsize=8)
axes[0].grid(color="0.92")

axes[1].axhline(0, color="black", lw=0.8)
axes[1].plot(temperatures, stresses[:, 0], "o-", label="Pxx")
axes[1].plot(temperatures, stresses[:, 1], "s-", label="Pyy")
axes[1].plot(temperatures, stresses[:, 2], "^-", label="Pzz")
axes[1].plot(temperatures, stresses[:, :2].mean(axis=1), "D--", color="#111827", label="in-plane mean")
axes[1].set_xlabel("Measured temperature (K)")
axes[1].set_ylabel("Mean pressure component (GPa)")
axes[1].set_title("Residual fixed-cell stress during production")
axes[1].legend(fontsize=8)
axes[1].grid(color="0.92")

fig.tight_layout()
out = ROOT / "soft_mode_diagnostics_energy_and_stress.png"
fig.savefig(out, dpi=250, bbox_inches="tight")
print(f"quadratic coefficient = {quadratic[0]:.8f} eV/A^2")
print(f"quadratic vertex = {vertex:.8f} A")
print(out)
