#!/usr/bin/env python3
"""Create paper-ready Fe2O3 analyses from completed EOS and V100 phonons."""

from __future__ import annotations

import csv
from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
FE2O3 = HERE.parent
EOS = FE2O3 / "eos_birch_murnaghan"
PHONONS = FE2O3 / "phonopy_V100_2x2x2"

THZ_TO_CM1 = 33.3564095198152
THZ_TO_MEV = 4.135667696
COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "gray": "#555555",
}


def load_eos() -> dict[str, np.ndarray]:
    path = EOS / "dftu_eos_diagnostics.csv"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    numeric = [
        "volume_A3",
        "energy_Ry",
        "pressure_GPa",
        "stress_anisotropy_kbar",
        "total_magnetization_muB",
        "absolute_magnetization_muB",
        "Fe_sphere_abs_moment_muB",
        "Fe_Hubbard_abs_moment_muB",
        "Hubbard_energy_eV",
        "O_internal_x",
        "Fe_internal_x",
        "Fe_O_min_A",
        "Fe_O_mean6_A",
        "Fe_O_max6_A",
    ]
    data: dict[str, np.ndarray] = {
        "volume_label": np.array([row["volume_label"] for row in rows])
    }
    for key in numeric:
        data[key] = np.array([float(row[key]) for row in rows])
    return data


def read_bm_summary() -> dict[str, float]:
    text = (EOS / "all_completed_birch_murnaghan_summary.txt").read_text()
    standard = text.split("STANDARD FIT:", 1)[1]
    values = {}
    for key in ("V0", "B0", "B0'"):
        match = re.search(rf"^{re.escape(key)}\s*=\s*([-+0-9.]+)", standard, re.MULTILINE)
        if not match:
            raise RuntimeError(f"Could not read {key} from Birch-Murnaghan summary")
        values[key] = float(match.group(1))
    return values


def linear_slope(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)


def plot_magnetostructural(data: dict[str, np.ndarray], bm: dict[str, float]) -> None:
    pressure = data["pressure_GPa"]
    volume = data["volume_A3"]
    order_v = np.argsort(volume)
    order_p = np.argsort(pressure)

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.7), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(volume[order_v], pressure[order_v], "o-", color=COLORS["blue"], lw=1.8)
    ax.axhline(0, color="0.35", lw=0.8)
    ax.axvline(bm["V0"], color=COLORS["red"], ls="--", lw=1.2,
               label=rf"BM $V_0={bm['V0']:.2f}$ Å$^3$")
    ax.set(xlabel=rf"Primitive-cell volume (Å$^3$)", ylabel="Pressure (GPa)")
    ax.legend(frameon=False, fontsize=9)
    ax.text(0.02, 0.95, "(a)", transform=ax.transAxes, va="top", fontweight="bold")

    ax = axes[0, 1]
    p = pressure[order_p]
    lower = data["Fe_O_min_A"][order_p]
    upper = data["Fe_O_max6_A"][order_p]
    mean = data["Fe_O_mean6_A"][order_p]
    ax.fill_between(p, lower, upper, color=COLORS["orange"], alpha=0.22,
                    label="Fe–O range (six neighbors)")
    ax.plot(p, mean, "o-", color=COLORS["red"], lw=1.8, label="Mean Fe–O")
    ax.set(xlabel="Pressure (GPa)", ylabel=rf"Fe–O distance (Å)")
    ax.legend(frameon=False, fontsize=9)
    ax.text(0.02, 0.95, "(b)", transform=ax.transAxes, va="top", fontweight="bold")

    ax = axes[1, 0]
    ax.plot(p, data["Fe_Hubbard_abs_moment_muB"][order_p], "o-",
            color=COLORS["blue"], lw=1.8, label="Hubbard-site moment")
    ax.plot(p, data["Fe_sphere_abs_moment_muB"][order_p], "s--",
            color=COLORS["green"], lw=1.5, label="Fe-sphere moment")
    ax.set(xlabel="Pressure (GPa)", ylabel=rf"Absolute Fe moment ($\mu_\mathrm{{B}}$)")
    ax.legend(frameon=False, fontsize=9)
    ax.text(0.02, 0.95, "(c)", transform=ax.transAxes, va="top", fontweight="bold")

    ax = axes[1, 1]
    ax.plot(p, data["Hubbard_energy_eV"][order_p], "o-", color=COLORS["purple"], lw=1.8)
    ax.set(xlabel="Pressure (GPa)", ylabel="Hubbard energy (eV/cell)")
    ax.text(0.02, 0.95, "(d)", transform=ax.transAxes, va="top", fontweight="bold")

    for ax in axes.flat:
        ax.tick_params(direction="in", top=True, right=True)
        ax.grid(alpha=0.16, lw=0.6)

    fig.suptitle(r"AFM $\alpha$-Fe$_2$O$_3$: coupled compression, bonding, and magnetism", fontsize=14)
    for suffix in ("png", "pdf"):
        fig.savefig(HERE / f"Fe2O3_eos_magnetostructural.{suffix}", dpi=350, bbox_inches="tight")
    plt.close(fig)


def load_projected_dos() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    raw = np.loadtxt(PHONONS / "projected_dos.dat")
    if raw.shape[1] != 11:
        raise RuntimeError(f"Expected frequency plus ten atomic PDOS columns, got {raw.shape[1]}")
    frequency = raw[:, 0]
    oxygen = raw[:, 1:7].sum(axis=1)
    iron = raw[:, 7:11].sum(axis=1)
    return frequency, oxygen, iron, oxygen + iron


def plot_projected_dos(
    frequency: np.ndarray, oxygen: np.ndarray, iron: np.ndarray, total: np.ndarray
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.fill_between(frequency, oxygen, color=COLORS["red"], alpha=0.42, label="O projection")
    ax.fill_between(frequency, iron, color=COLORS["blue"], alpha=0.48, label="Fe projection")
    ax.plot(frequency, total, color="black", lw=1.3, label="Total")
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_xlim(-0.15, 20.6)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Frequency (THz)")
    ax.set_ylabel("Phonon DOS (states THz$^{-1}$ cell$^{-1}$)")
    ax.set_title(r"AFM $\alpha$-Fe$_2$O$_3$ V100: species-projected phonon DOS")
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.tick_params(direction="in", top=True, right=True)
    ax.grid(axis="x", alpha=0.16, lw=0.6)
    for suffix in ("png", "pdf"):
        fig.savefig(HERE / f"Fe2O3_V100_projected_phonon_dos.{suffix}", dpi=350, bbox_inches="tight")
    plt.close(fig)


def gamma_modes() -> list[dict[str, float | int | str]]:
    with h5py.File(PHONONS / "band.hdf5", "r") as h5:
        frequencies = np.array(h5["frequency"][0, 0, :], dtype=float)
        eigenvectors = np.array(h5["eigenvector"][0, 0, :, :])

    if eigenvectors.shape != (30, 30):
        raise RuntimeError(f"Unexpected Gamma eigenvector shape: {eigenvectors.shape}")
    norms = np.sum(np.abs(eigenvectors) ** 2, axis=0)
    if not np.allclose(norms, 1.0, atol=1e-7):
        raise RuntimeError("Gamma eigenvectors are not normalized by mode columns")

    groups: list[list[int]] = []
    for index, frequency in enumerate(frequencies):
        if groups and abs(frequency - frequencies[groups[-1][0]]) < 1e-4:
            groups[-1].append(index)
        else:
            groups.append([index])
    degeneracy = {index: len(group) for group in groups for index in group}

    rows: list[dict[str, float | int | str]] = []
    for index, frequency in enumerate(frequencies):
        vector = eigenvectors[:, index]
        atom_weight = (np.abs(vector) ** 2).reshape(10, 3).sum(axis=1)
        oxygen = float(atom_weight[:6].sum())
        iron = float(atom_weight[6:].sum())
        if abs(frequency) < 0.05:
            character = "acoustic translation"
        elif oxygen >= 0.70:
            character = "O-dominated"
        elif iron >= 0.70:
            character = "Fe-dominated"
        else:
            character = "mixed Fe–O"
        rows.append(
            {
                "mode": index + 1,
                "frequency_THz": float(frequency),
                "frequency_cm-1": float(frequency * THZ_TO_CM1),
                "energy_meV": float(frequency * THZ_TO_MEV),
                "degeneracy": degeneracy[index],
                "O_weight_percent": 100 * oxygen,
                "Fe_weight_percent": 100 * iron,
                "species_character": character,
            }
        )
    return rows


def write_gamma_table(rows: list[dict[str, float | int | str]]) -> None:
    path = HERE / "Fe2O3_V100_Gamma_modes.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            for key in ("frequency_THz", "frequency_cm-1", "energy_meV"):
                formatted[key] = f"{float(row[key]):.8f}"
            for key in ("O_weight_percent", "Fe_weight_percent"):
                formatted[key] = f"{float(row[key]):.4f}"
            writer.writerow(formatted)


def write_summary(
    data: dict[str, np.ndarray],
    bm: dict[str, float],
    frequency: np.ndarray,
    oxygen: np.ndarray,
    iron: np.ndarray,
    total: np.ndarray,
    modes: list[dict[str, float | int | str]],
) -> None:
    pressure = data["pressure_GPa"]
    bond_slope, _ = linear_slope(pressure, data["Fe_O_mean6_A"])
    hubbard_moment_slope, _ = linear_slope(pressure, data["Fe_Hubbard_abs_moment_muB"])
    sphere_moment_slope, _ = linear_slope(pressure, data["Fe_sphere_abs_moment_muB"])
    hubbard_energy_slope, _ = linear_slope(pressure, data["Hubbard_energy_eV"])
    mask = frequency >= 0
    int_o = float(np.trapz(oxygen[mask], frequency[mask]))
    int_fe = float(np.trapz(iron[mask], frequency[mask]))
    int_total = float(np.trapz(total[mask], frequency[mask]))
    optical = [float(row["frequency_THz"]) for row in modes if abs(float(row["frequency_THz"])) >= 0.05]

    text = f"""Existing-results analysis for AFM alpha-Fe2O3
================================================

No new DFT calculations were used.

EOS and magnetostructural data
------------------------------
Pressure range: {pressure.min():.3f} to {pressure.max():.3f} GPa
Volume range: {data['volume_A3'].min():.6f} to {data['volume_A3'].max():.6f} A^3
Birch-Murnaghan V0: {bm['V0']:.6f} A^3
Birch-Murnaghan B0: {bm['B0']:.6f} GPa
Birch-Murnaghan B0': {bm["B0'"]:.6f}

Linear descriptive slopes over the sampled pressure interval:
Mean Fe-O distance: {bond_slope:+.6f} A/GPa
Fe Hubbard-site moment: {hubbard_moment_slope:+.6f} muB/GPa
Fe sphere moment: {sphere_moment_slope:+.6f} muB/GPa
Hubbard energy: {hubbard_energy_slope:+.6f} eV/GPa

These slopes summarize the sampled range and are not an EOS model.

V100 lattice dynamics
----------------------
Integrated projected DOS (non-negative frequency grid):
O contribution: {int_o:.5f} states/cell
Fe contribution: {int_fe:.5f} states/cell
Total: {int_total:.5f} states/cell (expected approximately 30)

Gamma acoustic residual: {min(float(row['frequency_THz']) for row in modes[:3]):+.8f} THz
Lowest optical Gamma frequency: {min(optical):.8f} THz
Highest Gamma frequency: {max(float(row['frequency_THz']) for row in modes):.8f} THz

The three tiny negative Gamma frequencies are numerical acoustic residuals.
Species weights are squared, normalized mass-weighted dynamical-matrix
eigenvector components. Crystallographic irreducible-representation and
Raman/IR activity labels are not assigned here.
"""
    (HERE / "existing_results_summary.txt").write_text(text)


def main() -> None:
    HERE.mkdir(exist_ok=True)
    data = load_eos()
    bm = read_bm_summary()
    plot_magnetostructural(data, bm)
    frequency, oxygen, iron, total = load_projected_dos()
    plot_projected_dos(frequency, oxygen, iron, total)
    modes = gamma_modes()
    write_gamma_table(modes)
    write_summary(data, bm, frequency, oxygen, iron, total, modes)
    print(f"Wrote publication analysis to {HERE}")


if __name__ == "__main__":
    main()
