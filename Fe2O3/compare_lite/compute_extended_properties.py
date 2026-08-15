#!/usr/bin/env python3
"""Extended three-volume phonon and preliminary QHA analysis for alpha-Fe2O3.

The QHA interpolation has exactly three phonon volumes and all calculated
zero-pressure minima lie above that interval.  Outputs therefore carry an
explicit extrapolation flag and should be regarded as convergence guidance,
not production QHA values.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml
from phonopy import load
from scipy.optimize import curve_fit, minimize_scalar


HERE = Path(__file__).resolve().parent
FE = HERE.parent
VOLUMES = np.array([100.62494687294131, 102.63744343944475, 104.6499400059482])
LABELS = ["V100", "V102", "V104"]
COLORS = ["#d62728", "#1f77b4", "#ff7f0e"]
EV_PER_KJMOL = 1.0 / 96.4853321233
GPA_A3_TO_EV = 1.0 / 160.21766208
R = 8.31446261815324
NA = 6.02214076e23
H_THz_KJMOL = 0.399031271
CM1_TO_THZ = 0.0299792458
BOHR_TO_ANGSTROM = 0.529177210903
VPH_MIN, VPH_MAX = VOLUMES[0], VOLUMES[-1]


def write_csv(name, fields, rows):
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def bm3(v, e0, v0, b0_ev_a3, bp):
    eta = (v0 / np.asarray(v)) ** (2.0 / 3.0)
    return e0 + 9.0 * v0 * b0_ev_a3 / 16.0 * (
        (eta - 1.0) ** 3 * bp + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )


def load_eos():
    rows = []
    with (FE / "eos_birch_murnaghan/all_completed_relaxation_steps.csv").open() as handle:
        for row in csv.DictReader(handle):
            if row["is_final"] == "True":
                rows.append((float(row["volume_A3"]), float(row["energy_eV"])))
    arr = np.array(rows)
    arr[:, 1] -= arr[:, 1].min()
    p0 = [0.0, arr[np.argmin(arr[:, 1]), 0], 190 / 160.21766208, 4.0]
    popt, pcov = curve_fit(bm3, arr[:, 0], arr[:, 1], p0=p0, maxfev=20000)
    return arr, popt, np.sqrt(np.diag(pcov))


def load_thermal(label):
    with (FE / f"phonopy_{label}_2x2x2/thermal_properties.yaml").open() as handle:
        d = yaml.safe_load(handle)
    rows = d["thermal_properties"]
    return {
        "T": np.array([r["temperature"] for r in rows], float),
        "F": np.array([r["free_energy"] for r in rows], float),
        "S": np.array([r["entropy"] for r in rows], float),
        "Cv": np.array([r["heat_capacity"] for r in rows], float),
        "E": np.array([r["energy"] for r in rows], float),
        "ZPE": float(d["zero_point_energy"]),
    }


def load_mesh(label):
    # Use a full common mesh because independently symmetry-reduced YAML meshes
    # can contain different representatives after small volume relaxations.
    folder = FE / f"phonopy_{label}_2x2x2"
    ph = load(str(folder / "phonopy_disp.yaml"), force_constants_filename=str(folder / "FORCE_CONSTANTS"))
    ph.run_mesh([12, 12, 12], is_mesh_symmetry=False)
    mesh = ph.get_mesh_dict()
    return np.asarray(mesh["qpoints"]), np.asarray(mesh["weights"]), np.asarray(mesh["frequencies"])


def mode_gruneisen():
    meshes = [load_mesh(label) for label in LABELS]
    for q, _, _ in meshes[1:]:
        if not np.allclose(q, meshes[0][0], atol=1e-8):
            raise RuntimeError("Volume meshes do not share q points")
    weights = meshes[1][1]
    f100, f102, f104 = meshes[0][2], meshes[1][2], meshes[2][2]
    good = (f100 > 0.05) & (f102 > 0.05) & (f104 > 0.05)
    gamma = np.full_like(f102, np.nan)
    gamma[good] = -np.log(f104[good] / f100[good]) / np.log(VOLUMES[2] / VOLUMES[0])
    weighted_modes = np.repeat(weights[:, None], f102.shape[1], axis=1)

    rows = []
    for temperature in [50, 100, 200, 300, 600, 1000]:
        x = 47.9924307 * f102 / temperature
        ex = np.exp(np.clip(x, None, 700))
        cv = R * x**2 * ex / np.expm1(x) ** 2
        select = good & np.isfinite(cv)
        wg = weighted_modes[select] * cv[select]
        rows.append({
            "temperature_K": temperature,
            "Cv_weighted_gamma": np.sum(wg * gamma[select]) / np.sum(wg),
            "valid_mode_fraction_percent": 100 * np.sum(weighted_modes[good]) / np.sum(weighted_modes),
        })
    write_csv("gruneisen_temperature.csv", list(rows[0]), rows)

    finite = gamma[np.isfinite(gamma)]
    summary = [{
        "quantity": "all_valid_modes",
        "mean": np.nanmean(finite),
        "median": np.nanmedian(finite),
        "p05": np.nanpercentile(finite, 5),
        "p95": np.nanpercentile(finite, 95),
        "minimum": np.nanmin(finite),
        "maximum": np.nanmax(finite),
        "note": "Branch-index matching; crossings can create outliers",
    }]
    write_csv("gruneisen_summary.csv", list(summary[0]), summary)

    raman = np.genfromtxt(HERE / "raman_modes.csv", delimiter=",", names=True, dtype=None, encoding="utf8")
    rrows = []
    for row in raman:
        g = -np.log(row["V104_cm1"] / row["V100_cm1"]) / np.log(VOLUMES[2] / VOLUMES[0])
        rrows.append({"mode": row["mode"], "gamma_V100_V104": g,
                      "V100_cm-1": row["V100_cm1"], "V102_cm-1": row["V102_cm1"],
                      "V104_cm-1": row["V104_cm1"]})
    write_csv("raman_mode_gruneisen.csv", list(rrows[0]), rrows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    axes[0].hist(np.clip(finite, -5, 10), bins=80, color="#4c78a8", alpha=.85)
    axes[0].axvline(np.nanmedian(finite), color="k", ls="--", label=f"median {np.nanmedian(finite):.2f}")
    axes[0].set(xlabel=r"Mode Gr\"uneisen parameter $\gamma$", ylabel="Mode count",
                title="Common 12x12x12 mesh (clipped for display)")
    axes[0].legend(frameon=False)
    axes[1].plot([r["temperature_K"] for r in rows], [r["Cv_weighted_gamma"] for r in rows], "o-")
    axes[1].set(xlabel="Temperature (K)", ylabel=r"$C_V$-weighted $\bar\gamma$",
                title="Thermodynamic Gr\"uneisen parameter")
    axes[1].grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(HERE / "gruneisen_parameters.pdf")
    fig.savefig(HERE / "gruneisen_parameters.png", dpi=220)
    plt.close(fig)


def projected_thermodynamics():
    rows = []
    for label in LABELS:
        a = np.loadtxt(FE / f"phonopy_{label}_2x2x2/projected_dos.dat")
        freq, pdos = a[:, 0], a[:, 1:]
        positive = freq > 1e-6
        freq, pdos = freq[positive], pdos[positive]
        groups = {"O": pdos[:, :6].sum(axis=1), "Fe": pdos[:, 6:].sum(axis=1)}
        for temperature in [0, 100, 300, 600, 1000]:
            for species, dos in groups.items():
                zpe = np.trapz(dos * 0.5 * H_THz_KJMOL * freq, freq) / 2
                if temperature == 0:
                    cv = entropy = 0.0
                else:
                    x = 47.9924307 * freq / temperature
                    ex = np.exp(np.clip(x, None, 700))
                    cv_mode = R * x**2 * ex / np.expm1(x) ** 2
                    s_mode = R * (x / np.expm1(x) - np.log1p(-np.exp(-x)))
                    cv = np.trapz(dos * cv_mode, freq) / 2
                    entropy = np.trapz(dos * s_mode, freq) / 2
                rows.append({"volume": label, "temperature_K": temperature, "species": species,
                             "ZPE_kJ_per_mol_Fe2O3": zpe, "Cv_J_per_molK": cv,
                             "S_J_per_molK": entropy})
    write_csv("species_projected_thermodynamics.csv", list(rows[0]), rows)

    selected = [r for r in rows if r["temperature_K"] == 300]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    for ax, field, ylabel in [(axes[0], "Cv_J_per_molK", r"$C_V$ (J mol$^{-1}$ K$^{-1}$)"),
                              (axes[1], "S_J_per_molK", r"$S$ (J mol$^{-1}$ K$^{-1}$)")]:
        x = np.arange(3)
        bottom = np.zeros(3)
        for species, color in [("Fe", "#b44"), ("O", "#4c78a8")]:
            vals = [next(r[field] for r in selected if r["volume"] == lab and r["species"] == species) for lab in LABELS]
            ax.bar(x, vals, bottom=bottom, label=species, color=color)
            bottom += vals
        ax.set_xticks(x, LABELS); ax.set_ylabel(ylabel); ax.grid(axis="y", alpha=.2)
    axes[0].legend(frameon=False); fig.suptitle("Species-projected harmonic thermodynamics at 300 K")
    fig.tight_layout(); fig.savefig(HERE / "species_projected_thermodynamics.pdf")
    fig.savefig(HERE / "species_projected_thermodynamics.png", dpi=220); plt.close(fig)


def displacements_and_velocities():
    temperatures = np.arange(0, 1001, 25)
    msd_rows, tensor_rows, velocity_rows = [], [], []
    for label in LABELS:
        folder = FE / f"phonopy_{label}_2x2x2"
        ph = load(str(folder / "phonopy_disp.yaml"), force_constants_filename=str(folder / "FORCE_CONSTANTS"))
        ph.run_mesh([12, 12, 12], is_mesh_symmetry=False, with_eigenvectors=True)
        ph.run_thermal_displacements(temperatures=temperatures)
        td = ph.get_thermal_displacements_dict()["thermal_displacements"].reshape(len(temperatures), 10, 3)
        symbols = np.array(ph.primitive.symbols)
        for it, temp in enumerate(temperatures):
            for species in ["Fe", "O"]:
                vals = td[it, symbols == species].sum(axis=1)
                msd_rows.append({"volume": label, "temperature_K": temp, "species": species,
                                 "mean_Uiso_A2": vals.mean() / 3, "mean_total_MSD_A2": vals.mean()})
        ph.run_thermal_displacement_matrices(temperatures=[300])
        mats = ph.get_thermal_displacement_matrices_dict()["thermal_displacement_matrices"][0]
        for ia, (species, mat) in enumerate(zip(symbols, mats), 1):
            tensor_rows.append({"volume": label, "atom_index": ia, "species": species,
                                "Uxx_A2": mat[0, 0], "Uyy_A2": mat[1, 1], "Uzz_A2": mat[2, 2],
                                "Uxy_A2": mat[0, 1], "Uxz_A2": mat[0, 2], "Uyz_A2": mat[1, 2],
                                "Uiso_A2": np.trace(mat) / 3})

        # Lattice vectors are rows; reciprocal basis vectors are rows of A^-T.
        rec = np.linalg.inv(ph.primitive.cell).T
        directions = [("a*", np.array([1., 0, 0])), ("b*", np.array([0., 1, 0])),
                      ("c*", np.array([0., 0, 1]))]
        delta = 0.001
        ph.run_qpoints([delta * x[1] for x in directions], with_group_velocities=True)
        qd = ph.get_qpoints_dict()
        for idir, (name, qred) in enumerate(directions):
            qcart = qred @ rec
            unit = qcart / np.linalg.norm(qcart)
            for branch in range(3):
                # QE phonopy cells have physical length unit au, so the raw
                # velocity is THz*bohr; 1 THz*bohr = 0.0529177 km/s.
                velocity_factor = 0.1 * BOHR_TO_ANGSTROM
                phase = qd["frequencies"][idir, branch] / (delta * np.linalg.norm(qcart)) * velocity_factor
                group = abs(np.dot(qd["group_velocities"][idir, branch], unit)) * velocity_factor
                velocity_rows.append({"volume": label, "propagation_direction": name,
                                      "acoustic_branch_sorted": branch + 1,
                                      "phase_velocity_km_s": phase, "group_velocity_projection_km_s": group,
                                      "q_reduced_magnitude": delta})
    write_csv("mean_square_displacements.csv", list(msd_rows[0]), msd_rows)
    write_csv("displacement_tensors_300K.csv", list(tensor_rows[0]), tensor_rows)
    write_csv("acoustic_velocities.csv", list(velocity_rows[0]), velocity_rows)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for label, color in zip(LABELS, COLORS):
        for species, ls in [("Fe", "-"), ("O", "--")]:
            ss = [r for r in msd_rows if r["volume"] == label and r["species"] == species]
            ax.plot([r["temperature_K"] for r in ss], [r["mean_total_MSD_A2"] for r in ss],
                    color=color, ls=ls, label=f"{label} {species}")
    ax.set(xlabel="Temperature (K)", ylabel=r"Mean total MSD ($\AA^2$)",
           title=r"Harmonic atomic mean-square displacement")
    ax.grid(alpha=.25); ax.legend(frameon=False, ncol=2); fig.tight_layout()
    fig.savefig(HERE / "mean_square_displacements.pdf"); fig.savefig(HERE / "mean_square_displacements.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11, 4), sharey=True)
    for ax, lab in zip(axes, LABELS):
        subset = [r for r in velocity_rows if r["volume"] == lab]
        for b, marker in zip([1, 2, 3], ["o", "s", "^"]):
            ss = [r for r in subset if r["acoustic_branch_sorted"] == b]
            ax.plot([r["propagation_direction"] for r in ss], [r["group_velocity_projection_km_s"] for r in ss], marker=marker, label=f"branch {b}")
        ax.set_title(lab); ax.set_xlabel("Propagation direction"); ax.grid(alpha=.2)
    axes[0].set_ylabel("Projected group velocity (km/s)"); axes[0].legend(frameon=False)
    fig.suptitle("Near-Gamma acoustic velocities"); fig.tight_layout()
    fig.savefig(HERE / "acoustic_velocities.pdf"); fig.savefig(HERE / "acoustic_velocities.png", dpi=220); plt.close(fig)


def qha():
    eos, eos_p, eos_err = load_eos()
    therm = [load_thermal(label) for label in LABELS]
    temps = therm[0]["T"]
    for t in therm[1:]:
        if not np.allclose(t["T"], temps): raise RuntimeError("Thermal temperature grids differ")

    qha_rows = []
    veqs, bts, cvs = [], [], []
    coeffs_all = []
    for i, temp in enumerate(temps):
        fvib = np.array([x["F"][i] for x in therm]) * EV_PER_KJMOL
        coef = np.polyfit(VOLUMES, fvib, 2)
        coeffs_all.append(coef)
        total = lambda v: bm3(v, *eos_p) + np.polyval(coef, v)
        opt = minimize_scalar(total, bounds=(94, 114), method="bounded", options={"xatol": 1e-10})
        v = opt.x
        h = 0.01
        curv = (total(v + h) - 2 * total(v) + total(v - h)) / h**2
        bt = v * curv * 160.21766208
        cv = np.polyval(np.polyfit(VOLUMES, [x["Cv"][i] for x in therm], 2), v) / 2
        veqs.append(v); bts.append(bt); cvs.append(cv)

    veqs, bts, cvs = map(np.asarray, (veqs, bts, cvs))
    alpha_v = np.gradient(np.log(veqs), temps, edge_order=2)
    alpha_l = alpha_v / 3
    molar_volume_formula = veqs / 2 * 1e-30 * NA
    cp_minus_cv = alpha_v**2 * (bts * 1e9) * molar_volume_formula * temps
    cps = cvs + cp_minus_cv
    for i, temp in enumerate(temps):
        qha_rows.append({"temperature_K": temp, "Veq_A3_per_10atom_cell": veqs[i],
                         "phonon_volume_extrapolated": not (VPH_MIN <= veqs[i] <= VPH_MAX),
                         "alpha_volume_1e-6_K-1": alpha_v[i] * 1e6,
                         "alpha_linear_1e-6_K-1": alpha_l[i] * 1e6,
                         "BT_GPa": bts[i], "Cv_J_per_molK_Fe2O3": cvs[i],
                         "Cp_J_per_molK_Fe2O3": cps[i], "Cp_minus_Cv_J_per_molK": cp_minus_cv[i]})
    write_csv("preliminary_qha_thermodynamics.csv", list(qha_rows[0]), qha_rows)

    pgrid_rows = []
    for temp in [0, 300, 600, 900]:
        i = int(temp / 10)
        coef = coeffs_all[i]
        for pressure in [0, 1, 2, 5, 10]:
            func = lambda v: bm3(v, *eos_p) + np.polyval(coef, v) + pressure * v * GPA_A3_TO_EV
            opt = minimize_scalar(func, bounds=(88, 114), method="bounded")
            pgrid_rows.append({"temperature_K": temp, "pressure_GPa": pressure,
                               "Veq_A3_per_10atom_cell": opt.x,
                               "G_relative_eV_per_cell": opt.fun,
                               "phonon_volume_extrapolated": not (VPH_MIN <= opt.x <= VPH_MAX)})
    base = min(r["G_relative_eV_per_cell"] for r in pgrid_rows)
    for r in pgrid_rows: r["G_relative_eV_per_cell"] -= base
    write_csv("preliminary_qha_pressure_grid.csv", list(pgrid_rows[0]), pgrid_rows)

    volume_rows = []
    for label, vol, th in zip(LABELS, VOLUMES, therm):
        for temp in [0, 300, 600, 1000]:
            i = int(temp / 10)
            volume_rows.append({"volume": label, "volume_A3": vol, "temperature_K": temp,
                                "Fvib_kJ_per_mol_Fe2O3": th["F"][i] / 2,
                                "entropy_J_per_molK_Fe2O3": th["S"][i] / 2,
                                "Cv_J_per_molK_Fe2O3": th["Cv"][i] / 2,
                                "ZPE_kJ_per_mol_Fe2O3": th["ZPE"] / 2})
    write_csv("volume_dependent_thermodynamics.csv", list(volume_rows[0]), volume_rows)

    static_v0 = eos_p[1]
    summary = [
        {"quantity": "static_V0", "value": static_v0, "units": "A3/10-atom cell", "note": "BM3, 9 static volumes"},
        {"quantity": "static_B0", "value": eos_p[2] * 160.21766208, "units": "GPa", "note": "isothermal 0 K BM3"},
        {"quantity": "static_B0_prime", "value": eos_p[3], "units": "dimensionless", "note": "BM3"},
        {"quantity": "QHA_V0_with_ZPE", "value": veqs[0], "units": "A3/10-atom cell", "note": "extrapolative"},
        {"quantity": "zero_point_volume_expansion", "value": 100 * (veqs[0] / static_v0 - 1), "units": "percent", "note": "extrapolative"},
        {"quantity": "QHA_V300", "value": veqs[30], "units": "A3/10-atom cell", "note": "extrapolative"},
        {"quantity": "QHA_B300", "value": bts[30], "units": "GPa", "note": "extrapolative"},
        {"quantity": "QHA_alpha_linear_300", "value": alpha_l[30] * 1e6, "units": "1e-6 K-1", "note": "extrapolative"},
        {"quantity": "QHA_Cp300", "value": cps[30], "units": "J mol-1 K-1", "note": "per Fe2O3; extrapolative"},
    ]
    write_csv("equation_of_state_qha_summary.csv", list(summary[0]), summary)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    axes[0,0].plot(temps, veqs); axes[0,0].axhspan(VPH_MIN, VPH_MAX, color="green", alpha=.12, label="phonon volume range")
    axes[0,0].set(ylabel=r"$V_{eq}$ ($\AA^3$/cell)", xlabel="Temperature (K)"); axes[0,0].legend(frameon=False)
    axes[0,1].plot(temps, alpha_l * 1e6); axes[0,1].set(ylabel=r"Linear $\alpha$ ($10^{-6}$ K$^{-1}$)", xlabel="Temperature (K)")
    axes[1,0].plot(temps, bts); axes[1,0].set(ylabel=r"$B_T$ (GPa)", xlabel="Temperature (K)")
    axes[1,1].plot(temps, cvs, label=r"$C_V$"); axes[1,1].plot(temps, cps, label=r"$C_P$")
    axes[1,1].set(ylabel=r"Heat capacity (J mol$^{-1}$ K$^{-1}$)", xlabel="Temperature (K)"); axes[1,1].legend(frameon=False)
    for ax in axes.flat: ax.grid(alpha=.25)
    fig.suptitle("Preliminary three-volume QHA (extrapolative)"); fig.tight_layout()
    fig.savefig(HERE / "preliminary_qha_thermodynamics.pdf"); fig.savefig(HERE / "preliminary_qha_thermodynamics.png", dpi=220); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.7))
    for temp in [0, 300, 600, 900]:
        ss = [r for r in pgrid_rows if r["temperature_K"] == temp]
        ax.plot([r["pressure_GPa"] for r in ss], [r["Veq_A3_per_10atom_cell"] for r in ss], "o-", label=f"{temp} K")
    ax.axhspan(VPH_MIN, VPH_MAX, color="green", alpha=.12); ax.set(xlabel="Pressure (GPa)", ylabel=r"$V_{eq}$ ($\AA^3$/cell)", title="Preliminary QHA pressure-temperature grid")
    ax.grid(alpha=.25); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(HERE / "preliminary_qha_pressure_grid.pdf"); fig.savefig(HERE / "preliminary_qha_pressure_grid.png", dpi=220); plt.close(fig)


def literature_table():
    qha = {r["quantity"]: r for r in csv.DictReader((HERE / "equation_of_state_qha_summary.csv").open())}
    thermal = list(csv.DictReader((HERE / "volume_dependent_thermodynamics.csv").open()))
    s102 = next(r for r in thermal if r["volume"] == "V102" and r["temperature_K"] == "300")
    rows = [
        {"characteristic": "entropy_near_300K", "this_work": s102["entropy_J_per_molK_Fe2O3"], "literature": 87.32, "units": "J mol-1 K-1", "comparison_basis": "V102 harmonic S(300 K) vs calorimetric S(298.15 K)", "source": "Snow et al. 2010, doi:10.1016/j.jct.2010.04.010"},
        {"characteristic": "Cp_300K", "this_work": qha["QHA_Cp300"]["value"], "literature": 104.155, "units": "J mol-1 K-1", "comparison_basis": "extrapolative QHA vs NIST-JANAF Shomate at 300 K", "source": "NIST-JANAF (Chase 1998), NIST WebBook SRD 69"},
        {"characteristic": "bulk_modulus", "this_work": qha["static_B0"]["value"], "literature": 206.6, "units": "GPa", "comparison_basis": "0 K static isothermal BM3 vs ambient ultrasonic elastic aggregate", "source": "Liebermann et al., summarized in Zhou et al. 2025"},
    ]
    write_csv("extended_literature_comparison.csv", list(rows[0]), rows)

    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    for ax, row, title in zip(axes, rows, ["Entropy near 300 K", r"$C_P$ at 300 K", "Bulk modulus"]):
        vals = [float(row["literature"]), float(row["this_work"])]
        bars = ax.bar(["Literature", "This work"], vals, color=["#666", "#1f77b4"])
        ax.bar_label(bars, fmt="%.1f"); ax.set_title(title); ax.set_ylabel(row["units"]); ax.grid(axis="y", alpha=.2)
    fig.suptitle(r"Extended benchmarks for $\alpha$-Fe$_2$O$_3$"); fig.tight_layout()
    fig.savefig(HERE / "extended_literature_comparison.pdf"); fig.savefig(HERE / "extended_literature_comparison.png", dpi=220); plt.close(fig)


if __name__ == "__main__":
    mode_gruneisen()
    projected_thermodynamics()
    displacements_and_velocities()
    qha()
    literature_table()
    print(f"Extended analysis written to {HERE}")
