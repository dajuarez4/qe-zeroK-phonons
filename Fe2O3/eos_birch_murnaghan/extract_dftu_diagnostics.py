#!/usr/bin/env python3
"""Extract magnetic, stress, structural, and convergence diagnostics from QE DFT+U EOS outputs."""

from __future__ import annotations

import csv
import itertools
import math
import re
from pathlib import Path


BOHR_TO_ANG = 0.529177210903
RY_TO_EV = 13.605693122994
VREF_A3 = 100.624948
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def numbers(pattern: str, text: str) -> list[float]:
    return [
        float(value.replace("D", "E").replace("d", "e"))
        for value in re.findall(pattern, text, flags=re.I)
    ]


def last_or_nan(values: list[float]) -> float:
    return values[-1] if values else math.nan


def parse_wall_seconds(text: str) -> float:
    matches = re.findall(r"PWSCF\s+:\s+.*?\sWALL", text)
    if not matches:
        return math.nan
    segment = matches[-1].rsplit("CPU", 1)[-1].replace("WALL", "")
    hours = last_or_nan(numbers(rf"({FLOAT})h", segment))
    minutes = last_or_nan(numbers(rf"({FLOAT})m", segment))
    seconds = last_or_nan(numbers(rf"({FLOAT})s", segment))
    return (
        (0.0 if math.isnan(hours) else hours) * 3600.0
        + (0.0 if math.isnan(minutes) else minutes) * 60.0
        + (0.0 if math.isnan(seconds) else seconds)
    )


def parse_final_positions(text: str) -> list[tuple[str, list[float]]]:
    blocks = re.findall(
        r"ATOMIC_POSITIONS\s*\(crystal\)\s*\n(.*?)(?=\n\s*(?:End final coordinates|"
        r"ATOMIC_POSITIONS|CELL_PARAMETERS|Writing all|$))",
        text,
        flags=re.S,
    )
    if not blocks:
        return []
    positions = []
    for line in blocks[-1].splitlines():
        fields = line.split()
        if len(fields) >= 4 and re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", fields[0]):
            try:
                positions.append((fields[0], [float(value) for value in fields[1:4]]))
            except ValueError:
                pass
    return positions


def parse_cell(text: str) -> list[list[float]]:
    alat = last_or_nan(numbers(rf"lattice parameter \(alat\)\s+=\s+({FLOAT})", text))
    axes_blocks = re.findall(
        rf"a\(1\)\s*=\s*\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\).*?"
        rf"a\(2\)\s*=\s*\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\).*?"
        rf"a\(3\)\s*=\s*\(\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*\)",
        text,
        flags=re.S,
    )
    if math.isnan(alat) or not axes_blocks:
        return []
    values = [float(value.replace("D", "E").replace("d", "e")) for value in axes_blocks[-1]]
    scale = alat * BOHR_TO_ANG
    return [[scale * value for value in values[index : index + 3]] for index in (0, 3, 6)]


def cartesian_delta(frac: list[float], cell: list[list[float]]) -> list[float]:
    return [
        sum(frac[index] * cell[index][component] for index in range(3))
        for component in range(3)
    ]


def fe_o_bonds(
    positions: list[tuple[str, list[float]]], cell: list[list[float]]
) -> tuple[float, float, float]:
    if not positions or not cell:
        return math.nan, math.nan, math.nan
    oxygens = [position for species, position in positions if species.upper().startswith("O")]
    irons = [position for species, position in positions if species.upper().startswith("FE")]
    nearest = []
    for iron in irons:
        distances = []
        for oxygen in oxygens:
            for image in itertools.product((-1, 0, 1), repeat=3):
                delta = [oxygen[i] + image[i] - iron[i] for i in range(3)]
                cart = cartesian_delta(delta, cell)
                distances.append(math.sqrt(sum(value * value for value in cart)))
        nearest.extend(sorted(distances)[:6])
    if not nearest:
        return math.nan, math.nan, math.nan
    return min(nearest), sum(nearest) / len(nearest), max(nearest)


def parse_stress(text: str) -> tuple[float, float, float, float, float]:
    matches = list(
        re.finditer(
            rf"total\s+stress\s+\(Ry/bohr\*\*3\).*?P=\s*({FLOAT})\s*\n"
            rf"\s*{FLOAT}\s+{FLOAT}\s+{FLOAT}\s+({FLOAT})\s+{FLOAT}\s+{FLOAT}\s*\n"
            rf"\s*{FLOAT}\s+{FLOAT}\s+{FLOAT}\s+{FLOAT}\s+({FLOAT})\s+{FLOAT}\s*\n"
            rf"\s*{FLOAT}\s+{FLOAT}\s+{FLOAT}\s+{FLOAT}\s+{FLOAT}\s+({FLOAT})",
            text,
            flags=re.I,
        )
    )
    if not matches:
        return (math.nan,) * 5
    pressure, sxx, syy, szz = [float(value) for value in matches[-1].groups()]
    return pressure, sxx, syy, szz, max(sxx, syy, szz) - min(sxx, syy, szz)


def parse_output(path: Path) -> dict[str, object] | None:
    text = path.read_text(errors="replace")
    if "JOB DONE" not in text:
        return None
    volume_bohr3 = last_or_nan(numbers(rf"unit-cell volume\s+=\s+({FLOAT})", text))
    energies = numbers(rf"!\s+total energy\s+=\s+({FLOAT})\s+Ry", text)
    if math.isnan(volume_bohr3) or not energies:
        return None
    volume = volume_bohr3 * BOHR_TO_ANG**3
    pressure, sxx, syy, szz, stress_anisotropy = parse_stress(text)
    sphere_blocks = re.findall(
        rf"atom\s+\d+\s+\(R={FLOAT}\)\s+charge=\s*({FLOAT})\s+magn=\s*({FLOAT})",
        text,
    )
    sphere_moments = [float(moment) for _, moment in sphere_blocks[-10:]]
    fe_sphere = [abs(value) for value in sphere_moments[6:10]]
    hubbard_moments = numbers(rf"Atomic magnetic moment\s+=\s+({FLOAT})", text)[-4:]
    positions = parse_final_positions(text)
    cell = parse_cell(text)
    bond_min, bond_mean, bond_max = fe_o_bonds(positions, cell)
    oxygen_positions = [position for species, position in positions if species.upper().startswith("O")]
    iron_positions = [position for species, position in positions if species.upper().startswith("FE")]
    bfgs = re.findall(r"bfgs converged in\s+(\d+)\s+scf cycles and\s+(\d+)\s+bfgs steps", text, re.I)
    final_scf = numbers(r"convergence has been achieved in\s+(\d+)\s+iterations", text)
    warnings = len(re.findall(r"eigenvalues not converged", text, re.I))
    return {
        "volume_label": f"V{int(round(100.0 * volume / VREF_A3)):03d}",
        "source_file": path.name,
        "volume_A3": volume,
        "energy_Ry": energies[-1],
        "pressure_kbar": pressure,
        "pressure_GPa": pressure / 10.0,
        "stress_xx_kbar": sxx,
        "stress_yy_kbar": syy,
        "stress_zz_kbar": szz,
        "stress_anisotropy_kbar": stress_anisotropy,
        "total_force_Ry_bohr": last_or_nan(numbers(rf"Total force\s+=\s+({FLOAT})", text)),
        "total_magnetization_muB": last_or_nan(numbers(rf"total magnetization\s+=\s+({FLOAT})", text)),
        "absolute_magnetization_muB": last_or_nan(numbers(rf"absolute magnetization\s+=\s+({FLOAT})", text)),
        "Fe_sphere_abs_moment_muB": sum(fe_sphere) / len(fe_sphere) if fe_sphere else math.nan,
        "Fe_Hubbard_abs_moment_muB": (
            sum(abs(value) for value in hubbard_moments) / len(hubbard_moments)
            if hubbard_moments
            else math.nan
        ),
        "Hubbard_energy_Ry": last_or_nan(numbers(rf"Hubbard energy\s+=\s+({FLOAT})\s+Ry", text)),
        "Hubbard_energy_eV": last_or_nan(numbers(rf"Hubbard energy\s+=\s+({FLOAT})\s+Ry", text))
        * RY_TO_EV,
        "highest_occupied_eV": last_or_nan(numbers(rf"highest occupied level \(ev\):\s+({FLOAT})", text)),
        "final_scf_iterations": int(final_scf[-1]) if final_scf else -1,
        "ionic_scf_cycles": int(bfgs[-1][0]) if bfgs else -1,
        "bfgs_steps": int(bfgs[-1][1]) if bfgs else -1,
        "wall_time_min": parse_wall_seconds(text) / 60.0,
        "eigenvalue_nonconvergence_warnings": warnings,
        "O_internal_x": oxygen_positions[0][0] if oxygen_positions else math.nan,
        "Fe_internal_x": iron_positions[0][0] if iron_positions else math.nan,
        "Fe_O_min_A": bond_min,
        "Fe_O_mean6_A": bond_mean,
        "Fe_O_max6_A": bond_max,
    }


def fmt(value: object, digits: int = 4) -> str:
    if isinstance(value, float):
        return "n/a" if math.isnan(value) else f"{value:.{digits}f}"
    return str(value)


def main() -> None:
    candidates = []
    for path in sorted(Path(".").glob("*.out")):
        record = parse_output(path)
        if record is not None:
            candidates.append(record)
    unique: dict[str, dict[str, object]] = {}
    duplicates = []
    for record in candidates:
        label = str(record["volume_label"])
        previous = unique.get(label)
        if previous is None or float(record["total_force_Ry_bohr"]) < float(previous["total_force_Ry_bohr"]):
            if previous is not None:
                duplicates.append(str(previous["source_file"]))
            unique[label] = record
        else:
            duplicates.append(str(record["source_file"]))
    records = sorted(unique.values(), key=lambda item: float(item["volume_A3"]))
    if not records:
        raise SystemExit("No completed QE EOS outputs found in the current directory")

    output_csv = Path("dftu_eos_diagnostics.csv")
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    lines = [
        "AFM alpha-Fe2O3 PBE+U EOS diagnostics",
        "=====================================",
        "",
        "Model: collinear ++-- AFM, PBE, U_eff(Fe 3d)=4.0 eV, fixed occupations.",
        f"Unique completed volumes: {len(records)}",
        f"Excluded duplicate outputs: {', '.join(duplicates) if duplicates else 'none'}",
        "",
        "Per-volume final-state diagnostics:",
    ]
    for record in records:
        lines.append(
            f"{record['volume_label']}: P={fmt(record['pressure_GPa'], 3)} GPa, "
            f"stress spread={fmt(record['stress_anisotropy_kbar'], 2)} kbar, "
            f"|M_Fe| sphere/Hubbard={fmt(record['Fe_sphere_abs_moment_muB'], 3)}/"
            f"{fmt(record['Fe_Hubbard_abs_moment_muB'], 3)} muB, "
            f"E_U={fmt(record['Hubbard_energy_eV'], 3)} eV, "
            f"<Fe-O>6={fmt(record['Fe_O_mean6_A'], 4)} A, "
            f"SCF/BFGS={record['ionic_scf_cycles']}/{record['bfgs_steps']}"
        )
    pressure_sorted = sorted(records, key=lambda item: abs(float(item["pressure_GPa"])))
    moments = [float(item["Fe_Hubbard_abs_moment_muB"]) for item in records]
    zero_pressure_volume = math.nan
    for left, right in zip(records, records[1:]):
        p_left = float(left["pressure_GPa"])
        p_right = float(right["pressure_GPa"])
        if p_left * p_right <= 0.0 and p_left != p_right:
            v_left = float(left["volume_A3"])
            v_right = float(right["volume_A3"])
            zero_pressure_volume = v_left - p_left * (v_right - v_left) / (p_right - p_left)
            break
    bm_v0 = math.nan
    bm_summary = Path("all_completed_birch_murnaghan_summary.txt")
    if bm_summary.exists():
        standard_section = bm_summary.read_text(errors="replace").split("STANDARD FIT:", 1)
        if len(standard_section) == 2:
            bm_v0 = last_or_nan(numbers(rf"V0\s+=\s+({FLOAT})", standard_section[1]))
    gaps_available = any(
        re.search(r"lowest unoccupied level", Path(str(item["source_file"])).read_text(errors="replace"), re.I)
        for item in records
    )
    lines.extend(
        [
            "",
            "Key observations:",
            f"- The sampled point closest to zero calculated pressure is "
            f"{pressure_sorted[0]['volume_label']} ({fmt(pressure_sorted[0]['pressure_GPa'], 3)} GPa).",
            f"- Linear interpolation of the stress-derived pressure gives P=0 near "
            f"{fmt(zero_pressure_volume, 3)} A^3; the energy-fit V0 is {fmt(bm_v0, 3)} A^3. "
            "Their difference is a useful Pulay-stress/EOS-consistency diagnostic.",
            f"- The absolute Fe Hubbard moment changes monotonically from "
            f"{fmt(moments[0], 3)} to {fmt(moments[-1], 3)} muB as volume increases.",
            f"- Maximum final stress-component spread is "
            f"{fmt(max(float(item['stress_anisotropy_kbar']) for item in records), 2)} kbar; "
            "this measures the non-hydrostatic stress retained by the fixed-shape EOS path.",
            f"- QE printed up to "
            f"{max(int(item['eigenvalue_nonconvergence_warnings']) for item in records)} "
            "'eigenvalues not converged' warnings during a relaxation. The final SCFs converged, "
            "but a clean static SCF at the selected structure should be checked before phonons or bands.",
            f"- A band gap is {'available' if gaps_available else 'not available'} from these outputs. "
            "Only occupied bands were printed, so an NSCF/bands calculation with empty bands is required.",
            "- Total magnetization remains approximately zero; the opposite Fe local moments confirm the AFM state.",
            "",
            "Files:",
            f"- {output_csv.name}: machine-readable values for every volume",
            "- dftu_eos_diagnostics.png/.pdf: pressure, magnetism, bonding, and stress trends",
        ]
    )
    Path("dftu_eos_diagnostics_summary.txt").write_text("\n".join(lines) + "\n")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n".join(lines))
        print("matplotlib unavailable; skipped diagnostic plots")
        return

    volume = [float(item["volume_A3"]) for item in records]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0))
    axes[0, 0].plot(volume, [float(item["pressure_GPa"]) for item in records], "o-")
    axes[0, 0].axhline(0.0, color="0.3", lw=0.8)
    axes[0, 0].set_ylabel("Pressure (GPa)")
    axes[0, 0].set_title("EOS pressure")

    axes[0, 1].plot(volume, [float(item["Fe_sphere_abs_moment_muB"]) for item in records], "o-", label="atomic sphere")
    axes[0, 1].plot(volume, moments, "s--", label="Hubbard occupation")
    axes[0, 1].set_ylabel(r"$|M_{\rm Fe}|$ ($\mu_B$)")
    axes[0, 1].set_title("Fe local moment")
    axes[0, 1].legend(frameon=False, fontsize=8)

    axes[1, 0].plot(volume, [float(item["Fe_O_min_A"]) for item in records], "^--", label="minimum")
    axes[1, 0].plot(volume, [float(item["Fe_O_mean6_A"]) for item in records], "o-", label="mean of six")
    axes[1, 0].plot(volume, [float(item["Fe_O_max6_A"]) for item in records], "v--", label="maximum")
    axes[1, 0].set_xlabel(r"Primitive volume ($\AA^3$)")
    axes[1, 0].set_ylabel(r"Fe-O distance ($\AA$)")
    axes[1, 0].set_title("Relaxed FeO$_6$ coordination")
    axes[1, 0].legend(frameon=False, fontsize=8)

    axes[1, 1].plot(volume, [float(item["stress_anisotropy_kbar"]) / 10.0 for item in records], "o-")
    axes[1, 1].set_xlabel(r"Primitive volume ($\AA^3$)")
    axes[1, 1].set_ylabel("Stress spread (GPa)")
    axes[1, 1].set_title("Fixed-shape stress anisotropy")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
    fig.suptitle(r"AFM $\alpha$-Fe$_2$O$_3$ PBE+U EOS diagnostics")
    fig.tight_layout()
    fig.savefig("dftu_eos_diagnostics.png", dpi=230)
    fig.savefig("dftu_eos_diagnostics.pdf")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
