#!/usr/bin/env python3
"""Extract QE EOS data and fit a third-order Birch-Murnaghan equation."""

from pathlib import Path
import csv
import re
import sys

import numpy as np


POINTS = ["V092", "V094", "V096", "V098", "V100", "V102", "V104", "V106", "V108"]
RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
EV_A3_TO_GPA = 160.21766208
V_REFERENCE = 100.624948
AHEX_REFERENCE = 5.0355
CHEX_REFERENCE = 13.7471


def last_float(pattern: str, text: str, name: str) -> float:
    matches = re.findall(pattern, text, flags=re.I)
    if not matches:
        raise ValueError(f"could not find {name}")
    return float(matches[-1])


def read_output(point: str, path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(path)
    text = path.read_text(errors="replace")
    if "JOB DONE" not in text:
        raise RuntimeError(f"{path} does not contain JOB DONE")
    energy_ry = last_float(r"!\s+total energy\s+=\s+([-+0-9.Ee]+)\s+Ry", text, "energy")
    volume_bohr3 = last_float(r"unit-cell volume\s+=\s+([-+0-9.Ee]+)", text, "volume")
    force = last_float(r"Total force\s+=\s+([-+0-9.Ee]+)", text, "final force")
    total_mag = last_float(r"total magnetization\s+=\s+([-+0-9.Ee]+)", text, "magnetization")
    absolute_mag = last_float(r"absolute magnetization\s+=\s+([-+0-9.Ee]+)", text, "absolute magnetization")
    volume_A3 = volume_bohr3 * BOHR_TO_ANG**3
    target_ratio = int(point[1:]) / 100.0
    target_volume = V_REFERENCE * target_ratio
    if abs(volume_A3 / target_volume - 1.0) > 0.005:
        raise ValueError(
            f"{path} volume {volume_A3:.6f} A^3 does not match "
            f"{point} target {target_volume:.6f} A^3"
        )
    return {
        "point": point,
        "source_file": str(path),
        "volume_A3": volume_A3,
        "energy_eV": energy_ry * RY_TO_EV,
        "force_Ry_bohr": force,
        "total_magnetization": total_mag,
        "absolute_magnetization": absolute_mag,
    }


def read_point(point: str) -> dict:
    candidates = [Path(point) / "alpha_Fe2O3.eos.out"]
    if point == "V100":
        candidates.append(Path("..") / "alpha_Fe2O3.relax.out")
    failures = []
    for path in candidates:
        try:
            return read_output(point, path)
        except Exception as exc:
            failures.append(str(exc))
    raise RuntimeError("; ".join(failures))


def birch_murnaghan(volume: np.ndarray, e0: float, v0: float, b0: float, bp: float) -> np.ndarray:
    x = (v0 / volume) ** (2.0 / 3.0)
    return e0 + 9.0 * v0 * b0 / 16.0 * (
        bp * (x - 1.0) ** 3 + (x - 1.0) ** 2 * (6.0 - 4.0 * x)
    )


def grid_fallback(volume: np.ndarray, energy: np.ndarray) -> tuple[np.ndarray, str]:
    """Numpy-only fallback if scipy is unavailable."""
    v_center = float(volume[np.argmin(energy)])
    bp_center = 4.0
    v_span = 0.08 * v_center
    bp_span = 3.0
    best = None
    for _ in range(5):
        for v0 in np.linspace(v_center - v_span, v_center + v_span, 301):
            x = (v0 / volume) ** (2.0 / 3.0)
            for bp in np.linspace(max(0.5, bp_center - bp_span), bp_center + bp_span, 241):
                f = 9.0 * v0 / 16.0 * (
                    bp * (x - 1.0) ** 3 + (x - 1.0) ** 2 * (6.0 - 4.0 * x)
                )
                denominator = float(np.dot(f - f.mean(), f - f.mean()))
                if denominator == 0.0:
                    continue
                b0 = float(np.dot(f - f.mean(), energy - energy.mean()) / denominator)
                if b0 <= 0.0:
                    continue
                e0 = float(energy.mean() - b0 * f.mean())
                residual = energy - (e0 + b0 * f)
                sse = float(np.dot(residual, residual))
                if best is None or sse < best[0]:
                    best = (sse, e0, v0, b0, bp)
        if best is None:
            raise RuntimeError("numpy fallback could not find a physical fit")
        _, _, v_center, _, bp_center = best
        v_span *= 0.18
        bp_span *= 0.18
    return np.array(best[1:], dtype=float), "numpy grid fallback"


def fit_eos(volume: np.ndarray, energy: np.ndarray) -> tuple[np.ndarray, str]:
    try:
        from scipy.optimize import curve_fit

        guess = [float(energy.min()), float(volume[np.argmin(energy)]), 1.2, 4.0]
        lower = [energy.min() - 10.0, volume.min() * 0.85, 1.0e-6, 0.5]
        upper = [energy.min() + 10.0, volume.max() * 1.15, 10.0, 12.0]
        parameters, _ = curve_fit(
            birch_murnaghan, volume, energy, p0=guess,
            bounds=(lower, upper), maxfev=100000,
        )
        return parameters, "scipy.optimize.curve_fit"
    except (ImportError, RuntimeError, ValueError) as exc:
        print(f"SciPy fit unavailable ({exc}); using numpy fallback.")
        return grid_fallback(volume, energy)


def main() -> None:
    rows = []
    errors = []
    for point in POINTS:
        try:
            rows.append(read_point(point))
        except Exception as exc:
            errors.append(f"{point}: {exc}")
    if errors:
        print("Skipping unavailable or incomplete points:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
    if len(rows) < 5:
        print(
            f"Need at least five completed points for a four-parameter fit; found {len(rows)}.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    rows.sort(key=lambda row: float(row["volume_A3"]))

    volume = np.array([float(row["volume_A3"]) for row in rows])
    energy = np.array([float(row["energy_eV"]) for row in rows])
    parameters, method = fit_eos(volume, energy)
    e0, v0, b0, bp = parameters
    fitted = birch_murnaghan(volume, *parameters)
    rmse_mev_cell = 1000.0 * float(np.sqrt(np.mean((energy - fitted) ** 2)))
    scale = (v0 / V_REFERENCE) ** (1.0 / 3.0)
    ahex0 = AHEX_REFERENCE * scale
    chex0 = CHEX_REFERENCE * scale
    minimum_index = int(np.argmin(energy))
    minimum_bracketed = 0 < minimum_index < len(energy) - 1
    fitted_v0_bracketed = volume.min() < v0 < volume.max()
    status = "FINAL nine-point dataset" if len(rows) == len(POINTS) else "PRELIMINARY partial dataset"
    warnings = []
    if not minimum_bracketed:
        warnings.append("The lowest calculated energy is at an endpoint; the minimum is not bracketed.")
    if not fitted_v0_bracketed:
        warnings.append("The fitted V0 is outside the completed volume range (extrapolation).")
    if len(rows) < 7:
        warnings.append("Fewer than seven points make a four-parameter third-order fit weakly constrained.")

    relative_mev_formula = (energy - energy.min()) * 1000.0 / 2.0
    with Path("eos_data.csv").open("w", newline="") as handle:
        fields = [
            "point", "source_file", "volume_A3", "energy_eV", "relative_meV_per_Fe2O3",
            "force_Ry_bohr", "total_magnetization", "absolute_magnetization",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row, relative in zip(rows, relative_mev_formula):
            output = dict(row)
            output["relative_meV_per_Fe2O3"] = relative
            writer.writerow(output)

    result = (
        "Third-order Birch-Murnaghan fit for AFM alpha-Fe2O3\n"
        f"Status                     : {status}\n"
        f"Fit method                 : {method}\n"
        f"Points                     : {len(rows)}\n"
        f"Completed point labels     : {', '.join(str(row['point']) for row in rows)}\n"
        f"Equilibrium energy/cell    : {e0:.10f} eV\n"
        f"Equilibrium primitive V0   : {v0:.6f} A^3\n"
        f"Equilibrium hexagonal a    : {ahex0:.6f} A\n"
        f"Equilibrium hexagonal c    : {chex0:.6f} A\n"
        f"Bulk modulus B0            : {b0 * EV_A3_TO_GPA:.4f} GPa\n"
        f"Pressure derivative B0'    : {bp:.6f}\n"
        f"Fit RMSE                   : {rmse_mev_cell:.6f} meV/primitive cell\n"
        "Primitive cell contains two Fe2O3 formula units.\n" +
        ("Warnings:\n- " + "\n- ".join(warnings) + "\n" if warnings else "Warnings                   : none\n")
    )
    Path("eos_results.txt").write_text(result)
    print(result)

    try:
        import matplotlib.pyplot as plt

        grid = np.linspace(min(volume.min(), v0) * 0.995, max(volume.max(), v0) * 1.005, 500)
        curve = birch_murnaghan(grid, *parameters)
        reference = energy.min()
        fig, ax = plt.subplots(figsize=(6.2, 4.8))
        ax.scatter(volume, (energy - reference) / 2.0, color="firebrick", zorder=3, label="completed QE points")
        ax.plot(grid, (curve - reference) / 2.0, color="black", lw=1.2, label="preliminary BM fit" if warnings else "3rd-order BM fit")
        ax.axvline(v0, color="0.5", ls="--", lw=0.8)
        ax.set_xlabel(r"Primitive-cell volume ($\mathrm{\AA}^3$)")
        ax.set_ylabel(r"Relative energy (eV / Fe$_2$O$_3$)")
        ax.set_title(r"AFM $\alpha$-Fe$_2$O$_3$ equation of state")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig("alpha_Fe2O3_birch_murnaghan.png", dpi=220)
        fig.savefig("alpha_Fe2O3_birch_murnaghan.pdf")
    except ImportError:
        print("matplotlib is unavailable; numerical fit files were still written.")


if __name__ == "__main__":
    main()
