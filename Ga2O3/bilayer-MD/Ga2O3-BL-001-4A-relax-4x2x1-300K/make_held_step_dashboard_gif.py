#!/usr/bin/env python3
"""Render a frame-by-frame HELD dispersion and MD diagnostic dashboard."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import imageio.v2 as imageio
import matplotlib
import numpy as np
from ase.data import atomic_masses, atomic_numbers

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ADAPTER_DIR = (Path(__file__).resolve().parent / "../Ga2O3-BL-001-300K").resolve()
PARSER_DIR = (Path(__file__).resolve().parent / "../Ga2O3-BL-001-4x2x1-300K").resolve()
sys.path.insert(0, str(PARSER_DIR))
sys.path.insert(0, str(ADAPTER_DIR))

from run_held_phonons import MultiSpeciesHeldModel, reciprocal_path  # noqa: E402
from qe_supercell_to_tdep import parse_qe_input, parse_qe_output  # noqa: E402
from prepare_tdep import parse_initial_supercell, primitive_from_replication  # noqa: E402


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "HELD_300K"
NAVY = "#10243e"
BLUE = "#2b6cb0"
RED = "#c23b2a"
LIGHT = "#f4f7fb"
GRAY = "#586575"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fps", type=float, default=8.0)
    parser.add_argument("--dpi", type=int, default=100)
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTDIR / "Ga2O3-BL-001-4A-HELD-step-dashboard.gif",
    )
    return parser.parse_args()


def build_model_and_data():
    output = ROOT / "Ga2O3-BL-001-4x2x1.md.out"
    unit_cell, _symbols, _positions = parse_qe_input(
        ROOT / "Ga2O3-BL-001.relax.in"
    )
    super_cell = [
        [value * 4 for value in unit_cell[0]],
        [value * 2 for value in unit_cell[1]],
        list(unit_cell[2]),
    ]
    super_symbols, super_reference = parse_initial_supercell(output)
    unit_symbols, unit_reference, _error = primitive_from_replication(
        super_symbols, super_reference
    )
    positions, forces, potential, _kinetic, temperatures, timestep, excluded = (
        parse_qe_output(output, len(super_symbols))
    )
    unit_cell = np.asarray(unit_cell, dtype=float)
    super_cell = np.asarray(super_cell, dtype=float)
    positions = np.asarray(positions, dtype=float)
    forces = np.asarray(forces, dtype=float)
    numbers = np.asarray([atomic_numbers[symbol] for symbol in unit_symbols])
    masses = np.asarray([atomic_masses[number] for number in numbers])
    model = MultiSpeciesHeldModel(
        primitive_cell=unit_cell,
        basis_frac=np.asarray(unit_reference),
        supercell_cell=super_cell,
        ideal_supercell_frac=np.asarray(super_reference),
        primitive_numbers=numbers,
        primitive_masses=masses,
        cutoff_ang=2.5,
    )
    return {
        "model": model,
        "unit_cell": unit_cell,
        "super_cell": super_cell,
        "symbols": np.asarray(super_symbols),
        "positions": positions,
        "forces": forces,
        "potential": np.asarray(potential),
        "temperatures": np.asarray(temperatures),
        "timestep_fs": float(timestep),
        "excluded": int(excluded),
    }


def prepare_cache(data, cache_path, force_recompute):
    coefficients = np.loadtxt(
        OUTDIR / "held_step_coefficients.csv", delimiter=",", skiprows=1
    )
    coefficients = np.atleast_2d(coefficients)
    nframes = len(data["positions"])
    if len(coefficients) != nframes:
        raise ValueError(
            f"Coefficient/trajectory mismatch: {len(coefficients)} versus {nframes}"
        )
    if cache_path.exists() and not force_recompute:
        cached = np.load(cache_path, allow_pickle=False)
        if int(cached["nframes"]) == nframes:
            return {name: np.asarray(cached[name]) for name in cached.files}

    model = data["model"]
    q_path, x_values, _labels, _ticks = reciprocal_path(data["unit_cell"])
    step_frequencies = np.empty((nframes, len(q_path), 30))
    force_rmse = np.empty(nframes)
    for index, (frame, frame_forces, values) in enumerate(
        zip(data["positions"], data["forces"], coefficients)
    ):
        step_frequencies[index] = model.dispersion_thz_from_reduced_path(
            values, q_path
        )
        design = model.build_design_matrix(model.frame_displacements_cart(frame))
        force_rmse[index] = np.sqrt(
            np.mean((design @ values - frame_forces.reshape(-1)) ** 2)
        )
        if index == 0 or (index + 1) % 20 == 0 or index + 1 == nframes:
            print(f"[HELD dashboard] evaluated {index + 1}/{nframes}")

    mean_frequencies = model.dispersion_thz_from_reduced_path(
        coefficients.mean(axis=0), q_path
    )
    tick_indices = np.asarray([0, 80, 160, 240, 320])
    payload = {
        "nframes": np.asarray(nframes),
        "x_values": x_values,
        "tick_positions": x_values[tick_indices],
        "step_frequencies_thz": step_frequencies,
        "mean_frequencies_thz": mean_frequencies,
        "force_rmse_eV_per_A": force_rmse,
    }
    np.savez_compressed(cache_path, **payload)
    return payload


def draw_structure(ax, fractional, cell, symbols):
    cart = (fractional - np.floor(fractional)) @ cell
    colors = {"Ga": "#386cb0", "O": "#8ecae6"}
    for symbol in sorted(set(symbols.tolist())):
        mask = symbols == symbol
        ax.scatter(
            cart[mask, 0], cart[mask, 1], cart[mask, 2],
            s=25, color=colors[symbol], edgecolors="white", linewidths=0.25,
            depthshade=True, label=symbol,
        )
    vertices = np.asarray([np.zeros(3), cell[0], cell[1], cell[0] + cell[1]])
    x_min, y_min = vertices[:, :2].min(axis=0)
    x_max, y_max = vertices[:, :2].max(axis=0)
    z_center = 0.5 * np.linalg.norm(cell[2])
    ax.set_xlim(x_min - 0.5, x_max + 0.5)
    ax.set_ylim(y_min - 0.5, y_max + 0.5)
    ax.set_zlim(z_center - 5.3, z_center + 5.3)
    ax.set_box_aspect((x_max - x_min + 1, y_max - y_min + 1, 10.6))
    ax.view_init(elev=22, azim=-67)
    ax.set_axis_off()
    ax.legend(frameon=False, fontsize=7, loc="lower left")


def draw_metrics(ax, frame, nframes, time_fs, temperature, rmse, frequencies):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#d9e1eb")
    negative = int(np.sum(frequencies < -0.1))
    lines = [
        ("MD / HELD STEP", f"{frame + 1} / {nframes}", BLUE),
        ("TIME", f"{time_fs:.2f} fs", BLUE),
        ("TEMPERATURE", f"{temperature:.1f} K", RED if abs(temperature - 300) > 50 else BLUE),
        ("FORCE RMSE", f"{rmse:.3f} eV/Å", RED),
        ("FREQUENCY RANGE", f"{frequencies.min():.1f} to {frequencies.max():.1f} THz", RED),
        ("VALUES < −0.1 THz", f"{negative:,}", RED if negative else BLUE),
    ]
    for row, (label, value, color) in enumerate(lines):
        y = 0.91 - row * 0.155
        ax.text(0.05, y, label, fontsize=7.4, color=GRAY, weight="bold")
        ax.text(0.95, y, value, fontsize=8.8, color=color, weight="bold", ha="right")


def render(args, data, cache):
    frequencies = cache["step_frequencies_thz"]
    mean_frequencies = cache["mean_frequencies_thz"]
    x_values = cache["x_values"]
    ticks = cache["tick_positions"]
    force_rmse = cache["force_rmse_eV_per_A"]
    temperatures = data["temperatures"]
    nframes = len(temperatures)
    times = np.arange(1, nframes + 1) * data["timestep_fs"]
    robust_min, robust_max = np.percentile(frequencies, [0.35, 99.65])
    y_min = max(-120.0, min(-5.0, float(robust_min) * 1.05))
    y_max = min(120.0, max(15.0, float(robust_max) * 1.05))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    last_png = args.output.with_name(args.output.stem + "-last-frame.png")

    with imageio.get_writer(
        args.output, mode="I", duration=int(round(1000 / args.fps)), loop=0
    ) as writer:
        for frame in range(nframes):
            fig = plt.figure(figsize=(14.4, 8.1), facecolor=LIGHT)
            header = fig.add_axes([0, 0.90, 1, 0.10])
            header.set_facecolor(NAVY)
            header.set_xticks([])
            header.set_yticks([])
            for spine in header.spines.values():
                spine.set_visible(False)
            header.text(
                0.035, 0.55, "Ga₂O₃ 4 Å bilayer · experimental HELD dashboard",
                color="white", fontsize=18, weight="bold", va="center",
            )
            header.text(
                0.965, 0.55, "80 atoms · 2.5 Å cutoff", color="#cdd9e8",
                fontsize=10, ha="right", va="center",
            )

            ax_disp = fig.add_axes([0.055, 0.10, 0.61, 0.75], facecolor="white")
            for position in ticks:
                ax_disp.axvline(position, color="#aab2bd", linewidth=0.7)
            for band in range(30):
                ax_disp.plot(
                    x_values, mean_frequencies[:, band], color="#7b8794",
                    linestyle="--", linewidth=0.65, alpha=0.48,
                )
                ax_disp.plot(
                    x_values, frequencies[frame, :, band], color=RED,
                    linewidth=1.0, alpha=0.94,
                )
            ax_disp.axhline(0, color="#202832", linewidth=0.9)
            ax_disp.set_xlim(x_values[0], x_values[-1])
            ax_disp.set_ylim(y_min, y_max)
            ax_disp.set_xticks(ticks, ["Γ", "X", "S", "Y", "Γ"])
            ax_disp.set_ylabel("Frequency (THz)")
            ax_disp.set_title(
                f"Per-step HELD phonons · frame {frame + 1}", loc="left",
                color=NAVY, fontsize=12, weight="bold",
            )
            ax_disp.grid(axis="y", alpha=0.14)
            ax_disp.legend(
                handles=[
                    Line2D([0], [0], color=RED, lw=1.4, label="Current step"),
                    Line2D([0], [0], color="#7b8794", lw=1, ls="--", label=f"{nframes}-frame mean"),
                ], frameon=False, loc="upper right", fontsize=9,
            )

            ax_metrics = fig.add_axes([0.70, 0.655, 0.265, 0.195])
            draw_metrics(
                ax_metrics, frame, nframes, times[frame], temperatures[frame],
                force_rmse[frame], frequencies[frame],
            )
            ax_temp = fig.add_axes([0.70, 0.43, 0.265, 0.17], facecolor="white")
            ax_temp.plot(times, temperatures, color="#c5cbd3", linewidth=1.2)
            ax_temp.plot(times[:frame + 1], temperatures[:frame + 1], color=RED, linewidth=1.8)
            ax_temp.scatter(times[frame], temperatures[frame], color=RED, s=22, zorder=3)
            ax_temp.axhline(300, color=BLUE, linestyle="--", linewidth=1)
            ax_temp.set_title("Temperature history", loc="left", fontsize=9, color=NAVY, weight="bold")
            ax_temp.set_ylabel("K", fontsize=8)
            ax_temp.tick_params(labelsize=7)
            ax_temp.grid(alpha=0.14)

            ax_rmse = fig.add_axes([0.70, 0.235, 0.265, 0.14], facecolor="white")
            ax_rmse.plot(times, force_rmse, color="#c5cbd3", linewidth=1.2)
            ax_rmse.plot(times[:frame + 1], force_rmse[:frame + 1], color=RED, linewidth=1.7)
            ax_rmse.scatter(times[frame], force_rmse[frame], color=RED, s=20, zorder=3)
            ax_rmse.set_title("Per-step force-fit RMSE", loc="left", fontsize=9, color=NAVY, weight="bold")
            ax_rmse.set_xlabel("Time (fs)", fontsize=8)
            ax_rmse.set_ylabel("eV/Å", fontsize=8)
            ax_rmse.tick_params(labelsize=7)
            ax_rmse.grid(alpha=0.14)

            ax_structure = fig.add_axes([0.735, 0.015, 0.19, 0.19], projection="3d")
            draw_structure(ax_structure, data["positions"][frame], data["super_cell"], data["symbols"])
            fig.canvas.draw()
            writer.append_data(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])
            if frame + 1 == nframes:
                fig.savefig(last_png, dpi=args.dpi, facecolor=fig.get_facecolor())
            plt.close(fig)
            if frame == 0 or (frame + 1) % 20 == 0 or frame + 1 == nframes:
                print(f"[HELD dashboard] rendered {frame + 1}/{nframes}")

    summary = {
        "status": "experimental_multispecies_HELD_step_dashboard",
        "frames": nframes,
        "duration_fs": float(times[-1]),
        "gif_fps": float(args.fps),
        "mean_temperature_K": float(np.mean(temperatures)),
        "mean_step_force_rmse_eV_per_A": float(np.mean(force_rmse)),
        "minimum_step_force_rmse_eV_per_A": float(np.min(force_rmse)),
        "maximum_step_force_rmse_eV_per_A": float(np.max(force_rmse)),
        "global_step_frequency_minimum_THz": float(np.min(frequencies)),
        "global_step_frequency_maximum_THz": float(np.max(frequencies)),
        "display_y_min_THz": y_min,
        "display_y_max_THz": y_max,
        "gif": args.output.name,
        "last_frame_png": last_png.name,
    }
    (OUTDIR / "held_step_dashboard_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    return summary


def main():
    args = parse_args()
    data = build_model_and_data()
    cache = prepare_cache(
        data, OUTDIR / "held_step_dashboard_cache.npz", args.force_recompute
    )
    print(json.dumps(render(args, data, cache), indent=2))


if __name__ == "__main__":
    main()
