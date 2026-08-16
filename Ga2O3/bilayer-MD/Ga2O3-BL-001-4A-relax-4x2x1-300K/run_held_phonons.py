#!/usr/bin/env python3
"""Experimental multi-species HELD diagnostic for the 4 Angstrom run."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import numpy as np
from ase.data import atomic_masses, atomic_numbers

ROOT = Path(__file__).resolve().parent
ADAPTER_DIR = (ROOT / "../Ga2O3-BL-001-300K").resolve()
PARSER_DIR = (ROOT / "../Ga2O3-BL-001-4x2x1-300K").resolve()
sys.path.insert(0, str(PARSER_DIR))
sys.path.insert(0, str(ADAPTER_DIR))

from run_held_phonons import MultiSpeciesHeldModel, reciprocal_path, write_coefficients  # noqa: E402
from qe_supercell_to_tdep import parse_qe_input, parse_qe_output  # noqa: E402
from prepare_tdep import parse_initial_supercell, primitive_from_replication  # noqa: E402


def main():
    output = ROOT / "Ga2O3-BL-001-4x2x1.md.out"
    unit_cell, _symbols, _positions = parse_qe_input(ROOT / "Ga2O3-BL-001.relax.in")
    super_cell = [
        [value * 4 for value in unit_cell[0]],
        [value * 2 for value in unit_cell[1]],
        list(unit_cell[2]),
    ]
    super_symbols, super_reference = parse_initial_supercell(output)
    unit_symbols, unit_reference, _error = primitive_from_replication(
        super_symbols, super_reference
    )
    positions, forces, _epot, _ekin, temperatures, timestep, excluded = parse_qe_output(
        output, len(super_symbols)
    )
    unit_cell = np.asarray(unit_cell)
    super_cell = np.asarray(super_cell)
    unit_reference = np.asarray(unit_reference)
    super_reference = np.asarray(super_reference)
    positions = np.asarray(positions)
    forces = np.asarray(forces)
    numbers = np.asarray([atomic_numbers[symbol] for symbol in unit_symbols])
    masses = np.asarray([atomic_masses[number] for number in numbers])

    model = MultiSpeciesHeldModel(
        primitive_cell=unit_cell,
        basis_frac=unit_reference,
        supercell_cell=super_cell,
        ideal_supercell_frac=super_reference,
        primitive_numbers=numbers,
        primitive_masses=masses,
        cutoff_ang=5.5,
    )
    designs = [
        model.build_design_matrix(model.frame_displacements_cart(frame))
        for frame in positions
    ]
    step_coefficients, step_ranks = [], []
    for design, frame_forces in zip(designs, forces):
        values, _residuals, rank, _singular = np.linalg.lstsq(
            design, frame_forces.reshape(-1), rcond=None
        )
        step_coefficients.append(values)
        step_ranks.append(int(rank))
    step_coefficients = np.asarray(step_coefficients)
    global_design = np.vstack(designs)
    global_forces = forces.reshape(-1)
    global_coefficients, _residuals, global_rank, singular_values = np.linalg.lstsq(
        global_design, global_forces, rcond=None
    )
    per_frame_full_rank = min(step_ranks) == step_coefficients.shape[1]
    coefficients = (
        step_coefficients.mean(axis=0)
        if per_frame_full_rank
        else global_coefficients
    )
    aggregate = "per_frame_mean" if per_frame_full_rank else "global_least_squares"
    residual = global_design @ coefficients - global_forces
    q_path, x_values, tick_labels, tick_positions = reciprocal_path(unit_cell)
    frequencies = model.dispersion_thz_from_reduced_path(coefficients, q_path)

    outdir = ROOT / "HELD_300K"
    outdir.mkdir(parents=True, exist_ok=True)
    write_coefficients(outdir / "held_coefficients.csv", model.basis_labels, coefficients)
    np.savetxt(
        outdir / "held_step_coefficients.csv",
        step_coefficients,
        delimiter=",",
        header=",".join(model.basis_labels),
        comments="",
    )
    np.savetxt(
        outdir / "held_dispersion.dat",
        np.column_stack([x_values, frequencies]),
        fmt="%.12e",
        header="x_Ainv " + " ".join(
            f"held_b{index + 1}_THz" for index in range(frequencies.shape[1])
        ),
    )
    summary = {
        "status": "experimental_multispecies_adapter",
        "n_unit_atoms": len(unit_symbols),
        "n_supercell_atoms": len(super_symbols),
        "n_frames": len(positions),
        "excluded_incomplete_frames": excluded,
        "timestep_fs": timestep,
        "trajectory_time_fs": len(positions) * timestep,
        "mean_temperature_K": float(np.mean(temperatures)),
        "latest_complete_temperature_K": float(temperatures[-1]),
        "cutoff_angstrom": 5.5,
        "aggregate": aggregate,
        "n_symmetry_operations": int(len(model._primitive_symmetry()[0])),
        "n_offsite_pairs": int(len(model.offsite_keys)),
        "n_coefficients": int(len(coefficients)),
        "equations_per_frame": int(designs[0].shape[0]),
        "minimum_per_frame_rank": min(step_ranks),
        "maximum_per_frame_rank": max(step_ranks),
        "global_equations": int(global_design.shape[0]),
        "global_design_rank": int(global_rank),
        "global_condition_number": float(singular_values[0] / singular_values[-1]),
        "force_component_rmse_eV_per_A": float(np.sqrt(np.mean(residual**2))),
        "force_vector_rmse_eV_per_A": float(
            np.sqrt(np.mean(np.sum(residual.reshape(-1, 3) ** 2, axis=1)))
        ),
        "minimum_frequency_THz": float(np.min(frequencies)),
        "maximum_frequency_THz": float(np.max(frequencies)),
        "negative_frequency_values": int(np.sum(frequencies < -1.0e-6)),
        "total_frequency_values": int(frequencies.size),
    }
    (outdir / "held_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10.5, 6.6), constrained_layout=True)
    for position in tick_positions:
        ax.axvline(position, color="#aaaaaa", linewidth=0.7, zorder=0)
    for band in range(frequencies.shape[1]):
        ax.plot(x_values, frequencies[:, band], color="#b13a2f", linewidth=1.0)
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.set_xlim(x_values[0], x_values[-1])
    ax.set_xticks(tick_positions, tick_labels)
    ax.set_ylabel("Frequency (THz)")
    ax.set_title(
        f"Ga₂O₃ 4 Å bilayer · experimental HELD ({len(positions)} MD frames)"
    )
    ax.grid(axis="y", alpha=0.2)
    fig.savefig(outdir / "Ga2O3-BL-001-4A-HELD-experimental.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
