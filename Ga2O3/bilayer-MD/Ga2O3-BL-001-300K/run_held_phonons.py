#!/usr/bin/env python3
"""Experimental multi-species HELD fit for the Ga2O3 (001) trajectory."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import numpy as np
import spglib
from ase.data import atomic_masses, atomic_numbers

HELD_ROOT = Path(
    os.environ.get("HELD_ROOT", "/Users/dajuarez4/Documents/Fe/HELD")
).resolve()
if not (HELD_ROOT / "held" / "model.py").is_file():
    raise FileNotFoundError(f"Could not find the HELD package under {HELD_ROOT}")
sys.path.insert(0, str(HELD_ROOT))

from held.model import SymmetryHeldModel  # noqa: E402

from qe_md_to_tdep import parse_output, parse_qe_input  # noqa: E402


class MultiSpeciesHeldModel(SymmetryHeldModel):
    """Minimal species/mass adapter around HELD's symmetry model."""

    def __init__(
        self,
        primitive_cell,
        basis_frac,
        supercell_cell,
        ideal_supercell_frac,
        primitive_numbers,
        primitive_masses,
        cutoff_ang,
        symprec=1.0e-5,
    ):
        self.primitive_numbers = np.asarray(primitive_numbers, dtype=int)
        self.primitive_masses = np.asarray(primitive_masses, dtype=float)
        if len(self.primitive_numbers) != len(basis_frac):
            raise ValueError("Primitive chemical-number count does not match the basis")
        super().__init__(
            primitive_cell=primitive_cell,
            basis_frac=basis_frac,
            supercell_cell=supercell_cell,
            ideal_supercell_frac=ideal_supercell_frac,
            atomic_number=int(self.primitive_numbers[0]),
            mass_amu=1.0,
            num_shells=1,
            cutoff_ang=cutoff_ang,
            symprec=symprec,
        )

    def _primitive_symmetry(self):
        symmetry = spglib.get_symmetry(
            (self.uc_cell, self.uc_frac, self.primitive_numbers.tolist()),
            symprec=self.symprec,
        )
        if symmetry is None:
            raise ValueError("spglib could not determine the species-aware symmetry")
        return (
            np.asarray(symmetry["rotations"], dtype=int),
            np.asarray(symmetry["translations"], dtype=float),
        )

    def _source_basis_index(self, frac_coord):
        for index, basis_frac in enumerate(self.uc_frac):
            difference = np.asarray(frac_coord) - basis_frac
            if np.allclose(difference - np.rint(difference), 0.0, atol=self.mapping_tol):
                return index
        raise ValueError(f"Could not identify source basis atom for {frac_coord}")

    def _map_primitive_atom(self, frac_coord, rotation, translation):
        source_index = self._source_basis_index(frac_coord)
        source_number = int(self.primitive_numbers[source_index])
        transformed = rotation @ frac_coord + translation
        for atom_j, basis_frac in enumerate(self.uc_frac):
            if int(self.primitive_numbers[atom_j]) != source_number:
                continue
            difference = transformed - basis_frac
            shift = np.rint(difference).astype(int)
            if np.allclose(difference - shift, 0.0, atol=self.mapping_tol):
                return atom_j, shift
        raise ValueError(
            f"Could not map species Z={source_number} at {np.asarray(frac_coord).tolist()}"
        )

    def _mass_factor(self, atom_i, atom_j):
        return float(
            np.sqrt(self.primitive_masses[int(atom_i)] * self.primitive_masses[int(atom_j)])
        )


def reciprocal_path(cell, points_per_segment=80):
    labels = ["Γ", "X", "S", "Y", "Γ"]
    special = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.0, 0.0]),
        np.array([0.5, 0.5, 0.0]),
        np.array([0.0, 0.5, 0.0]),
        np.array([0.0, 0.0, 0.0]),
    ]
    reciprocal = 2.0 * np.pi * np.linalg.inv(cell).T
    q_blocks = []
    tick_positions = [0.0]
    cumulative = 0.0
    for segment_index, (start, end) in enumerate(zip(special[:-1], special[1:])):
        fractions = np.linspace(0.0, 1.0, points_per_segment + 1)
        if segment_index:
            fractions = fractions[1:]
        block = start[None, :] + fractions[:, None] * (end - start)[None, :]
        q_blocks.append(block)
        segment_cart = block @ reciprocal.T
        cumulative += float(np.linalg.norm(np.diff(segment_cart, axis=0), axis=1).sum())
        tick_positions.append(cumulative)
    q_path = np.vstack(q_blocks)
    q_cart = q_path @ reciprocal.T
    x_values = np.zeros(len(q_path))
    x_values[1:] = np.cumsum(np.linalg.norm(np.diff(q_cart, axis=0), axis=1))
    return q_path, x_values, labels, np.asarray(tick_positions)


def write_coefficients(path, labels, values):
    path.write_text(
        ",".join(labels) + "\n" + ",".join(f"{value:.16e}" for value in values) + "\n"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qe-input", type=Path, default=Path("Ga2O3-BL-001.md.in"))
    parser.add_argument("--qe-output", type=Path, default=Path("Ga2O3-BL-001.md.out"))
    parser.add_argument("--outdir", type=Path, default=Path("HELD_300K"))
    parser.add_argument("--cutoff", type=float, default=2.5)
    parser.add_argument(
        "--aggregate",
        choices=("mean", "global"),
        default="mean",
        help="HELD coefficient aggregation. The official CLI defaults to per-frame mean.",
    )
    args = parser.parse_args()

    cell, symbols, ideal_positions = parse_qe_input(args.qe_input)
    (
        positions,
        forces,
        _potential,
        _kinetic,
        temperatures,
        timestep_fs,
        excluded,
    ) = parse_output(args.qe_output, len(symbols))

    cell = np.asarray(cell, dtype=float)
    ideal_positions = np.asarray(ideal_positions, dtype=float)
    positions = np.asarray(positions, dtype=float)
    forces = np.asarray(forces, dtype=float)
    numbers = np.asarray([atomic_numbers[symbol] for symbol in symbols], dtype=int)
    masses = np.asarray([atomic_masses[number] for number in numbers], dtype=float)

    model = MultiSpeciesHeldModel(
        primitive_cell=cell,
        basis_frac=ideal_positions,
        supercell_cell=cell,
        ideal_supercell_frac=ideal_positions,
        primitive_numbers=numbers,
        primitive_masses=masses,
        cutoff_ang=args.cutoff,
    )

    design_blocks = [
        model.build_design_matrix(model.frame_displacements_cart(frame))
        for frame in positions
    ]
    global_design = np.vstack(design_blocks)
    global_forces = forces.reshape(-1)
    global_coefficients, _residuals, global_rank, singular_values = np.linalg.lstsq(
        global_design, global_forces, rcond=None
    )
    step_coefficients = []
    step_ranks = []
    for design, frame_forces in zip(design_blocks, forces):
        values, _residuals, frame_rank, _singular_values = np.linalg.lstsq(
            design, frame_forces.reshape(-1), rcond=None
        )
        step_coefficients.append(values)
        step_ranks.append(int(frame_rank))
    step_coefficients = np.asarray(step_coefficients)
    coefficients = (
        step_coefficients.mean(axis=0)
        if args.aggregate == "mean"
        else global_coefficients
    )
    predicted = global_design @ coefficients
    residual = predicted - global_forces
    force_rmse = float(np.sqrt(np.mean(residual**2)))
    force_vector_rmse = float(
        np.sqrt(np.mean(np.sum(residual.reshape(-1, 3) ** 2, axis=1)))
    )

    q_path, x_values, tick_labels, tick_positions = reciprocal_path(cell)
    frequencies = model.dispersion_thz_from_reduced_path(coefficients, q_path)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_coefficients(
        args.outdir / "held_coefficients.csv", model.basis_labels, coefficients
    )
    np.savetxt(
        args.outdir / "held_step_coefficients.csv",
        step_coefficients,
        delimiter=",",
        header=",".join(model.basis_labels),
        comments="",
    )
    np.savetxt(
        args.outdir / "held_dispersion.dat",
        np.column_stack([x_values, frequencies]),
        fmt="%.12e",
        header="x_Ainv " + " ".join(f"held_b{i + 1}_THz" for i in range(frequencies.shape[1])),
    )

    summary = {
        "status": "experimental_multispecies_adapter",
        "held_root": str(HELD_ROOT),
        "n_atoms": len(symbols),
        "symbols": symbols,
        "n_frames": len(positions),
        "excluded_incomplete_frames": excluded,
        "timestep_fs": timestep_fs,
        "trajectory_time_fs": len(positions) * timestep_fs,
        "mean_temperature_K": float(np.mean(temperatures)),
        "latest_complete_temperature_K": float(temperatures[-1]),
        "cutoff_angstrom": args.cutoff,
        "n_symmetry_operations": int(len(model._primitive_symmetry()[0])),
        "n_offsite_pairs": int(len(model.offsite_keys)),
        "n_coefficients": int(len(coefficients)),
        "n_equations": int(global_design.shape[0]),
        "aggregate": args.aggregate,
        "per_frame_equations": int(design_blocks[0].shape[0]),
        "minimum_per_frame_rank": int(min(step_ranks)),
        "maximum_per_frame_rank": int(max(step_ranks)),
        "global_design_rank": int(global_rank),
        "global_condition_number": float(
            singular_values[0] / singular_values[-1]
            if singular_values.size and singular_values[-1] > 0
            else np.inf
        ),
        "force_component_rmse_eV_per_A": force_rmse,
        "force_vector_rmse_eV_per_A": force_vector_rmse,
        "minimum_frequency_THz": float(np.min(frequencies)),
        "maximum_frequency_THz": float(np.max(frequencies)),
        "negative_frequency_values": int(np.sum(frequencies < -1.0e-6)),
        "total_frequency_values": int(frequencies.size),
    }
    (args.outdir / "held_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

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
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_ylabel("Frequency (THz)")
    ax.set_title(
        f"Ga₂O₃ (001) experimental HELD adapter ({len(positions)} MD frames)"
    )
    ax.grid(axis="y", alpha=0.2)
    fig.savefig(args.outdir / "Ga2O3-BL-001-HELD-experimental.png", dpi=180)
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
