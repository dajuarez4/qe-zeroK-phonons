#!/usr/bin/env python3
"""Convert the complete QE MD frames to the NPZ format used by No-Vito."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

from qe_supercell_to_tdep import parse_qe_input, parse_qe_output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=Path)
    parser.add_argument("output_file", type=Path)
    parser.add_argument("-o", "--npz", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument(
        "--md-only",
        action="store_true",
        help="Write exactly the complete MD frames, without prepending the input geometry.",
    )
    args = parser.parse_args()

    cell, species, reference = parse_qe_input(args.input_file)
    positions, _forces, _energies, _kinetic, temperatures, timestep_fs, excluded = (
        parse_qe_output(args.output_file, len(species))
    )
    reference = np.asarray(reference, dtype=float)
    positions = np.asarray(positions, dtype=float)

    output_text = args.output_file.read_text(errors="replace")
    pressures_kbar = [
        float(value)
        for value in re.findall(
            r"^\s*total\s+stress.*P=\s*([-+0-9.Ee]+)",
            output_text,
            re.MULTILINE,
        )
    ]
    md_pressures_gpa = np.full(len(positions), np.nan)
    if len(pressures_kbar) >= len(positions) + 1:
        md_pressures_gpa = 0.1 * np.asarray(
            pressures_kbar[1 : len(positions) + 1], dtype=float
        )

    if args.md_only:
        render_positions = positions
        iteration = np.arange(1, len(positions) + 1, dtype=int)
        time_ps = iteration.astype(float) * timestep_fs / 1000.0
        temperature = np.asarray(temperatures, dtype=float)
        pressure_gpa = md_pressures_gpa
        description = "Complete QE MD frames"
    else:
        # Prepend the QE input geometry so the animation shows displacement
        # from the reference before the complete MD frames.
        render_positions = np.concatenate((reference[None, :, :], positions), axis=0)
        iteration = np.arange(len(render_positions), dtype=int)
        time_ps = iteration.astype(float) * timestep_fs / 1000.0
        temperature = np.concatenate(
            ([np.nan], np.asarray(temperatures, dtype=float))
        )
        initial_pressure = (
            0.1 * pressures_kbar[0] if pressures_kbar else np.nan
        )
        pressure_gpa = np.concatenate(([initial_pressure], md_pressures_gpa))
        description = "Reference QE input geometry followed by complete MD frames"
    nframes = render_positions.shape[0]

    metadata = {
        "source_input": args.input_file.name,
        "source_output": args.output_file.name,
        "description": description,
        "complete_md_frames": int(positions.shape[0]),
        "excluded_incomplete_frames": int(excluded),
        "timestep_fs": float(timestep_fs),
    }

    np.savez_compressed(
        args.npz,
        metadata_json=np.array(json.dumps(metadata)),
        species=np.asarray(species),
        positions=render_positions,
        positions_unit=np.asarray(["crystal"] * nframes),
        iteration=iteration,
        time_ps=time_ps,
        temperature_K=temperature,
        pressure_GPa=pressure_gpa,
        input_cell_parameters=np.asarray(cell, dtype=float),
        input_cell_unit=np.array("angstrom"),
    )

    frame_stats = []
    for frame_index, fractional in enumerate(positions, start=1):
        delta = fractional - reference
        delta -= np.rint(delta)
        cartesian_delta = delta @ cell
        magnitudes = np.linalg.norm(cartesian_delta, axis=1)
        frame_stats.append(
            {
                "md_frame": frame_index,
                "time_fs": frame_index * float(timestep_fs),
                "temperature_K": float(temperatures[frame_index - 1]),
                "rms_displacement_A": float(np.sqrt(np.mean(magnitudes**2))),
                "maximum_displacement_A": float(np.max(magnitudes)),
            }
        )

    summary = {
        **metadata,
        "rendered_structures": int(nframes),
        "displacements_relative_to_input_geometry": frame_stats,
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")

    print(f"Wrote {args.npz} with {nframes} structures")
    print(f"Wrote {args.summary}")


if __name__ == "__main__":
    main()
