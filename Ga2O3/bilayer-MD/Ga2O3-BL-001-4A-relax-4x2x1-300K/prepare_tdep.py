#!/usr/bin/env python3
"""Prepare TDEP files from the interrupted nine-frame 4 Angstrom QE run."""

from __future__ import annotations

import json
import math
import os
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ADAPTER_DIR = (ROOT / "../Ga2O3-BL-001-4x2x1-300K").resolve()
sys.path.insert(0, str(ADAPTER_DIR))

from qe_supercell_to_tdep import (  # noqa: E402
    parse_qe_input,
    parse_qe_output,
    species_order,
    write_poscar,
)


def parse_initial_supercell(path: Path, nat: int = 80):
    """Read the actual initial crystal coordinates printed by QE."""
    text = path.read_text(errors="replace")
    marker = text.find("Crystallographic axes")
    if marker < 0:
        raise ValueError("Could not find QE's initial crystallographic axes")
    atom_pattern = re.compile(
        r"^\s*\d+\s+([A-Za-z]+)\s+tau\(\s*\d+\)\s*=\s*\(\s*"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*\)",
        re.MULTILINE,
    )
    matches = list(atom_pattern.finditer(text, marker))[:nat]
    if len(matches) != nat:
        raise ValueError(f"Expected {nat} initial atoms, found {len(matches)}")
    symbols = [match.group(1) for match in matches]
    positions = [[float(value) for value in match.groups()[1:]] for match in matches]
    return symbols, positions


def primitive_from_replication(super_symbols, super_positions, replication=(4, 2, 1)):
    """Recover the ten-atom primitive reference from QE's 4x2x1 ordering."""
    rx, ry, rz = replication
    nunit = len(super_symbols) // (rx * ry * rz)
    unit_symbols = super_symbols[:nunit]
    unit_positions = [
        [
            (super_positions[index][0] * rx) % 1.0,
            (super_positions[index][1] * ry) % 1.0,
            (super_positions[index][2] * rz) % 1.0,
        ]
        for index in range(nunit)
    ]
    maximum_error = 0.0
    index = 0
    for ix in range(rx):
        for iy in range(ry):
            for iz in range(rz):
                for symbol, position in zip(unit_symbols, unit_positions):
                    expected = [
                        (position[0] + ix) / rx,
                        (position[1] + iy) / ry,
                        (position[2] + iz) / rz,
                    ]
                    if super_symbols[index] != symbol:
                        raise ValueError(f"Species replication mismatch at atom {index + 1}")
                    error = max(
                        abs(super_positions[index][axis] - expected[axis])
                        for axis in range(3)
                    )
                    maximum_error = max(maximum_error, error)
                    index += 1
    if maximum_error > 2.0e-6:
        raise ValueError(f"Initial structure is not a 4x2x1 replication: {maximum_error}")
    return unit_symbols, unit_positions, maximum_error


def prepare(outdir: Path):
    output = ROOT / "Ga2O3-BL-001-4x2x1.md.out"
    relax_input = ROOT / "Ga2O3-BL-001.relax.in"
    unit_cell, _input_symbols, _input_positions = parse_qe_input(relax_input)
    super_cell = [
        [value * 4 for value in unit_cell[0]],
        [value * 2 for value in unit_cell[1]],
        list(unit_cell[2]),
    ]
    super_symbols, super_reference = parse_initial_supercell(output)
    unit_symbols, unit_reference, replication_error = primitive_from_replication(
        super_symbols, super_reference
    )
    (
        positions,
        forces,
        potential,
        kinetic,
        temperatures,
        timestep,
        excluded,
    ) = parse_qe_output(output, len(super_symbols))
    _species, order = species_order(super_symbols)

    outdir.mkdir(parents=True, exist_ok=True)
    write_poscar(
        outdir / "infile.ucposcar",
        "Ga2O3 (001) 4 Angstrom relaxed bilayer reference",
        unit_cell,
        unit_symbols,
        unit_reference,
    )
    write_poscar(
        outdir / "infile.ssposcar",
        "Ga2O3 (001) 4 Angstrom 4x2x1 MD reference",
        super_cell,
        super_symbols,
        super_reference,
    )
    with (outdir / "infile.positions").open("w") as handle:
        for frame in positions:
            for atom in order:
                handle.write(" ".join(f"{value: .16e}" for value in frame[atom]) + "\n")
    with (outdir / "infile.forces").open("w") as handle:
        for frame in forces:
            for atom in order:
                handle.write(" ".join(f"{value: .16e}" for value in frame[atom]) + "\n")

    mean_temperature = sum(temperatures) / len(temperatures)
    if not math.isfinite(mean_temperature):
        raise ValueError("Non-finite mean temperature")
    (outdir / "infile.meta").write_text(
        f"{len(super_symbols)}\n{len(positions)}\n{timestep:.10f}\n"
        f"{mean_temperature:.10f}\n"
    )
    with (outdir / "infile.stat").open("w") as handle:
        for index, (epot, ekin, temperature) in enumerate(
            zip(potential, kinetic, temperatures), start=1
        ):
            handle.write(
                f"{index:d} {index * timestep:.10f} {epot + ekin:.16e} "
                f"{epot:.16e} {ekin:.16e} {temperature:.10f} 0 0 0 0 0 0 0\n"
            )
    (outdir / "infile.qpoints_dispersion").write_text(
        "CUSTOM\n80\n4\n"
        "0.0 0.0 0.0  0.5 0.0 0.0  GM X\n"
        "0.5 0.0 0.0  0.5 0.5 0.0  X S\n"
        "0.5 0.5 0.0  0.0 0.5 0.0  S Y\n"
        "0.0 0.5 0.0  0.0 0.0 0.0  Y GM\n"
    )
    summary = {
        "n_unit_atoms": len(unit_symbols),
        "n_supercell_atoms": len(super_symbols),
        "replication": [4, 2, 1],
        "reference_replication_max_error_fractional": replication_error,
        "n_frames": len(positions),
        "excluded_incomplete_position_blocks": excluded,
        "timestep_fs": timestep,
        "trajectory_time_fs": len(positions) * timestep,
        "mean_temperature_K": mean_temperature,
        "latest_complete_temperature_K": temperatures[-1],
        "qe_job_completed": "JOB DONE." in output.read_text(errors="replace"),
    }
    (outdir / "reference_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    prepare(ROOT / os.environ.get("TDEP_OUTDIR", "TDEP_300K"))
