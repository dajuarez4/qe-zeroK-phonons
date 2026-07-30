#!/usr/bin/env python3
"""Convert fixed-cell QE supercell MD output to TDEP input files."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
RY_BOHR_TO_EV_ANG = RY_TO_EV / BOHR_TO_ANG


def parse_qe_input(path: Path):
    text = path.read_text()
    nat_match = re.search(r"\bnat\s*=\s*(\d+)", text, re.IGNORECASE)
    if not nat_match:
        raise ValueError(f"Could not parse nat from {path}")
    nat = int(nat_match.group(1))
    lines = text.splitlines()

    cell_start = next(i for i, line in enumerate(lines) if line.startswith("CELL_PARAMETERS")) + 1
    cell = [
        [float(value) for value in lines[cell_start + offset].split()[:3]]
        for offset in range(3)
    ]

    pos_start = next(i for i, line in enumerate(lines) if line.startswith("ATOMIC_POSITIONS")) + 1
    symbols = []
    positions = []
    for line in lines[pos_start : pos_start + nat]:
        fields = line.split()
        symbols.append(fields[0])
        positions.append([float(value) for value in fields[1:4]])
    return cell, symbols, positions


def parse_qe_output(path: Path, nat: int):
    text = path.read_text(errors="replace")
    position_pattern = re.compile(
        r"^ATOMIC_POSITIONS\s+\(crystal\)\s*$((?:\n\s*[A-Za-z]+\s+"
        r"[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+[^\n]*){" + str(nat) + r"})",
        re.MULTILINE,
    )
    positions = []
    for match in position_pattern.finditer(text):
        frame = []
        for line in match.group(1).strip().splitlines():
            fields = line.split()
            frame.append([float(value) for value in fields[1:4]])
        positions.append(frame)

    forces = []
    for marker in re.finditer(r"^\s*Forces acting on atoms .*:$", text, re.MULTILINE):
        frame = []
        for match in re.finditer(
            r"^\s*atom\s+\d+\s+type\s+\d+\s+force\s*=\s*"
            r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
            text[marker.end() :],
            re.MULTILINE,
        ):
            frame.append(
                [float(value) * RY_BOHR_TO_EV_ANG for value in match.groups()]
            )
            if len(frame) == nat:
                break
        if len(frame) == nat:
            forces.append(frame)

    energies = [
        float(value) * RY_TO_EV
        for value in re.findall(
            r"^\s*!\s+total energy\s*=\s*([-+0-9.Ee]+)\s+Ry",
            text,
            re.MULTILINE,
        )
    ]
    kinetic = [
        float(value) * RY_TO_EV
        for value in re.findall(
            r"^\s*kinetic energy \(Ekin\)\s*=\s*([-+0-9.Ee]+)\s+Ry",
            text,
            re.MULTILINE,
        )
    ]
    temperatures = [
        float(value)
        for value in re.findall(
            r"^\s*temperature\s*=\s*([-+0-9.Ee]+)\s+K",
            text,
            re.MULTILINE,
        )
    ]
    timestep_match = re.search(
        r"Time step\s*=\s*[-+0-9.Ee]+\s+a\.u\.,\s*([-+0-9.Ee]+)\s+femto-seconds",
        text,
    )
    if not timestep_match:
        raise ValueError("Could not parse the QE timestep")

    nframes = min(
        len(positions),
        max(0, len(forces) - 1),
        max(0, len(energies) - 1),
        len(kinetic),
        len(temperatures),
    )
    if nframes == 0:
        raise ValueError("No complete position/force MD pairs were found")
    return (
        positions[:nframes],
        forces[1 : nframes + 1],
        energies[1 : nframes + 1],
        kinetic[:nframes],
        temperatures[:nframes],
        float(timestep_match.group(1)),
        len(positions) - nframes,
    )


def species_order(symbols):
    species = list(dict.fromkeys(symbols))
    return species, [
        index
        for element in species
        for index, symbol in enumerate(symbols)
        if symbol == element
    ]


def write_poscar(path, title, cell, symbols, positions):
    species, order = species_order(symbols)
    lines = [title, "1.0"]
    lines.extend(" ".join(f"{value: .16f}" for value in row) for row in cell)
    lines.append(" ".join(species))
    lines.append(" ".join(str(symbols.count(element)) for element in species))
    lines.append("Direct")
    lines.extend(
        " ".join(f"{value: .16f}" for value in positions[index])
        for index in order
    )
    path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--unit-input", type=Path, required=True)
    parser.add_argument("--supercell-input", type=Path, required=True)
    parser.add_argument("--qe-output", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("TDEP_300K"))
    args = parser.parse_args()

    unit_cell, unit_symbols, unit_positions = parse_qe_input(args.unit_input)
    super_cell, super_symbols, super_positions = parse_qe_input(args.supercell_input)
    (
        positions,
        forces,
        potential,
        kinetic,
        temperatures,
        timestep,
        excluded,
    ) = parse_qe_output(args.qe_output, len(super_symbols))
    _species, order = species_order(super_symbols)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_poscar(
        args.outdir / "infile.ucposcar",
        "Ga2O3 (001) bilayer 10-atom unit cell",
        unit_cell,
        unit_symbols,
        unit_positions,
    )
    write_poscar(
        args.outdir / "infile.ssposcar",
        "Ga2O3 (001) bilayer 4x2x1 MD supercell",
        super_cell,
        super_symbols,
        super_positions,
    )

    with (args.outdir / "infile.positions").open("w") as handle:
        for frame in positions:
            for index in order:
                handle.write(
                    " ".join(f"{value: .16e}" for value in frame[index]) + "\n"
                )
    with (args.outdir / "infile.forces").open("w") as handle:
        for frame in forces:
            for index in order:
                handle.write(
                    " ".join(f"{value: .16e}" for value in frame[index]) + "\n"
                )

    mean_temperature = sum(temperatures) / len(temperatures)
    (args.outdir / "infile.meta").write_text(
        f"{len(super_symbols)}\n{len(positions)}\n{timestep:.10f}\n"
        f"{mean_temperature:.10f}\n"
    )
    with (args.outdir / "infile.stat").open("w") as handle:
        for index, (epot, ekin, temperature) in enumerate(
            zip(potential, kinetic, temperatures), start=1
        ):
            handle.write(
                f"{index:d} {index * timestep:.10f} {epot + ekin:.16e} "
                f"{epot:.16e} {ekin:.16e} {temperature:.10f} 0 0 0 0 0 0 0\n"
            )
    (args.outdir / "infile.qpoints_dispersion").write_text(
        "CUSTOM\n80\n4\n"
        "0.0 0.0 0.0  0.5 0.0 0.0  GM X\n"
        "0.5 0.0 0.0  0.5 0.5 0.0  X S\n"
        "0.5 0.5 0.0  0.0 0.5 0.0  S Y\n"
        "0.0 0.5 0.0  0.0 0.0 0.0  Y GM\n"
    )
    if not math.isfinite(mean_temperature):
        raise ValueError("Non-finite mean temperature")
    print(f"Wrote {len(positions)} complete 80-atom configurations")
    print(f"Mean temperature: {mean_temperature:.3f} K")
    print(f"Incomplete position blocks excluded: {excluded}")


if __name__ == "__main__":
    main()
