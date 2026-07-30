#!/usr/bin/env python3
"""Convert this fixed-cell QE pw.x MD output to TDEP text input files."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
RY_BOHR_TO_EV_ANG = RY_TO_EV / BOHR_TO_ANG


def card(lines: list[str], name: str, count: int) -> list[list[str]]:
    for index, line in enumerate(lines):
        if line.strip().upper().startswith(name):
            result = []
            for candidate in lines[index + 1 :]:
                fields = candidate.split()
                if len(fields) < 3:
                    continue
                result.append(fields)
                if len(result) == count:
                    return result
    raise ValueError(f"Could not read {count} rows after {name}")


def parse_qe_input(path: Path):
    lines = path.read_text().splitlines()
    cell = [[float(x) for x in row[:3]] for row in card(lines, "CELL_PARAMETERS", 3)]
    atoms = card(lines, "ATOMIC_POSITIONS", 10)
    symbols = [row[0] for row in atoms]
    positions = [[float(x) for x in row[1:4]] for row in atoms]
    return cell, symbols, positions


def parse_output(path: Path, nat: int):
    text = path.read_text(errors="replace")

    position_blocks = []
    pattern = re.compile(
        r"^ATOMIC_POSITIONS\s+\(crystal\)\s*$((?:\n\s*[A-Za-z]+\s+"
        r"[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+[^\n]*){" + str(nat) + r"})",
        re.MULTILINE,
    )
    for match in pattern.finditer(text):
        rows = []
        for line in match.group(1).strip().splitlines():
            fields = line.split()
            rows.append([float(x) for x in fields[1:4]])
        position_blocks.append(rows)

    force_blocks = []
    for marker in re.finditer(r"^\s*Forces acting on atoms .*:$", text, re.MULTILINE):
        rows = []
        tail = text[marker.end() :]
        for match in re.finditer(
            r"^\s*atom\s+\d+\s+type\s+\d+\s+force\s*=\s*"
            r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
            tail,
            re.MULTILINE,
        ):
            rows.append([float(x) * RY_BOHR_TO_EV_ANG for x in match.groups()])
            if len(rows) == nat:
                break
        if len(rows) == nat:
            force_blocks.append(rows)

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
    timestep = re.search(
        r"Time step\s*=\s*[-+0-9.Ee]+\s+a\.u\.,\s*([-+0-9.Ee]+)\s+femto-seconds",
        text,
    )
    if not timestep:
        raise ValueError("Could not parse the MD timestep")

    # force_blocks[0]/energies[0] are the pre-MD evaluation. Each subsequent
    # force block belongs to the preceding printed MD position block.
    nframes = min(
        len(position_blocks),
        max(0, len(force_blocks) - 1),
        max(0, len(energies) - 1),
        len(kinetic),
        len(temperatures),
    )
    if nframes == 0:
        raise ValueError("No complete QE MD position/force pairs were found")

    return (
        position_blocks[:nframes],
        force_blocks[1 : nframes + 1],
        energies[1 : nframes + 1],
        kinetic[:nframes],
        temperatures[:nframes],
        float(timestep.group(1)),
        len(position_blocks) - nframes,
    )


def write_poscar(
    path: Path,
    title: str,
    cell: list[list[float]],
    symbols: list[str],
    positions: list[list[float]],
    order: list[int],
) -> None:
    species = list(dict.fromkeys(symbols))
    lines = [title, "1.0"]
    lines.extend(" ".join(f"{value: .16f}" for value in row) for row in cell)
    lines.append(" ".join(species))
    lines.append(" ".join(str(symbols.count(element)) for element in species))
    lines.append("Direct")
    for index in order:
        row = positions[index]
        lines.append(" ".join(f"{value: .16f}" for value in row))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("qe_input", type=Path)
    parser.add_argument("qe_output", type=Path)
    parser.add_argument("-o", "--outdir", type=Path, default=Path("TDEP_300K"))
    args = parser.parse_args()

    cell, symbols, ideal_positions = parse_qe_input(args.qe_input)
    (
        positions,
        forces,
        potential_energies,
        kinetic_energies,
        temperatures,
        timestep_fs,
        incomplete_positions,
    ) = parse_output(args.qe_output, len(symbols))

    species = list(dict.fromkeys(symbols))
    order = [i for element in species for i, symbol in enumerate(symbols) if symbol == element]
    args.outdir.mkdir(parents=True, exist_ok=True)

    for name, title in (
        ("infile.ucposcar", "Ga2O3 (001) bilayer unit cell"),
        ("infile.ssposcar", "Ga2O3 (001) bilayer MD simulation cell"),
    ):
        write_poscar(args.outdir / name, title, cell, symbols, ideal_positions, order)

    with (args.outdir / "infile.positions").open("w") as handle:
        for frame in positions:
            for index in order:
                handle.write(" ".join(f"{value: .16e}" for value in frame[index]) + "\n")

    with (args.outdir / "infile.forces").open("w") as handle:
        for frame in forces:
            for index in order:
                handle.write(" ".join(f"{value: .16e}" for value in frame[index]) + "\n")

    mean_temperature = sum(temperatures) / len(temperatures)
    (args.outdir / "infile.meta").write_text(
        f"{len(symbols)}\n{len(positions)}\n{timestep_fs:.10f}\n{mean_temperature:.10f}\n"
    )

    with (args.outdir / "infile.stat").open("w") as handle:
        for i, (epot, ekin, temperature) in enumerate(
            zip(potential_energies, kinetic_energies, temperatures), start=1
        ):
            etot = epot + ekin
            handle.write(
                f"{i:d} {i * timestep_fs:.10f} {etot:.16e} {epot:.16e} "
                f"{ekin:.16e} {temperature:.10f} 0 0 0 0 0 0 0\n"
            )

    (args.outdir / "infile.qpoints_dispersion").write_text(
        "CUSTOM\n"
        "80\n"
        "4\n"
        "0.0 0.0 0.0  0.5 0.0 0.0  GM X\n"
        "0.5 0.0 0.0  0.5 0.5 0.0  X S\n"
        "0.5 0.5 0.0  0.0 0.5 0.0  S Y\n"
        "0.0 0.5 0.0  0.0 0.0 0.0  Y GM\n"
    )

    if not math.isfinite(mean_temperature):
        raise ValueError("Non-finite mean temperature")
    print(f"Wrote {len(positions)} complete configurations to {args.outdir}")
    print(f"Timestep: {timestep_fs:.4f} fs")
    print(f"Mean sampled temperature: {mean_temperature:.3f} K")
    print(f"Unpaired/incomplete printed position blocks excluded: {incomplete_positions}")


if __name__ == "__main__":
    main()
