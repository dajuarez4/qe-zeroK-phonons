#!/usr/bin/env python3
"""Collect all Al69 finite-displacement QE calculations into one compressed NPZ.

The archive contains fixed-size numerical arrays and never modifies the QE files.
Missing or unfinished outputs are represented by NaN values and status flags, so
the same command can be rerun while the Slurm array is still progressing.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


RY_TO_EV = 13.605693122994
BOHR_TO_ANGSTROM = 0.529177210903
RY_BOHR_TO_EV_ANGSTROM = RY_TO_EV / BOHR_TO_ANGSTROM
KBAR_TO_GPA = 0.1

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
ENERGY_RE = re.compile(rf"^\s*!\s+total energy\s+=\s+({FLOAT})\s+Ry", re.MULTILINE)
FORCE_RE = re.compile(
    rf"atom\s+\d+\s+type\s+\d+\s+force\s+=\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})"
)
STRESS_HEADER_RE = re.compile(r"total\s+stress\s+\(Ry/bohr\*\*3\).*?\(kbar\)", re.I)
NAT_RE = re.compile(r"\bnat\s*=\s*(\d+)", re.I)


def number(value: str) -> float:
    """Convert QE's E/D exponent notation to a Python float."""
    return float(value.replace("D", "E").replace("d", "e"))


def parse_qe_input(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return symbols, fractional positions, and cell (bohr) from a QE input."""
    text = path.read_text(errors="replace")
    nat_match = NAT_RE.search(text)
    if not nat_match:
        raise ValueError(f"nat was not found in {path}")
    nat = int(nat_match.group(1))
    lines = text.splitlines()

    cell = None
    symbols = None
    positions = None
    for i, line in enumerate(lines):
        words = line.split()
        if words and words[0].upper() == "CELL_PARAMETERS":
            unit = words[1].lower().strip("{}()") if len(words) > 1 else "alat"
            values = np.array([[number(x) for x in lines[i + j].split()[:3]] for j in range(1, 4)])
            if unit == "angstrom":
                values /= BOHR_TO_ANGSTROM
            elif unit != "bohr":
                raise ValueError(f"unsupported CELL_PARAMETERS unit {unit!r} in {path}")
            cell = values
        elif words and words[0].upper() == "ATOMIC_POSITIONS":
            unit = words[1].lower().strip("{}()") if len(words) > 1 else "alat"
            rows = [lines[i + j].split() for j in range(1, nat + 1)]
            symbols = np.array([row[0] for row in rows], dtype="U4")
            values = np.array([[number(x) for x in row[1:4]] for row in rows])
            if unit == "crystal":
                positions = values
            elif unit == "angstrom":
                if cell is None:
                    raise ValueError(f"cell must precede Cartesian positions in {path}")
                positions = (values / BOHR_TO_ANGSTROM) @ np.linalg.inv(cell)
            elif unit == "bohr":
                if cell is None:
                    raise ValueError(f"cell must precede Cartesian positions in {path}")
                positions = values @ np.linalg.inv(cell)
            else:
                raise ValueError(f"unsupported ATOMIC_POSITIONS unit {unit!r} in {path}")

    if cell is None or symbols is None or positions is None:
        raise ValueError(f"cell or atomic positions were not found in {path}")
    return symbols, positions, cell


def parse_qe_output(path: Path, nat: int) -> dict[str, object]:
    """Extract the final SCF energy, force block, stress, and completion state."""
    text = path.read_text(errors="replace")
    energies = [number(x) for x in ENERGY_RE.findall(text)]
    force_rows = FORCE_RE.findall(text)
    forces = np.full((nat, 3), np.nan)
    if len(force_rows) >= nat:
        forces = np.array([[number(x) for x in row] for row in force_rows[-nat:]])

    stress = np.full((3, 3), np.nan)
    headers = list(STRESS_HEADER_RE.finditer(text))
    if headers:
        tail = text[headers[-1].end() :].splitlines()
        rows = []
        for line in tail:
            vals = re.findall(FLOAT, line)
            if len(vals) >= 6:
                rows.append([number(x) for x in vals[-3:]])
                if len(rows) == 3:
                    stress = np.array(rows)
                    break

    return {
        "energy_ry": energies[-1] if energies else np.nan,
        "forces_ry_bohr": forces,
        "stress_kbar": stress,
        "job_done": "JOB DONE." in text,
        "scf_converged": "convergence has been achieved" in text,
        "output_size_bytes": path.stat().st_size,
    }


def build_archive(root: Path, destination: Path, expected: int) -> None:
    simulation_ids = np.arange(1, expected + 1, dtype=np.int32)
    input_paths = [root / "displacements" / f"disp-{i:03d}" / "Al69_Ga2O3.fd.scf.in" for i in simulation_ids]
    output_paths = [p.with_suffix(".out") for p in input_paths]

    reference = next((p for p in input_paths if p.is_file()), None)
    if reference is None:
        raise FileNotFoundError(f"no displacement inputs found below {root}")
    ref_symbols, _, _ = parse_qe_input(reference)
    nat = len(ref_symbols)

    symbols = np.full((expected, nat), "", dtype="U4")
    positions = np.full((expected, nat, 3), np.nan)
    cells = np.full((expected, 3, 3), np.nan)
    energies = np.full(expected, np.nan)
    forces = np.full((expected, nat, 3), np.nan)
    stresses = np.full((expected, 3, 3), np.nan)
    input_exists = np.zeros(expected, dtype=bool)
    output_exists = np.zeros(expected, dtype=bool)
    job_done = np.zeros(expected, dtype=bool)
    scf_converged = np.zeros(expected, dtype=bool)
    output_sizes = np.zeros(expected, dtype=np.int64)
    errors = np.full(expected, "", dtype="U512")

    for index, (input_path, output_path) in enumerate(zip(input_paths, output_paths)):
        try:
            if input_path.is_file():
                input_exists[index] = True
                symbols[index], positions[index], cells[index] = parse_qe_input(input_path)
            if output_path.is_file():
                output_exists[index] = True
                result = parse_qe_output(output_path, nat)
                energies[index] = result["energy_ry"]
                forces[index] = result["forces_ry_bohr"]
                stresses[index] = result["stress_kbar"]
                job_done[index] = result["job_done"]
                scf_converged[index] = result["scf_converged"]
                output_sizes[index] = result["output_size_bytes"]
        except (OSError, ValueError) as exc:
            errors[index] = str(exc)[:512]

    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        destination,
        format_version=np.array("1.0"),
        description=np.array("Al69_0K Quantum ESPRESSO finite-displacement calculations"),
        simulation_id=simulation_ids,
        relative_directory=np.array([f"displacements/disp-{i:03d}" for i in simulation_ids]),
        symbols=symbols,
        positions_fractional=positions,
        cell_bohr=cells,
        cell_angstrom=cells * BOHR_TO_ANGSTROM,
        energy_ry=energies,
        energy_ev=energies * RY_TO_EV,
        forces_ry_bohr=forces,
        forces_ev_angstrom=forces * RY_BOHR_TO_EV_ANGSTROM,
        stress_kbar=stresses,
        stress_gpa=stresses * KBAR_TO_GPA,
        input_exists=input_exists,
        output_exists=output_exists,
        job_done=job_done,
        scf_converged=scf_converged,
        output_size_bytes=output_sizes,
        parse_error=errors,
    )

    print(f"Archive: {destination}")
    print(f"Simulations expected: {expected}")
    print(f"Inputs found: {input_exists.sum()}")
    print(f"Outputs found: {output_exists.sum()}")
    print(f"Jobs completed: {job_done.sum()}")
    print(f"Rows with parse errors: {np.count_nonzero(errors != '')}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parent,
        help="Al69_0K_phonopy directory (default: script directory)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="output NPZ (default: ROOT/Al69_0K_all_simulations.npz)",
    )
    parser.add_argument("--expected", type=int, default=240, help="expected number of displacements")
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve() if args.output else root / "Al69_0K_all_simulations.npz"
    build_archive(root, output, args.expected)


if __name__ == "__main__":
    main()
