#!/usr/bin/env python3
"""Create ready-to-run Phonopy displacement folders for Fe2O3 QHA volumes."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent
TEMPLATE = ROOT / "phonopy_V100_2x2x2"
EOS_ROOT = ROOT / "eos_birch_murnaghan"
VOLUMES = {
    "V102": EOS_ROOT / "alpha_Fe2O3.eosv102.out",
    "V104": EOS_ROOT / "alpha_Fe2O3.eosv104.out",
    "V106": EOS_ROOT / "alpha_Fe2O3.eosv106.out",
    "V108": EOS_ROOT / "alpha_Fe2O3.eosv108.out",
}


def final_positions(output: Path) -> list[str]:
    """Return the ten crystal-coordinate rows in QE's final-coordinate block."""
    text = output.read_text(errors="replace")
    if "JOB DONE" not in text:
        raise RuntimeError(f"EOS output is incomplete: {output}")
    blocks = re.findall(
        r"Begin final coordinates\s+ATOMIC_POSITIONS \(crystal\)\s+"
        r"(.*?)\s+End final coordinates",
        text,
        flags=re.DOTALL,
    )
    if not blocks:
        raise RuntimeError(f"No final crystal coordinates found in {output}")
    rows = [line.split() for line in blocks[-1].splitlines() if line.strip()]
    if len(rows) != 10 or any(len(row) < 4 for row in rows):
        raise RuntimeError(f"Expected ten final atomic positions in {output}")
    return [f"{r[0]:<4s}  {float(r[1]):.10f}  {float(r[2]):.10f}  {float(r[3]):.10f}" for r in rows]


def template_cell() -> list[list[float]]:
    text = (TEMPLATE / "unitcell.in").read_text()
    match = re.search(
        r"CELL_PARAMETERS angstrom\s+((?:[^\n]+\n){3})",
        text,
    )
    if not match:
        raise RuntimeError("Could not read the V100 template cell")
    return [list(map(float, line.split())) for line in match.group(1).splitlines()]


def determinant(cell: list[list[float]]) -> float:
    a, b, c = cell
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def make_unitcell(label: str, positions: list[str]) -> str:
    text = (TEMPLATE / "unitcell.in").read_text().replace("V100", label)
    ratio = int(label[1:]) / 100.0
    scale = ratio ** (1.0 / 3.0)
    cell = [[value * scale for value in row] for row in template_cell()]
    cell_text = "\n".join("  " + "  ".join(f"{x: .12f}" for x in row) for row in cell)
    text, ncell = re.subn(
        r"(?<=CELL_PARAMETERS angstrom\n)(?:[^\n]+\n){3}",
        cell_text + "\n",
        text,
    )
    position_text = "\n".join(positions)
    text, npos = re.subn(
        r"(?<=ATOMIC_POSITIONS crystal\n)(?:[^\n]+\n){10}",
        position_text + "\n",
        text,
    )
    if ncell != 1 or npos != 1:
        raise RuntimeError(f"Failed to update the unit-cell template for {label}")
    expected = 100.624948 * ratio
    actual = abs(determinant(cell))
    if not math.isclose(actual, expected, rel_tol=2e-6):
        raise RuntimeError(f"{label} cell volume is {actual:.6f}, expected {expected:.6f} A^3")
    return text


def readme(label: str, source: Path) -> str:
    volume = 100.624948 * int(label[1:]) / 100.0
    return f"""# Fe2O3 {label} finite-displacement phonons

This folder is one volume point in the five-volume Phonopy/QHA series
`V100, V102, V104, V106, V108`.

- EOS source: `../eos_birch_murnaghan/{source.name}`
- primitive-cell volume: {volume:.6f} A^3
- structure: final fixed-volume relaxed coordinates from the EOS output
- model: collinear AFM PBE+U, U(Fe 3d) = 4 eV, ortho-atomic projectors
- supercell: 2x2x2 (80 atoms)
- displacements: three central plus/minus pairs, 0.02 bohr
- independent QE force calculations: 6
- supercell k mesh: 4x4x4

Submit the force calculations from this directory:

```bash
sbatch run_displacements.sbatch
```

After all six outputs contain `JOB DONE`, activate Phonopy and run:

```bash
bash collect_forces.sh
phonopy-load phonopy_disp.yaml --config band.conf --save
phonopy-load phonopy_disp.yaml --config mesh.conf
```

Keep all electronic, magnetic, displacement, and sampling settings identical
between volumes. Do not relax a displaced supercell.
"""


def create_source_folder(label: str, source: Path) -> Path:
    target = ROOT / f"phonopy_{label}_2x2x2"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite existing folder: {target}")
    target.mkdir()

    positions = final_positions(source)
    (target / "unitcell.in").write_text(make_unitcell(label, positions))

    for name in ("band.conf", "mesh.conf", "collect_forces.sh"):
        shutil.copy2(TEMPLATE / name, target / name)

    prepare = (TEMPLATE / "prepare_qe_inputs.py").read_text().replace("V100", label)
    (target / "prepare_qe_inputs.py").write_text(prepare)

    batch = (TEMPLATE / "run_displacements.sbatch").read_text()
    batch = batch.replace("#SBATCH --job-name=Fe2O3FD", f"#SBATCH --job-name=Fe2O3-{label}-FD")
    (target / "run_displacements.sbatch").write_text(batch)
    shutil.copymode(TEMPLATE / "run_displacements.sbatch", target / "run_displacements.sbatch")

    (target / "README.md").write_text(readme(label, source))
    return target


def generate_displacements(target: Path, phonopy: str) -> None:
    command = [
        phonopy,
        "--qe",
        "-d",
        "--dim", "2", "2", "2",
        "--amplitude", "0.02",
        "--pm",
        "--magmom", "0", "0", "0", "0", "0", "0", "4", "4", "-4", "-4",
        "-c", "unitcell.in",
    ]
    subprocess.run(command, cwd=target, check=True)
    subprocess.run(["python3", "prepare_qe_inputs.py"], cwd=target, check=True)


def validate(target: Path) -> None:
    supercells = sorted(target.glob("supercell-[0-9][0-9][0-9].in"))
    inputs = sorted(target.glob("displacements/disp-*/alpha_Fe2O3.fd.scf.in"))
    if len(supercells) != 6 or len(inputs) != 6:
        raise RuntimeError(
            f"{target.name}: expected six supercells and six QE inputs, "
            f"found {len(supercells)} and {len(inputs)}"
        )
    for path in inputs:
        text = path.read_text()
        required = ("nat               = 80", "4 4 4 0 0 0", "Hubbard_U(2)      = 4.0")
        if any(item not in text for item in required):
            raise RuntimeError(f"Generated input failed validation: {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phonopy",
        default=shutil.which("phonopy"),
        help="Path to the Phonopy executable",
    )
    args = parser.parse_args()
    if not args.phonopy:
        raise SystemExit("Phonopy was not found; pass its path with --phonopy")

    created: list[Path] = []
    for label, source in VOLUMES.items():
        target = create_source_folder(label, source)
        generate_displacements(target, args.phonopy)
        validate(target)
        created.append(target)
        print(f"Ready: {target.relative_to(ROOT.parent)}")

    print(f"Created and validated {len(created)} QHA phonon folders.")


if __name__ == "__main__":
    main()
