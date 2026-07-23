#!/usr/bin/env python3
"""Transfer the final QE relaxed coordinates into the phonon SCF input."""

from pathlib import Path
import re
import shutil
import sys


NAT = 10
RELAX_OUT = Path("alpha_Fe2O3.relax.out")
SCF_IN = Path("alpha_Fe2O3.scf.in")
BACKUP = Path("alpha_Fe2O3.scf.in.before_relax")


def final_positions(text: str) -> tuple[str, list[str]]:
    markers = list(re.finditer(r"Begin final coordinates", text, re.I))
    if not markers:
        raise RuntimeError("No 'Begin final coordinates' block found")
    tail = text[markers[-1].end():]
    match = re.search(r"ATOMIC_POSITIONS\s*\(([^)]+)\)\s*\n", tail, re.I)
    if not match:
        raise RuntimeError("No ATOMIC_POSITIONS block after final-coordinate marker")
    unit = match.group(1).strip()
    lines = tail[match.end():].splitlines()[:NAT]
    if len(lines) != NAT or any(len(line.split()) < 4 for line in lines):
        raise RuntimeError(f"Expected {NAT} valid atomic-position lines")
    return unit, [line.rstrip() for line in lines]


def replace_positions(scf: str, unit: str, positions: list[str]) -> str:
    block = "ATOMIC_POSITIONS " + unit + "\n" + "\n".join(positions) + "\n\n"
    pattern = re.compile(
        r"ATOMIC_POSITIONS[^\n]*\n(?:[^\n]*\n){" + str(NAT) + r"}\s*",
        re.I,
    )
    updated, count = pattern.subn(block, scf, count=1)
    if count != 1:
        raise RuntimeError("Could not replace the SCF ATOMIC_POSITIONS block")
    return updated


def main() -> None:
    if not RELAX_OUT.exists() or not SCF_IN.exists():
        raise FileNotFoundError("Run from the Fe2O3 directory after relaxation")
    relax_text = RELAX_OUT.read_text(errors="replace")
    if "JOB DONE" not in relax_text:
        raise RuntimeError("Relaxation output does not contain JOB DONE")
    unit, positions = final_positions(relax_text)
    scf_text = SCF_IN.read_text()
    if not BACKUP.exists():
        shutil.copy2(SCF_IN, BACKUP)
    SCF_IN.write_text(replace_positions(scf_text, unit, positions))
    print(f"Updated {SCF_IN} with {NAT} final positions in {unit} units")
    print(f"Original input preserved as {BACKUP}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)

