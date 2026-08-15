#!/usr/bin/env python3
"""Extract the final QE vc-relax geometry into a Phonopy-ready QE input."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "01_vc-relax.out"
TARGET = ROOT / "unitcell.in"
NAT = 40

HEADER = """&CONTROL
  calculation = 'scf'
  prefix      = 'Al69-Ga2O3-phonon'
  outdir      = './tmp_scf'
  pseudo_dir  = './pseudo'
  verbosity   = 'high'
  tstress     = .true.
  tprnfor     = .true.
/

&SYSTEM
  ibrav       = 0
  nat         = 40
  ntyp        = 3
  input_dft   = 'PBE'
  ecutwfc     = 80
  ecutrho     = 640
  occupations = 'fixed'
/

&ELECTRONS
  electron_maxstep = 300
  conv_thr         = 1.0d-12
  mixing_beta      = 0.20
  mixing_mode      = 'local-TF'
  diagonalization  = 'david'
/

ATOMIC_SPECIES
O   15.999  O.pbe-n-kjpaw_psl.1.0.0.UPF
Ga  69.723  Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
Al  26.982  Al.pbe-n-kjpaw_psl.1.0.0.UPF
"""


def main() -> None:
    if not OUTPUT.exists():
        raise SystemExit(f"Missing {OUTPUT.name}; run the vc-relax first.")
    text = OUTPUT.read_text(errors="replace")
    if "JOB DONE" not in text:
        raise SystemExit("The vc-relax output does not contain JOB DONE.")
    blocks = re.findall(
        r"CELL_PARAMETERS\s*\(([^)]+)\)\s*\n"
        r"((?:\s*[-+0-9.EeDd]+\s+[-+0-9.EeDd]+\s+[-+0-9.EeDd]+\s*\n){3})"
        r"\s*ATOMIC_POSITIONS\s*\(([^)]+)\)\s*\n"
        r"((?:\s*[A-Za-z][A-Za-z0-9]*\s+[-+0-9.EeDd]+\s+[-+0-9.EeDd]+\s+[-+0-9.EeDd]+[^\n]*\n){40})",
        text,
    )
    if not blocks:
        raise SystemExit("Could not locate the final CELL_PARAMETERS/ATOMIC_POSITIONS block.")
    cell_unit, cell, pos_unit, positions = blocks[-1]
    if "angstrom" not in cell_unit.lower() or "crystal" not in pos_unit.lower():
        raise SystemExit(f"Unexpected final units: cell={cell_unit}, positions={pos_unit}")

    clean_positions = []
    for line in positions.splitlines()[:NAT]:
        fields = line.split()
        clean_positions.append(f"{fields[0]:2s}  {float(fields[1]): .12f}  {float(fields[2]): .12f}  {float(fields[3]): .12f}")
    if len(clean_positions) != NAT:
        raise SystemExit(f"Expected {NAT} atoms, found {len(clean_positions)}")

    target_text = (
        HEADER
        + "\nCELL_PARAMETERS angstrom\n"
        + cell.strip()
        + "\n\nATOMIC_POSITIONS crystal\n"
        + "\n".join(clean_positions)
        + "\n\nK_POINTS automatic\n6 12 6 0 0 0\n"
    )
    TARGET.write_text(target_text)
    print(f"Wrote relaxed 40-atom QE cell to {TARGET}")


if __name__ == "__main__":
    main()
