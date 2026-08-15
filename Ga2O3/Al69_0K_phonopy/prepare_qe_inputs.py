#!/usr/bin/env python3
"""Build QE force inputs from Phonopy's displaced 1x2x1 supercells."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
DISP_ROOT = ROOT / "displacements"

CONTROL = """&CONTROL
  calculation = 'scf'
  prefix      = 'Al69-Ga2O3-fd'
  outdir      = './tmp'
  pseudo_dir  = '../../pseudo'
  verbosity   = 'high'
  tstress     = .true.
  tprnfor     = .true.
  disk_io     = 'low'
/

&SYSTEM
  ibrav       = 0
  nat         = 80
  ntyp        = 3
  input_dft   = 'PBE'
  ecutwfc     = 80
  ecutrho     = 640
  occupations = 'fixed'
  nosym       = .true.
  noinv       = .true.
/

&ELECTRONS
  electron_maxstep = 300
  conv_thr         = 1.0d-12
  mixing_beta      = 0.20
  mixing_mode      = 'local-TF'
  diagonalization  = 'david'
/
"""

K_POINTS = """
K_POINTS automatic
6 6 6 0 0 0
"""


def main() -> None:
    sources = sorted(ROOT.glob("supercell-[0-9][0-9][0-9].in"))
    if not sources:
        raise SystemExit("No supercell-###.in files found; run generate_displacements.sh.")
    DISP_ROOT.mkdir(exist_ok=True)
    for index, source in enumerate(sources, start=1):
        structure = source.read_text()
        match = re.search(r"nat\s*=\s*(\d+)", structure)
        if not match or int(match.group(1)) != 80:
            raise SystemExit(f"{source.name} is not the expected 80-atom supercell")
        structure = "\n".join(
            line for line in structure.splitlines() if not line.lstrip().startswith("!")
        ).strip()
        target_dir = DISP_ROOT / f"disp-{index:03d}"
        target_dir.mkdir(exist_ok=True)
        (target_dir / "Al69_Ga2O3.fd.scf.in").write_text(f"{CONTROL}\n{structure}\n{K_POINTS}")
    print(f"Prepared {len(sources)} QE force inputs under {DISP_ROOT}")


if __name__ == "__main__":
    main()
