#!/usr/bin/env python3
"""Build complete QE force inputs from Phonopy's displaced supercells."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
DISP_ROOT = ROOT / "displacements"

CONTROL = """&CONTROL
  calculation = 'scf'
  prefix      = 'alpha-Fe2O3-V104-fd'
  outdir      = './tmp'
  pseudo_dir  = '../../../pseudo'
  verbosity   = 'high'
  tstress     = .true.
  tprnfor     = .true.
  disk_io     = 'low'
/

&SYSTEM
  ibrav             = 0
  nat               = 80
  ntyp              = 3
  input_dft         = 'PBE'
  ecutwfc           = 80
  ecutrho           = 640
  occupations       = 'fixed'
  nspin             = 2
  tot_magnetization = 0.0
  starting_magnetization(2) =  0.80
  starting_magnetization(3) = -0.80
  lda_plus_u        = .true.
  Hubbard_U(2)      = 4.0
  Hubbard_U(3)      = 4.0
  U_projection_type = 'ortho-atomic'
  nosym             = .true.
  noinv             = .true.
/

&ELECTRONS
  electron_maxstep = 300
  conv_thr         = 1.0d-12
  mixing_beta      = 0.15
  mixing_mode      = 'local-TF'
  diagonalization  = 'david'
/
"""

K_POINTS = """
K_POINTS automatic
4 4 4 0 0 0
"""


def main() -> None:
    sources = sorted(ROOT.glob("supercell-[0-9][0-9][0-9].in"))
    if not sources:
        raise SystemExit("No supercell-###.in files were found.")

    DISP_ROOT.mkdir(exist_ok=True)
    for index, source in enumerate(sources, start=1):
        if source.name != f"supercell-{index:03d}.in":
            raise SystemExit(f"Unexpected displacement sequence at {source.name}")

        structure = source.read_text()
        match = re.search(r"nat\s*=\s*(\d+)", structure)
        if not match or int(match.group(1)) != 80:
            raise SystemExit(f"{source.name} is not an 80-atom supercell")

        # The leading Phonopy comment contains QE variables only as a comment.
        structure = "\n".join(
            line for line in structure.splitlines() if not line.lstrip().startswith("!")
        ).strip()

        target_dir = DISP_ROOT / f"disp-{index:03d}"
        target_dir.mkdir(exist_ok=True)
        target = target_dir / "alpha_Fe2O3.fd.scf.in"
        target.write_text(f"{CONTROL}\n{structure}\n{K_POINTS}")

    print(f"Prepared {len(sources)} QE inputs under {DISP_ROOT}")


if __name__ == "__main__":
    main()
