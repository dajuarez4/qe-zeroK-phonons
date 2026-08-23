#!/usr/bin/env python3
"""Build AFM QE force inputs from the common displaced supercells."""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
HEADER = "&CONTROL\n  calculation = 'scf'\n  prefix      = 'alpha-Fe2O3-BMeq-AFM-fd'\n  outdir      = './tmp'\n  pseudo_dir  = '../../../../pseudo'\n  verbosity   = 'high'\n  tstress     = .true.\n  tprnfor     = .true.\n  disk_io     = 'low'\n/\n\n&SYSTEM\n  ibrav             = 0\n  nat               = 80\n  ntyp              = 3\n  input_dft         = 'PBE'\n  ecutwfc           = 80\n  ecutrho           = 640\n  occupations       = 'fixed'\n  nspin             = 2\n  tot_magnetization = 0.0\n  starting_magnetization(2) =  0.80\n  starting_magnetization(3) = -0.80\n  lda_plus_u        = .true.\n  Hubbard_U(2)      = 4.0\n  Hubbard_U(3)      = 4.0\n  U_projection_type = 'ortho-atomic'\n  nosym             = .true.\n  noinv             = .true.\n/\n\n&ELECTRONS\n  electron_maxstep = 400\n  conv_thr         = 1.0d-12\n  mixing_beta      = 0.10\n  mixing_mode      = 'local-TF'\n  mixing_ndim      = 12\n  diagonalization  = 'david'\n/\n"
K_POINTS = "\nK_POINTS automatic\n4 4 4 0 0 0\n"

sources = sorted(ROOT.glob("supercell-[0-9][0-9][0-9].in"))
if len(sources) != 6:
    raise SystemExit(f"Expected six displaced supercells, found {len(sources)}")
for index, source in enumerate(sources, 1):
    source_text = source.read_text()
    if not re.search(r"nat\s*=\s*80", source_text):
        raise SystemExit(f"{source} is not an 80-atom supercell")
    structure = "\n".join(
        line for line in source_text.splitlines()
        if not line.lstrip().startswith("!")
    ).strip()
    directory = ROOT / "displacements" / f"disp-{index:03d}"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "alpha_Fe2O3.fd.scf.in").write_text(
        HEADER + "\n" + structure + K_POINTS
    )
print(f"Prepared {len(sources)} AFM displacement inputs")
