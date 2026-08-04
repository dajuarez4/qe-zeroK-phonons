#!/usr/bin/env python3
"""Build the 80-atom MD input from the final relaxed 10-atom geometry."""

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RELAX_OUT = ROOT / "Ga2O3-BL-001.relax.out"
RELAX_IN = ROOT / "Ga2O3-BL-001.relax.in"
MD_IN = ROOT / "Ga2O3-BL-001-4x2x1.md.in"
NX, NY, NZ = 4, 2, 1
NAT = 10

CELL = [
    [3.246626537, 0.140350018, 0.0],
    [1.253133791, 5.730779279, 0.0],
    [0.0, 0.0, 30.0],
]


def final_atoms():
    if not RELAX_OUT.exists():
        raise SystemExit(f"Missing {RELAX_OUT.name}")
    text = RELAX_OUT.read_text(errors="replace")
    if "JOB DONE." not in text or "End of BFGS Geometry Optimization" not in text:
        raise SystemExit("Relaxation did not finish; MD input was not generated")
    if "maximum number of steps has been reached" in text:
        raise SystemExit("Relaxation reached nstep without convergence; MD input was not generated")
    lines = text.splitlines()
    blocks = [i for i, line in enumerate(lines) if line.strip().startswith("ATOMIC_POSITIONS")]
    if not blocks:
        raise SystemExit("No final ATOMIC_POSITIONS block found")
    start = blocks[-1]
    if "crystal" not in lines[start].lower():
        raise SystemExit(f"Expected crystal coordinates, found: {lines[start]}")
    atoms = []
    for line in lines[start + 1 : start + 1 + NAT]:
        fields = line.split()
        if len(fields) < 4:
            raise SystemExit("The final relaxed coordinate block is incomplete")
        atoms.append((fields[0], *(float(v) for v in fields[1:4])))
    return atoms


def initial_atoms():
    """Read the starting unit-cell atoms to create a visible MD input template."""
    lines = RELAX_IN.read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line.strip().startswith("ATOMIC_POSITIONS")]
    if not starts:
        raise SystemExit(f"No ATOMIC_POSITIONS block in {RELAX_IN.name}")
    start = starts[-1]
    atoms = []
    for line in lines[start + 1 : start + 1 + NAT]:
        fields = line.split()
        atoms.append((fields[0], *(float(v) for v in fields[1:4])))
    return atoms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--initialize",
        action="store_true",
        help="create the MD input from initial coordinates for inspection",
    )
    args = parser.parse_args()
    atoms = initial_atoms() if args.initialize else final_atoms()
    supercell = [
        [NX * v for v in CELL[0]],
        [NY * v for v in CELL[1]],
        [NZ * v for v in CELL[2]],
    ]
    super_atoms = []
    for ix in range(NX):
        for iy in range(NY):
            for symbol, x, y, z in atoms:
                super_atoms.append(
                    (symbol, (x + ix) / NX, (y + iy) / NY, z / NZ)
                )

    header = """&CONTROL
  calculation  = 'md'
  prefix       = 'Ga2O3-BL-001-4x2x1-300K'
  outdir       = './tmp/md'
  pseudo_dir   = '../../pseudo'
  verbosity    = 'high'
  restart_mode = 'from_scratch'
  nstep        = 10000
  dt           = 20.0
  tstress      = .true.
  tprnfor      = .true.
  disk_io      = 'low'
/

&SYSTEM
  ibrav           = 0
  nat             = 80
  ntyp            = 2
  input_dft       = 'PBE'
  ecutwfc         = 80
  ecutrho         = 640
  occupations     = 'fixed'
  assume_isolated = '2D'
  nosym           = .true.
/

&ELECTRONS
  electron_maxstep = 200
  conv_thr         = 1.0d-8
  mixing_beta      = 0.2
  diagonalization  = 'david'
/

&IONS
  ion_dynamics      = 'verlet'
  ion_temperature   = 'svr'
  tempw             = 300.0
  nraise            = 200
  pot_extrapolation = 'second_order'
  wfc_extrapolation = 'second_order'
/

ATOMIC_SPECIES
Ga  69.723  Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
O   15.999  O.pbe-n-kjpaw_psl.1.0.0.UPF

CELL_PARAMETERS angstrom
"""
    output = [header.rstrip()]
    output += ["  " + "  ".join(f"{v:.12f}" for v in row) for row in supercell]
    output += ["", "ATOMIC_POSITIONS crystal"]
    output += [
        f"{s:<2}  {x:.12f}  {y:.12f}  {z:.12f}"
        for s, x, y, z in super_atoms
    ]
    output += ["", "K_POINTS automatic", "2 2 1 0 0 0", ""]
    MD_IN.write_text("\n".join(output))
    print(f"Wrote {MD_IN.name} with {len(super_atoms)} atoms")
    if args.initialize:
        print("The batch job will replace these coordinates with the relaxed geometry")


if __name__ == "__main__":
    main()
