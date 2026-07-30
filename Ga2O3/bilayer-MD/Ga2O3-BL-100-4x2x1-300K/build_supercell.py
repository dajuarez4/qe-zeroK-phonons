#!/usr/bin/env python3
"""Build the 4x2x1 (100) bilayer QE MD input from the 10-atom input."""

from pathlib import Path

SOURCE = Path("../Ga2O3-BL-100-300K/Ga2O3-BL-100.md.in")
OUTPUT = Path("Ga2O3-BL-100-4x2x1.md.in")
NX, NY, NZ = 4, 2, 1


def rows_after(lines: list[str], marker: str, count: int) -> list[list[str]]:
    start = next(i for i, line in enumerate(lines) if line.startswith(marker)) + 1
    return [lines[i].split() for i in range(start, start + count)]


lines = SOURCE.read_text().splitlines()
cell = [[float(value) for value in row] for row in rows_after(lines, "CELL_PARAMETERS", 3)]
atoms = rows_after(lines, "ATOMIC_POSITIONS", 10)

supercell = [
    [NX * value for value in cell[0]],
    [NY * value for value in cell[1]],
    [NZ * value for value in cell[2]],
]

supercell_atoms = []
for ix in range(NX):
    for iy in range(NY):
        for symbol, x, y, z in atoms:
            supercell_atoms.append(
                (symbol, (float(x) + ix) / NX, (float(y) + iy) / NY, float(z) / NZ)
            )

header = """&CONTROL
  calculation  = 'md'
  prefix       = 'Ga2O3-BL-100-4x2x1-300K'
  outdir       = './tmp'
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

output_lines = [header.rstrip()]
output_lines.extend("  " + "  ".join(f"{value:.12f}" for value in row) for row in supercell)
output_lines.extend(["", "ATOMIC_POSITIONS crystal"])
output_lines.extend(
    f"{symbol:<2}  {x:.12f}  {y:.12f}  {z:.12f}"
    for symbol, x, y, z in supercell_atoms
)
output_lines.extend(["", "K_POINTS automatic", "2 2 1 0 0 0", ""])
OUTPUT.write_text("\n".join(output_lines))

print(f"Wrote {OUTPUT} with {len(supercell_atoms)} atoms")
