"""Build symmetry-preserving Quantum ESPRESSO inputs for ideal beta-Ga2O3.

The experimental C2/m structural parameters are from the refinement reported
in Table 3.1 of Playford (2012), based on the classic Geller structure.
The primitive cell and Setyawan-Curtarolo monoclinic k path are generated
with pymatgen so that the structure and reciprocal coordinates use the same
basis.
"""

from pathlib import Path

from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo


ROOT = Path(__file__).resolve().parent

# Refined experimental conventional cell, C2/m (unique axis b).
LATTICE = Lattice.monoclinic(
    a=12.22530,
    b=3.03636,
    c=5.80645,
    beta=103.7540,
)

# Five independent 4i sites. C2/m symmetry expands these to 20 atoms.
UNIQUE_SPECIES = ["Ga", "Ga", "O", "O", "O"]
UNIQUE_FRAC_COORDS = [
    [0.090141, 0.0, 0.79541],
    [0.158459, 0.5, 0.31417],
    [0.166230, 0.0, 0.11044],
    [0.173216, 0.0, 0.56313],
    [0.995590, 0.5, 0.25620],
]

PSEUDOS = {
    "Ga": "Ga.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
}


def make_structures() -> tuple[Structure, Structure, dict]:
    conventional = Structure.from_spacegroup(
        "C2/m", LATTICE, UNIQUE_SPECIES, UNIQUE_FRAC_COORDS
    )
    standard_primitive = SpacegroupAnalyzer(
        conventional, symprec=1.0e-3
    ).get_primitive_standard_structure(international_monoclinic=False)
    kpath = KPathSetyawanCurtarolo(standard_primitive, symprec=1.0e-3)
    primitive = kpath.prim.copy()
    primitive.sort(key=lambda site: 0 if site.species_string == "Ga" else 1)
    return conventional, primitive, kpath.kpath


def format_cell(structure: Structure) -> str:
    return "\n".join(
        "  " + " ".join(f"{value:16.10f}" for value in vector)
        for vector in structure.lattice.matrix
    )


def format_positions(structure: Structure) -> str:
    return "\n".join(
        f"{site.species_string:2s} "
        + " ".join(f"{value:16.12f}" for value in site.frac_coords)
        for site in structure
    )


def pw_input(structure: Structure, calculation: str) -> str:
    ions = """
&IONS
  ion_dynamics = 'bfgs'
/""" if calculation == "relax" else ""
    return f"""&CONTROL
  calculation = '{calculation}'
  prefix      = 'beta-Ga2O3'
  outdir      = './tmp'
  pseudo_dir  = '../pseudo'
  verbosity   = 'high'
  tstress     = .true.
  tprnfor     = .true.
/

&SYSTEM
  ibrav       = 0
  nat         = {len(structure)}
  ntyp        = 2
  input_dft   = 'PBE'
  ecutwfc     = 80
  ecutrho     = 640
  occupations = 'fixed'
/

&ELECTRONS
  electron_maxstep = 200
  conv_thr         = 1.0d-10
  mixing_beta      = 0.2
  diagonalization  = 'david'
/
{ions}

ATOMIC_SPECIES
Ga  69.723  {PSEUDOS['Ga']}
O   15.999  {PSEUDOS['O']}

CELL_PARAMETERS angstrom
{format_cell(structure)}

ATOMIC_POSITIONS crystal
{format_positions(structure)}

K_POINTS automatic
12 12 12 0 0 0
"""


def ph_input(grid: int) -> str:
    return f"""Phonons of ideal beta-Ga2O3 on a {grid}x{grid}x{grid} q grid

&INPUTPH
  prefix   = 'beta-Ga2O3'
  outdir   = './tmp'
  fildyn   = 'beta_Ga2O3.dyn'
  tr2_ph   = 1.0d-14
  ldisp    = .true.
  nq1      = {grid}
  nq2      = {grid}
  nq3      = {grid}
  epsil    = .true.
  trans    = .true.
  recover  = .true.
/
"""


def matdyn_input(labels: list[str], kpoints: dict) -> str:
    lines = [
        "&input",
        "  asr              = 'crystal',",
        "  flfrc            = 'beta_Ga2O3.fc',",
        "  flfrq            = 'PLACEHOLDER.freq',",
        "  q_in_cryst_coord = .true.,",
        "  q_in_band_form   = .true.,",
        "/",
        str(len(labels)),
    ]
    for index, label in enumerate(labels):
        qpoint = kpoints[label]
        npoints = 30 if index < len(labels) - 1 else 1
        lines.append(
            "  "
            + " ".join(f"{value:14.10f}" for value in qpoint)
            + f"  {npoints:d}"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    conventional, primitive, kpath = make_structures()
    analyzer = SpacegroupAnalyzer(primitive, symprec=1.0e-3)

    if len(conventional) != 20 or len(primitive) != 10:
        raise RuntimeError("Unexpected beta-Ga2O3 cell size")
    if analyzer.get_space_group_symbol() != "C2/m":
        raise RuntimeError(
            f"Expected C2/m, found {analyzer.get_space_group_symbol()}"
        )

    CifWriter(conventional, symprec=1.0e-3).write_file(
        ROOT / "beta_Ga2O3_ideal_conventional.cif", mode="wt"
    )
    # Keep the actual 10-atom primitive vectors/sites in this file. Asking the
    # CIF writer to emit symmetry operations conventionalizes a C-centered cell.
    CifWriter(primitive).write_file(
        ROOT / "beta_Ga2O3_ideal_primitive.cif", mode="wt"
    )
    (ROOT / "beta_Ga2O3.relax.in").write_text(
        pw_input(primitive, "relax"), encoding="utf-8"
    )
    (ROOT / "beta_Ga2O3.scf_ideal.in").write_text(
        pw_input(primitive, "scf"), encoding="utf-8"
    )
    (ROOT / "beta_Ga2O3.ph_grid_2x2x2.in").write_text(
        ph_input(2), encoding="utf-8"
    )
    (ROOT / "beta_Ga2O3.ph_grid_4x4x4.in").write_text(
        ph_input(4), encoding="utf-8"
    )
    (ROOT / "beta_Ga2O3.q2r.in").write_text(
        """&input
  fildyn = 'beta_Ga2O3.dyn',
  flfrc  = 'beta_Ga2O3.fc',
  zasr   = 'crystal',
/
""",
        encoding="utf-8",
    )
    (ROOT / "beta_Ga2O3.phdos.in").write_text(
        """&input
  asr    = 'crystal',
  dos    = .true.,
  flfrc  = 'beta_Ga2O3.fc',
  fldos  = 'beta_Ga2O3.phdos.dat',
  nk1    = 16,
  nk2    = 16,
  nk3    = 16,
  deltaE = 1.0,
/
""",
        encoding="utf-8",
    )

    for index, labels in enumerate(kpath["path"], start=1):
        stem = f"beta_Ga2O3.path{index}"
        content = matdyn_input(labels, kpath["kpoints"]).replace(
            "PLACEHOLDER.freq", f"{stem}.freq"
        )
        (ROOT / f"{stem}.matdyn.in").write_text(content, encoding="utf-8")

    lengths = conventional.lattice.abc
    angles = conventional.lattice.angles
    validation = f"""Ideal beta-Ga2O3 structure validation
=======================================
Conventional atoms: {len(conventional)}
Primitive atoms: {len(primitive)}
Formula: {primitive.composition.reduced_formula}
Space group: {analyzer.get_space_group_symbol()} ({analyzer.get_space_group_number()})
Point group: {analyzer.get_point_group_symbol()}
Conventional a,b,c (angstrom): {lengths[0]:.6f} {lengths[1]:.6f} {lengths[2]:.6f}
Conventional alpha,beta,gamma (degree): {angles[0]:.6f} {angles[1]:.6f} {angles[2]:.6f}
Conventional volume (angstrom^3): {conventional.volume:.6f}
Primitive volume (angstrom^3): {primitive.volume:.6f}
Path: Gamma-Y-F-L-I | I1-Z-F1 | Y-X1 | X-Gamma-N | M-Gamma
"""
    (ROOT / "structure_validation.txt").write_text(validation, encoding="utf-8")
    print(validation, end="")


if __name__ == "__main__":
    main()
