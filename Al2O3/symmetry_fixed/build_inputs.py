"""Build a symmetry-preserving QE workflow from the existing Al2O3 geometry."""

from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.io.pwscf import PWInput
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo


ROOT = Path(__file__).resolve().parent
SOURCE_INPUT = ROOT.parent / "Al2O3.scf.in"
BOHR_ANGSTROM = 0.529177210903
PREFIX = "alpha-Al2O3-sym"


def source_structure() -> tuple[Structure, float, float]:
    structure = PWInput.from_file(SOURCE_INPUT).structure
    analyzer = SpacegroupAnalyzer(structure, symprec=1.0e-5)
    if analyzer.get_space_group_number() != 167:
        raise RuntimeError(
            f"Expected R-3c (167), found {analyzer.get_space_group_symbol()}"
        )

    lengths = np.asarray(structure.lattice.abc)
    angles = np.asarray(structure.lattice.angles)
    if np.ptp(lengths) > 1.0e-5 or np.ptp(angles) > 1.0e-5:
        raise RuntimeError("Source cell is not metrically rhombohedral")

    rhombohedral_a = float(lengths.mean())
    rhombohedral_angle = float(angles.mean())
    exact_lattice = Lattice.rhombohedral(rhombohedral_a, rhombohedral_angle)
    source_coords = np.asarray([site.frac_coords for site in structure])
    u_candidates = np.asarray(
        [
            1.0 - source_coords[0, 0],
            source_coords[0, 1] - 0.5,
            0.5 - source_coords[1, 0],
            source_coords[1, 2],
            source_coords[2, 0] - 0.5,
            1.0 - source_coords[2, 2],
            1.0 - source_coords[3, 1],
            source_coords[3, 2] - 0.5,
            source_coords[4, 1],
            0.5 - source_coords[4, 2],
            source_coords[5, 0],
            0.5 - source_coords[5, 1],
        ]
    )
    al_candidates = np.concatenate(
        [
            source_coords[6],
            0.5 - source_coords[7],
            source_coords[8] - 0.5,
            1.0 - source_coords[9],
        ]
    )
    oxygen_u = float(u_candidates.mean())
    aluminum_z = float(al_candidates.mean())
    exact_coords = [
        [1.0 - oxygen_u, 0.5 + oxygen_u, 0.25],
        [0.5 - oxygen_u, 0.75, oxygen_u],
        [0.5 + oxygen_u, 0.25, 1.0 - oxygen_u],
        [0.25, 1.0 - oxygen_u, 0.5 + oxygen_u],
        [0.75, oxygen_u, 0.5 - oxygen_u],
        [oxygen_u, 0.5 - oxygen_u, 0.75],
        [aluminum_z, aluminum_z, aluminum_z],
        [0.5 - aluminum_z, 0.5 - aluminum_z, 0.5 - aluminum_z],
        [0.5 + aluminum_z, 0.5 + aluminum_z, 0.5 + aluminum_z],
        [1.0 - aluminum_z, 1.0 - aluminum_z, 1.0 - aluminum_z],
    ]
    exact_structure = Structure(
        exact_lattice,
        [site.species for site in structure],
        exact_coords,
    )
    if np.max(np.abs(source_coords - np.asarray(exact_coords))) > 1.0e-6:
        raise RuntimeError("Unexpectedly large displacement during symmetrization")
    return exact_structure, rhombohedral_a, rhombohedral_angle


def format_positions(structure: Structure) -> str:
    return "\n".join(
        f"{site.species_string:2s} "
        + " ".join(f"{value:16.12f}" for value in site.frac_coords)
        for site in structure
    )


def scf_input(
    structure: Structure, rhombohedral_a: float, rhombohedral_angle: float
) -> str:
    celldm1 = rhombohedral_a / BOHR_ANGSTROM
    celldm4 = np.cos(np.radians(rhombohedral_angle))
    return f"""&CONTROL
  calculation = 'scf'
  prefix      = '{PREFIX}'
  outdir      = './tmp'
  pseudo_dir  = '../pseudo'
  verbosity   = 'high'
  tstress     = .true.
  tprnfor     = .true.
/

&SYSTEM
  ibrav       = 5
  celldm(1)   = {celldm1:.10f}
  celldm(4)   = {celldm4:.10f}
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

ATOMIC_SPECIES
O   15.999  O.pbe-n-kjpaw_psl.1.0.0.UPF
Al  26.982  Al.pbe-n-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
{format_positions(structure)}

K_POINTS automatic
12 12 12 0 0 0
"""


def ph_input(grid: int) -> str:
    return f"""Phonons of symmetry-fixed alpha-Al2O3 on a {grid}x{grid}x{grid} q grid

&INPUTPH
  prefix   = '{PREFIX}'
  outdir   = './tmp'
  fildyn   = 'alpha_Al2O3.dyn'
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


def matdyn_input(labels: list[str], kpoints: dict, stem: str) -> str:
    lines = [
        "&input",
        "  asr              = 'crystal',",
        "  flfrc            = 'alpha_Al2O3.fc',",
        f"  flfrq            = '{stem}.freq',",
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
            + f"  {npoints}"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    structure, rhombohedral_a, rhombohedral_angle = source_structure()
    analyzer = SpacegroupAnalyzer(structure, symprec=1.0e-5)
    primitive = analyzer.get_primitive_standard_structure()
    kpath = KPathSetyawanCurtarolo(primitive, symprec=1.0e-5).kpath

    expected_path = [
        ["\\Gamma", "L", "B_1"],
        ["B", "Z", "\\Gamma", "X"],
        ["Q", "F", "P_1", "Z"],
        ["L", "P"],
    ]
    if kpath["path"] != expected_path:
        raise RuntimeError(f"Unexpected rhombohedral path: {kpath['path']}")

    CifWriter(structure).write_file(
        ROOT / "alpha_Al2O3_symmetry_fixed.cif", mode="wt"
    )
    (ROOT / "alpha_Al2O3.scf.in").write_text(
        scf_input(structure, rhombohedral_a, rhombohedral_angle),
        encoding="utf-8",
    )
    for grid in (2, 4):
        (ROOT / f"alpha_Al2O3.ph_grid_{grid}x{grid}x{grid}.in").write_text(
            ph_input(grid), encoding="utf-8"
        )
    (ROOT / "alpha_Al2O3.q2r.in").write_text(
        """&input
  fildyn = 'alpha_Al2O3.dyn',
  flfrc  = 'alpha_Al2O3.fc',
  zasr   = 'crystal',
/
""",
        encoding="utf-8",
    )
    (ROOT / "alpha_Al2O3.phdos.in").write_text(
        """&input
  asr    = 'crystal',
  dos    = .true.,
  flfrc  = 'alpha_Al2O3.fc',
  fldos  = 'alpha_Al2O3.phdos.dat',
  nk1    = 16,
  nk2    = 16,
  nk3    = 16,
  deltaE = 1.0,
/
""",
        encoding="utf-8",
    )
    for index, labels in enumerate(kpath["path"], start=1):
        stem = f"alpha_Al2O3.path{index}"
        (ROOT / f"{stem}.matdyn.in").write_text(
            matdyn_input(labels, kpath["kpoints"], stem), encoding="utf-8"
        )

    validation = f"""Symmetry-fixed alpha-Al2O3 validation
=========================================
Source: ../Al2O3.scf.in
Atoms: {len(structure)}
Formula: {structure.composition.reduced_formula}
Space group: {analyzer.get_space_group_symbol()} ({analyzer.get_space_group_number()})
Point group: {analyzer.get_point_group_symbol()}
Symmetry operations: {len(analyzer.get_symmetry_operations())}
Rhombohedral a (angstrom): {rhombohedral_a:.8f}
Rhombohedral angle (degree): {rhombohedral_angle:.8f}
Primitive volume (angstrom^3): {structure.volume:.8f}
Minimum interatomic distance (angstrom): {structure.distance_matrix[structure.distance_matrix > 0].min():.8f}
Path: Gamma-L-B1 | B-Z-Gamma-X | Q-F-P1-Z | L-P

Existing-run audit:
Optical branches are positive.
Minimum interpolated frequency: -9.6361 cm^-1 (acoustic artifact).
Existing QE symmetry operations: 2 (identity and inversion only).
Expected R-3c symmetry operations: 12.
Existing matdyn maximum non-Hermiticity: 0.017658.
"""
    (ROOT / "validation.txt").write_text(validation, encoding="utf-8")
    print(validation, end="")


if __name__ == "__main__":
    main()
