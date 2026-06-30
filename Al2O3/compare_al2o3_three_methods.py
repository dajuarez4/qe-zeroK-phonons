#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
QE_SCF_FILE = ROOT / "Al2O3.scf.in"
QE_PH_GRID_FILE = ROOT / "Al2O3.ph_grid.in"
QE_PH_GAMMA_FILE = ROOT / "Al2O3.phG.in"
MP_BAND_FILE = ROOT / "mp-1143-band.json"
MP_DOS_FILE = ROOT / "mp-1143-dos.json"
MP_DFPT_FILE = ROOT / "mp-1143-phonon-0.json"
MP_PHEASY_FILE = ROOT / "mp-1143-phonon-1.json"
OUTPUT_CSV = ROOT / "Al2O3_three_method_summary.csv"


def determinant_3x3(matrix: list[list[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def parse_namelist_assignments(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("&") or line == "/":
                continue
            if "=" not in line:
                continue

            key, value = line.split("=", 1)
            values[key.strip().lower()] = value.strip().rstrip(",")

    return values


def parse_qe_cell(path: Path) -> list[list[float]]:
    with path.open() as handle:
        lines = handle.readlines()

    for index, raw_line in enumerate(lines):
        if raw_line.strip().startswith("CELL_PARAMETERS"):
            return [[float(item) for item in lines[index + offset].split()[:3]] for offset in range(1, 4)]

    raise ValueError(f"CELL_PARAMETERS not found in {path}")


def parse_qe_pseudos(path: Path) -> list[str]:
    pseudos: list[str] = []
    in_species_block = False

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "ATOMIC_SPECIES":
                in_species_block = True
                continue
            if in_species_block:
                if not line:
                    break
                parts = line.split()
                if len(parts) >= 3:
                    pseudos.append(parts[2])

    return pseudos


def parse_qe_kmesh(path: Path) -> str:
    with path.open() as handle:
        lines = handle.readlines()

    for index, raw_line in enumerate(lines):
        if raw_line.strip().startswith("K_POINTS"):
            return "x".join(lines[index + 1].split()[:3])

    return "n/a"


def parse_qe_qgrid(path: Path) -> str:
    values = parse_namelist_assignments(path)
    nq1 = values.get("nq1")
    nq2 = values.get("nq2")
    nq3 = values.get("nq3")
    if nq1 and nq2 and nq3:
        return f"{nq1}x{nq2}x{nq3}"
    return "n/a"


def parse_qe_epsil(path: Path) -> str:
    values = parse_namelist_assignments(path)
    return values.get("epsil", "n/a")


def load_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def format_value(value: object) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.2f}"
    return str(value)


def build_qe_summary() -> dict[str, object]:
    scf_values = parse_namelist_assignments(QE_SCF_FILE)
    cell = parse_qe_cell(QE_SCF_FILE)
    volume = abs(determinant_3x3(cell))
    pseudos = ", ".join(parse_qe_pseudos(QE_SCF_FILE))

    return {
        "label": "This work (QE)",
        "kind": "local",
        "code": "Quantum ESPRESSO",
        "phonon_method": "dfpt",
        "functional": scf_values.get("input_dft", "n/a").strip("'"),
        "basis_or_pseudos": pseudos,
        "volume_a3": volume,
        "delta_vs_qe_pct": 0.0,
        "k_mesh_or_supercell": parse_qe_kmesh(QE_SCF_FILE),
        "phonon_sampling": parse_qe_qgrid(QE_PH_GRID_FILE),
        "nac_or_epsil": parse_qe_epsil(QE_PH_GAMMA_FILE),
        "born_or_dielectric": "not computed in current inputs",
        "band_curve_cached": "yes",
        "dos_curve_cached": "yes",
        "last_updated": "n/a",
        "source_note": "Al2O3.scf.in + Al2O3.ph_grid.in",
    }


def build_mp_summary(label: str, path: Path, qe_volume: float, curve_match: str) -> dict[str, object]:
    data = load_json(path)
    structure = data.get("structure", {})
    lattice = structure.get("lattice", {})
    volume = lattice.get("volume") or data.get("volume")
    epsilon_static = data.get("epsilon_static")
    born = data.get("born")
    code = data.get("code")
    phonon_method = data.get("phonon_method")
    supercell_matrix = data.get("supercell_matrix")
    primitive_matrix = data.get("primitive_matrix")

    if supercell_matrix:
        diagonal = [str(int(row[index])) for index, row in enumerate(supercell_matrix)]
        sampling = "supercell " + "x".join(diagonal)
    else:
        sampling = "n/a"

    notes: list[str] = []
    if epsilon_static is not None and born is not None:
        notes.append("born+dielectric available")
    if primitive_matrix is not None:
        notes.append("primitive matrix stored")
    if data.get("sum_rules_breaking") is not None:
        notes.append("sum-rule diagnostics stored")

    if curve_match:
        notes.append(curve_match)

    return {
        "label": label,
        "kind": "materials-project",
        "code": code or "n/a",
        "phonon_method": phonon_method or "n/a",
        "functional": "n/a in local JSON",
        "basis_or_pseudos": "n/a in local JSON",
        "volume_a3": volume,
        "delta_vs_qe_pct": (100.0 * (volume - qe_volume) / qe_volume) if volume is not None else None,
        "k_mesh_or_supercell": "n/a",
        "phonon_sampling": sampling,
        "nac_or_epsil": "n/a in phonon doc",
        "born_or_dielectric": "yes" if epsilon_static is not None and born is not None else "unknown",
        "band_curve_cached": "yes" if curve_match else "no",
        "dos_curve_cached": "yes" if curve_match else "no",
        "last_updated": data.get("last_updated"),
        "source_note": "; ".join(notes) if notes else "n/a",
    }


def build_band_file_summary(dfpt_volume: float, pheasy_volume: float) -> str:
    band = load_json(MP_BAND_FILE)
    band_volume = band.get("structure", {}).get("lattice", {}).get("volume")
    has_nac = band.get("has_nac")
    path_convention = band.get("path_convention")

    match_note = "band volume unmatched"
    if band_volume is not None:
        if abs(band_volume - dfpt_volume) < 1.0e-6:
            match_note = "band curve matches DFPT record"
        elif abs(band_volume - pheasy_volume) < 1.0e-6:
            match_note = "band curve matches PHEASY record"

    return f"{match_note}; has_nac={has_nac}; path={path_convention}"


def print_table(rows: list[dict[str, object]]) -> None:
    columns = [
        ("Label", "label"),
        ("Code", "code"),
        ("Method", "phonon_method"),
        ("Volume(A^3)", "volume_a3"),
        ("dV_vs_QE(%)", "delta_vs_qe_pct"),
        ("K/Supercell", "k_mesh_or_supercell"),
        ("PhononSampling", "phonon_sampling"),
        ("NAC/Epsil", "nac_or_epsil"),
        ("BandCached", "band_curve_cached"),
        ("DOSCached", "dos_curve_cached"),
    ]

    table_rows = [[format_value(row[key]) for _, key in columns] for row in rows]
    widths = [len(header) for header, _ in columns]

    for row in table_rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))

    header_line = "  ".join(columns[index][0].ljust(widths[index]) for index in range(len(columns)))
    separator_line = "  ".join("-" * widths[index] for index in range(len(columns)))
    print(header_line)
    print(separator_line)

    for row in table_rows:
        print("  ".join(row[index].ljust(widths[index]) for index in range(len(row))))


def write_csv(rows: list[dict[str, object]]) -> None:
    fieldnames = list(rows[0].keys())
    with OUTPUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    qe_row = build_qe_summary()
    qe_volume = float(qe_row["volume_a3"])

    dfpt_row = build_mp_summary("MP DFPT / Petretto-style", MP_DFPT_FILE, qe_volume, "")
    pheasy_row = build_mp_summary("MP PHEASY / VASP", MP_PHEASY_FILE, qe_volume, "")

    band_note = build_band_file_summary(float(dfpt_row["volume_a3"]), float(pheasy_row["volume_a3"]))
    dfpt_row["source_note"] = f"{dfpt_row['source_note']}; {band_note}"
    if "matches DFPT record" in band_note:
        dfpt_row["band_curve_cached"] = "yes"
        dfpt_row["dos_curve_cached"] = "yes"
    if "matches PHEASY record" in band_note:
        pheasy_row["band_curve_cached"] = "yes"
        pheasy_row["dos_curve_cached"] = "yes"

    rows = [qe_row, dfpt_row, pheasy_row]

    print("Three-way phonon-method comparison for Al2O3")
    print()
    print_table(rows)
    print()
    print("Details:")
    for row in rows:
        print(
            f"  {row['label']}: functional={row['functional']}, "
            f"pseudos/basis={row['basis_or_pseudos']}, "
            f"born/dielectric={row['born_or_dielectric']}, "
            f"last_updated={row['last_updated']}, "
            f"notes={row['source_note']}"
        )
    print()
    print("Takeaways:")
    print(
        f"  The cached MP band/DOS files match the older DFPT record volume "
        f"({float(dfpt_row['volume_a3']):.2f} A^3), not the newer PHEASY/VASP record "
        f"({float(pheasy_row['volume_a3']):.2f} A^3)."
    )
    print(
        f"  Your QE cell is larger than both MP records: "
        f"{abs(float(dfpt_row['delta_vs_qe_pct'])):.2f}% larger than the DFPT record and "
        f"{abs(float(pheasy_row['delta_vs_qe_pct'])):.2f}% larger than the PHEASY record."
    )
    print(
        "  The newer PHEASY/VASP record is present as metadata only in this repo, "
        "so it cannot yet be added as a third plotted branch/DOS curve without its band/DOS JSON."
    )

    write_csv(rows)
    print()
    print(f"Saved {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
