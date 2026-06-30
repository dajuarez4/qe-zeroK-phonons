#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "sapphire_20K_experiment_points.csv"
DEFAULT_PDF_CANDIDATES = [
    ROOT / "Lattice dynamics of sapphire (A1203).pdf",
    Path("/Users/dajuarez4/Documents/Tax-papers/Lattice dynamics of sapphire (A1203).pdf"),
]

TABLES = {
    5: {"segment": "gamma_z", "direction_label": "Gamma-Z", "representation": "A1"},
    6: {"segment": "gamma_z", "direction_label": "Gamma-Z", "representation": "A2"},
    7: {"segment": "gamma_z", "direction_label": "Gamma-Z", "representation": "A3"},
    8: {"segment": "gamma_a", "direction_label": "Gamma-A", "representation": "Sigma1"},
    9: {"segment": "gamma_a", "direction_label": "Gamma-A", "representation": "Sigma2"},
    10: {"segment": "gamma_d", "direction_label": "Gamma-D", "representation": "Sigma1_prime"},
    11: {"segment": "gamma_d", "direction_label": "Gamma-D", "representation": "Sigma2_prime"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract 20 K sapphire experimental phonon points from the PDF tables.")
    parser.add_argument("--pdf", type=Path, default=None, help="Path to the sapphire lattice-dynamics PDF.")
    parser.add_argument("--csv", type=Path, default=DEFAULT_OUTPUT, help="Output CSV for the extracted points.")
    return parser.parse_args()


def resolve_pdf_path(path: Path | None) -> Path:
    if path is not None:
        if not path.is_file():
            raise FileNotFoundError(f"PDF not found: {path}")
        return path

    for candidate in DEFAULT_PDF_CANDIDATES:
        if candidate.is_file():
            return candidate

    raise FileNotFoundError("Could not locate the sapphire PDF. Pass it explicitly with --pdf.")


def extract_text(pdf_path: Path) -> str:
    result = subprocess.run(
        ["pdftotext", "-layout", str(pdf_path), "-"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.replace("\f", "\n")


def normalize_row_prefix(line: str) -> str:
    line = re.sub(r"^\s*0\.(\d)\.(\d)", lambda m: f"0.{m.group(1)}{m.group(2)}", line)
    line = re.sub(r"(\d\.\d+)\s+\((\d+)\)", r"\1(\2)", line)
    line = re.sub(r"(?<=\d)[lI](?=\D|$)", "1", line)
    return line


def parse_q_token(token: str) -> float:
    digits = "".join(ch for ch in token if ch.isdigit())
    if len(digits) < 2:
        raise ValueError(f"Unexpected q token: {token}")
    return float(digits) / (10 ** (len(digits) - 1))


def parse_row(row_text: str) -> tuple[float, list[float]] | None:
    row_text = normalize_row_prefix(row_text)
    tokens = re.findall(r"\d+(?:\.\d+)?(?:\(\d+\))?", row_text)
    if not tokens:
        return None

    q_value = parse_q_token(tokens[0])
    frequencies: list[float] = []

    for token in tokens[1:]:
        frequencies.append(float(re.sub(r"\(\d+\)", "", token)))

    return q_value, frequencies


def parse_side_by_side_rows(lines: list[str], start: int, stop: int) -> list[dict[str, str | float | int]]:
    rows: list[dict[str, str | float | int]] = []

    for line in lines[start:stop]:
        if not re.match(r"^\s*0\.", line):
            continue

        cleaned = normalize_row_prefix(line)
        q_matches = list(re.finditer(r"(?<!\S)0\.[0-9.]+", cleaned))
        if len(q_matches) >= 2:
            split_index = q_matches[1].start()
            halves = [(5, cleaned[:split_index]), (6, cleaned[split_index:])]
        else:
            halves = [(5, cleaned)]

        for table_id, half in halves:
            parsed = parse_row(half)
            if parsed is None:
                continue

            q_value, frequencies = parsed
            table_info = TABLES[table_id]

            for branch_index, frequency_thz in enumerate(frequencies, start=1):
                rows.append(
                    {
                        "segment": table_info["segment"],
                        "direction_label": table_info["direction_label"],
                        "representation": table_info["representation"],
                        "table": table_id,
                        "branch_in_table": branch_index,
                        "xi": q_value,
                        "frequency_thz": frequency_thz,
                    }
                )

    return rows


def parse_single_table_rows(
    lines: list[str], start: int, stop: int, table_id: int
) -> list[dict[str, str | float | int]]:
    rows: list[dict[str, str | float | int]] = []
    table_info = TABLES[table_id]

    for line in lines[start:stop]:
        if not re.match(r"^\s*0\.", line):
            continue

        parsed = parse_row(line)
        if parsed is None:
            continue

        q_value, frequencies = parsed

        for branch_index, frequency_thz in enumerate(frequencies, start=1):
            rows.append(
                {
                    "segment": table_info["segment"],
                    "direction_label": table_info["direction_label"],
                    "representation": table_info["representation"],
                    "table": table_id,
                    "branch_in_table": branch_index,
                    "xi": q_value,
                    "frequency_thz": frequency_thz,
                }
            )

    return rows


def collect_rows(text: str) -> list[dict[str, str | float | int]]:
    lines = text.splitlines()
    table_line_indices: dict[int, int] = {}

    for index, line in enumerate(lines):
        for table_id in TABLES:
            if f"Table {table_id}." in line and table_id not in table_line_indices:
                table_line_indices[table_id] = index

    missing = [table_id for table_id in TABLES if table_id not in table_line_indices]
    if missing:
        raise ValueError(f"Could not locate tables in PDF text: {missing}")

    rows = []
    rows.extend(parse_side_by_side_rows(lines, table_line_indices[5] + 1, table_line_indices[7]))
    rows.extend(parse_single_table_rows(lines, table_line_indices[7] + 1, table_line_indices[8], 7))
    rows.extend(parse_single_table_rows(lines, table_line_indices[8] + 1, table_line_indices[9], 8))
    rows.extend(parse_single_table_rows(lines, table_line_indices[9] + 1, table_line_indices[10], 9))
    rows.extend(parse_single_table_rows(lines, table_line_indices[10] + 1, table_line_indices[11], 10))
    rows.extend(parse_single_table_rows(lines, table_line_indices[11] + 1, len(lines), 11))
    return rows


def write_csv(path: Path, rows: list[dict[str, str | float | int]]) -> None:
    fieldnames = [
        "segment",
        "direction_label",
        "representation",
        "table",
        "branch_in_table",
        "xi",
        "frequency_thz",
    ]

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    pdf_path = resolve_pdf_path(args.pdf)
    rows = collect_rows(extract_text(pdf_path))
    rows.sort(key=lambda row: (str(row["segment"]), float(row["xi"]), int(row["table"]), int(row["branch_in_table"])))
    write_csv(args.csv, rows)

    segment_counts: dict[str, int] = {}
    for row in rows:
        segment = str(row["segment"])
        segment_counts[segment] = segment_counts.get(segment, 0) + 1

    print(f"source PDF: {pdf_path}")
    print(f"saved {args.csv}")
    for segment, count in sorted(segment_counts.items()):
        print(f"  {segment}: {count} experimental points")


if __name__ == "__main__":
    main()
