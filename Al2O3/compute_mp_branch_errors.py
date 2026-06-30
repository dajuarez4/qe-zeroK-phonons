#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parent

QE_SEGMENTS = [
    {"file": ROOT / "Al2O3.mp_seg1.freq.gp", "labels": [r"\Gamma", "L", "B_1"], "tick_indices": [0, 30, 60]},
    {
        "file": ROOT / "Al2O3.mp_seg2.freq.gp",
        "labels": ["B", "Z", r"\Gamma", "X"],
        "tick_indices": [0, 30, 60, 90],
    },
    {
        "file": ROOT / "Al2O3.mp_seg3.freq.gp",
        "labels": ["Q", "F", "P_1", "Z"],
        "tick_indices": [0, 30, 60, 90],
    },
    {"file": ROOT / "Al2O3.mp_seg4.freq.gp", "labels": ["L", "P"], "tick_indices": [0, 30]},
]

MP_BAND_FILE = ROOT / "mp-1143-band.json"
MP_LABEL_SEQUENCES = [
    [r"\Gamma", "L", "B_1"],
    ["B", "Z", r"\Gamma", "X"],
    ["Q", "F", "P_1", "Z"],
    ["L", "P"],
]

CM1_PER_THZ = 33.35640952
DEFAULT_MIN_REFERENCE_CM1 = 10.0
DEFAULT_CSV = ROOT / "Al2O3_mp_branch_errors.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare QE phonon branches against the Materials Project / "
            "Petretto et al. (2018) branch data on the shared segmented path."
        )
    )
    parser.add_argument(
        "--min-reference-cm1",
        type=float,
        default=DEFAULT_MIN_REFERENCE_CM1,
        help=(
            "Only Petretto points at or above this frequency are included in "
            "the MAPE calculation. Absolute-error metrics still use all points."
        ),
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help="CSV output path for the branch-by-branch summary.",
    )
    return parser.parse_args()


def load_qe_segment(path: Path) -> tuple[list[float], list[list[float]]]:
    x_values: list[float] = []
    rows: list[list[float]] = []

    with path.open() as handle:
        for raw_line in handle:
            parts = raw_line.split()
            if not parts:
                continue

            values = [float(part) for part in parts]
            x_values.append(values[0])
            rows.append([max(value, 0.0) for value in values[1:]])

    if not rows:
        raise ValueError(f"No data found in {path}")

    origin = x_values[0]
    x_values = [value - origin for value in x_values]
    return x_values, rows


def transpose_rows(rows: list[list[float]]) -> list[list[float]]:
    return [list(column) for column in zip(*rows)]


def cumulative_distance(qpoints: list[list[float]]) -> list[float]:
    if not qpoints:
        return []

    distances = [0.0]
    total = 0.0

    for previous, current in zip(qpoints, qpoints[1:]):
        total += math.dist(previous, current)
        distances.append(total)

    return distances


def find_label_indices(
    qpoints: list[list[float]],
    labels_dict: dict[str, list[float]],
    label_sequence: list[str],
    start_index: int,
) -> list[int]:
    indices: list[int] = []
    start = start_index

    for label in label_sequence:
        target = labels_dict[label]
        match = None

        for i in range(start, len(qpoints)):
            if all(abs(a - b) < 1.0e-8 for a, b in zip(qpoints[i], target)):
                match = i
                break

        if match is None:
            raise ValueError(f"Could not locate label {label} in MP band data")

        indices.append(match)
        start = match + 1

    return indices


def load_mp_segments() -> list[dict[str, list[list[float]] | list[float] | list[str]]]:
    with MP_BAND_FILE.open() as handle:
        data = json.load(handle)

    qpoints = data["qpoints"]
    labels_dict = data["labels_dict"]
    raw_frequencies = data["frequencies"]
    branch_count = len(raw_frequencies)
    segments = []
    search_start = 0

    for label_sequence in MP_LABEL_SEQUENCES:
        indices = find_label_indices(qpoints, labels_dict, label_sequence, search_start)
        start = indices[0]
        stop = indices[-1] + 1
        segment_qpoints = qpoints[start:stop]
        segment_x = cumulative_distance(segment_qpoints)
        tick_positions = [segment_x[index - start] for index in indices]
        segment_rows: list[list[float]] = []

        for qpoint_index in range(start, stop):
            row = []

            for branch_index in range(branch_count):
                frequency = float(raw_frequencies[branch_index][qpoint_index]) * CM1_PER_THZ
                row.append(max(frequency, 0.0))

            segment_rows.append(row)

        segments.append(
            {
                "x": segment_x,
                "rows": segment_rows,
                "labels": label_sequence,
                "ticks": tick_positions,
            }
        )
        search_start = indices[-1] + 1

    return segments


def interpolate_series(x_source: list[float], y_source: list[float], x_target: list[float]) -> list[float]:
    if len(x_source) != len(y_source):
        raise ValueError("x_source and y_source must have the same length")
    if not x_source:
        return []
    if len(x_source) == 1:
        return [y_source[0] for _ in x_target]

    result: list[float] = []
    left = 0

    for x_value in x_target:
        if x_value <= x_source[0]:
            result.append(y_source[0])
            continue
        if x_value >= x_source[-1]:
            result.append(y_source[-1])
            continue

        while left + 1 < len(x_source) and x_source[left + 1] < x_value:
            left += 1

        x0 = x_source[left]
        x1 = x_source[left + 1]
        y0 = y_source[left]
        y1 = y_source[left + 1]

        if x1 == x0:
            result.append(y0)
            continue

        fraction = (x_value - x0) / (x1 - x0)
        result.append(y0 + fraction * (y1 - y0))

    return result


def map_qe_x_to_mp(qe_x: list[float], qe_tick_indices: list[int], mp_tick_positions: list[float]) -> list[float]:
    qe_tick_positions = [qe_x[index] for index in qe_tick_indices]
    return interpolate_series(qe_tick_positions, mp_tick_positions, qe_x)


def mean_or_nan(values: list[float]) -> float:
    return fmean(values) if values else float("nan")


def symmetric_percent_error(qe_value: float, reference_value: float) -> float | None:
    denominator = abs(qe_value) + abs(reference_value)
    if denominator == 0.0:
        return None
    return 200.0 * abs(qe_value - reference_value) / denominator


def compute_metrics(
    qe_values: list[float], reference_values: list[float], min_reference_cm1: float
) -> dict[str, float | int]:
    if len(qe_values) != len(reference_values):
        raise ValueError("QE and reference series must have the same length")

    absolute_errors: list[float] = []
    percent_errors: list[float] = []
    symmetric_percent_errors: list[float] = []

    for qe_value, reference_value in zip(qe_values, reference_values):
        absolute_error = abs(qe_value - reference_value)
        absolute_errors.append(absolute_error)

        if abs(reference_value) >= min_reference_cm1:
            percent_errors.append(100.0 * absolute_error / abs(reference_value))

        symmetric_error = symmetric_percent_error(qe_value, reference_value)
        if symmetric_error is not None:
            symmetric_percent_errors.append(symmetric_error)

    rmse = math.sqrt(sum(error * error for error in absolute_errors) / len(absolute_errors))

    return {
        "points": len(reference_values),
        "percent_points": len(percent_errors),
        "mae_cm1": mean_or_nan(absolute_errors),
        "rmse_cm1": rmse,
        "max_abs_cm1": max(absolute_errors),
        "mape_percent": mean_or_nan(percent_errors),
        "max_ape_percent": max(percent_errors) if percent_errors else float("nan"),
        "smape_percent": mean_or_nan(symmetric_percent_errors),
    }


def build_branch_series(
    mp_segments: list[dict[str, list[list[float]] | list[float] | list[str]]],
) -> tuple[list[list[float]], list[list[float]], list[dict[str, float | str | int]]]:
    if len(mp_segments) != len(QE_SEGMENTS):
        raise ValueError("Mismatch between QE and MP segment definitions")

    qe_by_branch: list[list[float]] | None = None
    mp_by_branch: list[list[float]] | None = None
    segment_summaries: list[dict[str, float | str | int]] = []

    for qe_segment, mp_segment in zip(QE_SEGMENTS, mp_segments):
        qe_x, qe_rows = load_qe_segment(qe_segment["file"])
        qe_branches = transpose_rows(qe_rows)
        mp_x = list(mp_segment["x"])
        mp_rows = list(mp_segment["rows"])
        mp_branches = transpose_rows(mp_rows)
        mp_ticks = list(mp_segment["ticks"])
        segment_labels = list(mp_segment["labels"])

        if len(qe_branches) != len(mp_branches):
            raise ValueError(f"Branch-count mismatch in segment {qe_segment['file'].name}")

        qe_x_mapped = map_qe_x_to_mp(qe_x, qe_segment["tick_indices"], mp_ticks)

        if qe_by_branch is None:
            qe_by_branch = [[] for _ in range(len(qe_branches))]
            mp_by_branch = [[] for _ in range(len(mp_branches))]

        segment_qe_all: list[float] = []
        segment_mp_all: list[float] = []

        for branch_index, (qe_branch, mp_branch) in enumerate(zip(qe_branches, mp_branches)):
            qe_interpolated = interpolate_series(qe_x_mapped, qe_branch, mp_x)
            qe_by_branch[branch_index].extend(qe_interpolated)
            mp_by_branch[branch_index].extend(mp_branch)
            segment_qe_all.extend(qe_interpolated)
            segment_mp_all.extend(mp_branch)

        segment_summaries.append(
            {
                "segment": f"{segment_labels[0]}->{segment_labels[-1]}",
                "branch_count": len(qe_branches),
                "qpoint_count": len(mp_x),
                "qe_values": segment_qe_all,
                "mp_values": segment_mp_all,
            }
        )

    if qe_by_branch is None or mp_by_branch is None:
        raise ValueError("No branch data were loaded")

    return qe_by_branch, mp_by_branch, segment_summaries


def format_metric(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    if math.isnan(value):
        return "n/a"
    return f"{value:.2f}"


def print_branch_table(branch_rows: list[dict[str, float | int]]) -> None:
    headers = ["Branch", "Points", "PctPts", "MAE(cm^-1)", "RMSE(cm^-1)", "MAPE(%)", "sMAPE(%)"]
    rows = [
        [
            str(int(row["branch"])),
            format_metric(row["points"]),
            format_metric(row["percent_points"]),
            format_metric(row["mae_cm1"]),
            format_metric(row["rmse_cm1"]),
            format_metric(row["mape_percent"]),
            format_metric(row["smape_percent"]),
        ]
        for row in branch_rows
    ]
    widths = [len(header) for header in headers]

    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(cell))

    header_line = "  ".join(header.ljust(widths[index]) for index, header in enumerate(headers))
    separator_line = "  ".join("-" * widths[index] for index in range(len(headers)))

    print(header_line)
    print(separator_line)

    for row in rows:
        print("  ".join(cell.ljust(widths[index]) for index, cell in enumerate(row)))


def write_csv(path: Path, branch_rows: list[dict[str, float | int]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "branch",
                "points",
                "percent_points",
                "mae_cm1",
                "rmse_cm1",
                "max_abs_cm1",
                "mape_percent",
                "max_ape_percent",
                "smape_percent",
            ]
        )

        for row in branch_rows:
            writer.writerow(
                [
                    int(row["branch"]),
                    int(row["points"]),
                    int(row["percent_points"]),
                    row["mae_cm1"],
                    row["rmse_cm1"],
                    row["max_abs_cm1"],
                    row["mape_percent"],
                    row["max_ape_percent"],
                    row["smape_percent"],
                ]
            )


def main() -> None:
    args = parse_args()
    mp_segments = load_mp_segments()
    qe_by_branch, mp_by_branch, segment_summaries = build_branch_series(mp_segments)
    branch_rows: list[dict[str, float | int]] = []

    print(
        "Percent errors are reported as MAPE relative to the Petretto branch values, "
        f"using only reference points >= {args.min_reference_cm1:.1f} cm^-1."
    )
    print("sMAPE uses all nonzero (QE, Petretto) pairs and is more stable near acoustic modes.")
    print()

    for branch_index, (qe_values, mp_values) in enumerate(zip(qe_by_branch, mp_by_branch), start=1):
        metrics = compute_metrics(qe_values, mp_values, args.min_reference_cm1)
        metrics["branch"] = branch_index
        branch_rows.append(metrics)

    print_branch_table(branch_rows)
    print()

    overall_metrics = compute_metrics(
        [value for series in qe_by_branch for value in series],
        [value for series in mp_by_branch for value in series],
        args.min_reference_cm1,
    )
    print(
        "Overall: "
        f"MAE = {format_metric(overall_metrics['mae_cm1'])} cm^-1, "
        f"RMSE = {format_metric(overall_metrics['rmse_cm1'])} cm^-1, "
        f"MAPE = {format_metric(overall_metrics['mape_percent'])}%, "
        f"sMAPE = {format_metric(overall_metrics['smape_percent'])}%"
    )
    print()

    print("Segment summary:")
    for summary in segment_summaries:
        metrics = compute_metrics(summary["qe_values"], summary["mp_values"], args.min_reference_cm1)
        print(
            f"  {summary['segment']}: "
            f"q-points = {summary['qpoint_count']}, "
            f"MAE = {format_metric(metrics['mae_cm1'])} cm^-1, "
            f"MAPE = {format_metric(metrics['mape_percent'])}%"
        )

    write_csv(args.csv, branch_rows)
    print()
    print(f"Saved branch summary to {args.csv}")


if __name__ == "__main__":
    main()
