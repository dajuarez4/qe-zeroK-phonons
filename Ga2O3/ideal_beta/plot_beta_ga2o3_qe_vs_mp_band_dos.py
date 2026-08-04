"""Compare ideal beta-Ga2O3 QE phonons with the mp-886 PhononDB dataset."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from phonopy import load


ROOT = Path(__file__).resolve().parent
REFERENCE_DIR = ROOT / "mp-886-phonondb"
REFERENCE_YAML = REFERENCE_DIR / "phonopy_params.yaml"

QE_SEGMENTS = [
    {"file": ROOT / "beta_Ga2O3.path1.freq.gp", "labels": [r"\Gamma", "Y", "F", "L", "I"], "ticks": [0, 30, 60, 90, 120]},
    {"file": ROOT / "beta_Ga2O3.path2.freq.gp", "labels": ["I_1", "Z", "F_1"], "ticks": [0, 30, 60]},
    {"file": ROOT / "beta_Ga2O3.path3.freq.gp", "labels": ["Y", "X_1"], "ticks": [0, 30]},
    {"file": ROOT / "beta_Ga2O3.path4.freq.gp", "labels": ["X", r"\Gamma", "N"], "ticks": [0, 30, 60]},
    {"file": ROOT / "beta_Ga2O3.path5.freq.gp", "labels": ["M", r"\Gamma"], "ticks": [0, 30]},
]

# Setyawan-Curtarolo points generated for the relaxed mp-886 reference cell.
REFERENCE_KPOINTS = {
    r"\Gamma": [0.0, 0.0, 0.0],
    "Y": [0.5, 0.5, 0.0],
    "F": [0.60288841, 0.60288841, 0.41087216],
    "L": [0.5, 0.5, 0.5],
    "I": [0.74184182, 0.25815818, 0.5],
    "I_1": [0.25815818, -0.25815818, 0.5],
    "Z": [0.0, 0.0, 0.5],
    "F_1": [0.39711159, 0.39711159, 0.58912784],
    "X_1": [0.73365654, 0.26634346, 0.0],
    "X": [0.26634346, -0.26634346, 0.0],
    "N": [0.5, 0.0, 0.0],
    "M": [0.5, 0.0, 0.5],
}

CM1_PER_THZ = 33.35640952
GAP = 0.10
QE_COLOR = "black"
REFERENCE_COLOR = "#1f5aa6"
PROJECT_LABEL = "This work"
REFERENCE_LABEL = "PhononDB mp-886 (Togo)"

OUTPUT_PNG = ROOT / "beta_Ga2O3_qe_vs_mp_band_dos.png"
OUTPUT_PDF = ROOT / "beta_Ga2O3_qe_vs_mp_band_dos.pdf"
SUMMARY_FILE = ROOT / "beta_Ga2O3_qe_vs_mp_summary.txt"
GAMMA_FILE = ROOT / "beta_Ga2O3_qe_vs_mp_Gamma.csv"
BRANCH_FILE = ROOT / "beta_Ga2O3_qe_vs_mp_branch_errors.csv"
REFERENCE_BAND_JSON = ROOT / "mp-886-reference-band.json"
REFERENCE_DOS_JSON = ROOT / "mp-886-reference-dos.json"


def format_label(label: str) -> str:
    return rf"${label}$"


def format_joint_label(left: str, right: str) -> str:
    return rf"${left}\,|\,{right}$"


def interpolate_path(labels: list[str], points_per_leg: int = 31) -> np.ndarray:
    pieces = []
    for index, (left, right) in enumerate(zip(labels[:-1], labels[1:])):
        q0 = np.asarray(REFERENCE_KPOINTS[left], dtype=float)
        q1 = np.asarray(REFERENCE_KPOINTS[right], dtype=float)
        piece = np.linspace(q0, q1, points_per_leg)
        if index:
            piece = piece[1:]
        pieces.append(piece)
    return np.vstack(pieces)


def calculate_reference() -> tuple[list[dict], tuple[np.ndarray, np.ndarray], object]:
    phonon = load(REFERENCE_YAML)
    paths = [interpolate_path(segment["labels"]) for segment in QE_SEGMENTS]
    phonon.run_band_structure(paths)
    band = phonon.get_band_structure_dict()

    segments = []
    serializable = {"material_id": "mp-886", "source": str(REFERENCE_YAML), "segments": []}
    for source, qpoints, distance, frequency in zip(
        QE_SEGMENTS, band["qpoints"], band["distances"], band["frequencies"]
    ):
        segment = {
            "labels": list(source["labels"]),
            "ticks": list(source["ticks"]),
            "qpoints": np.asarray(qpoints, dtype=float),
            "distance": np.asarray(distance, dtype=float),
            "frequencies": np.asarray(frequency, dtype=float),
        }
        segments.append(segment)
        serializable["segments"].append(
            {
                "labels": segment["labels"],
                "tick_indices": segment["ticks"],
                "qpoints": segment["qpoints"].tolist(),
                "distance": segment["distance"].tolist(),
                "frequencies_THz": segment["frequencies"].tolist(),
            }
        )
    REFERENCE_BAND_JSON.write_text(json.dumps(serializable, indent=2))

    phonon.run_mesh([20, 20, 20], is_gamma_center=True)
    phonon.run_total_dos()
    dos = phonon.get_total_dos_dict()
    dos_frequency = np.asarray(dos["frequency_points"], dtype=float)
    dos_density = np.asarray(dos["total_dos"], dtype=float)
    REFERENCE_DOS_JSON.write_text(
        json.dumps(
            {
                "material_id": "mp-886",
                "mesh": [20, 20, 20],
                "frequencies_THz": dos_frequency.tolist(),
                "densities_states_per_THz_cell": dos_density.tolist(),
            },
            indent=2,
        )
    )
    return segments, (dos_frequency, dos_density), phonon


def load_qe_segment(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] != 31:
        raise ValueError(f"Expected x plus 30 branches in {path}, found {data.shape}")
    x = data[:, 0].astype(float)
    x -= x[0]
    return x, data[:, 1:].astype(float) / CM1_PER_THZ


def load_qe_dos() -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(ROOT / "beta_Ga2O3.phdos.dat", comments="#")
    frequency = data[:, 0].astype(float) / CM1_PER_THZ
    density = data[:, 1].astype(float) * CM1_PER_THZ
    return frequency, density


def map_qe_x(qe_x: np.ndarray, qe_ticks: list[int], ref_ticks: np.ndarray) -> np.ndarray:
    return np.interp(qe_x, qe_x[qe_ticks], ref_ticks)


def build_ticks(segment_ticks: list[dict]) -> tuple[list[float], list[str]]:
    positions = list(segment_ticks[0]["positions"][:-1])
    labels = [format_label(label) for label in segment_ticks[0]["labels"][:-1]]
    for left, right in zip(segment_ticks[:-1], segment_ticks[1:]):
        positions.append(0.5 * (left["positions"][-1] + right["positions"][0]))
        labels.append(format_joint_label(left["labels"][-1], right["labels"][0]))
        for position, label in zip(right["positions"][1:-1], right["labels"][1:-1]):
            positions.append(position)
            labels.append(format_label(label))
    positions.append(segment_ticks[-1]["positions"][-1])
    labels.append(format_label(segment_ticks[-1]["labels"][-1]))
    return positions, labels


def primitive_volume(cell: np.ndarray) -> float:
    return float(abs(np.linalg.det(np.asarray(cell, dtype=float))))


def write_comparison_tables(qe_arrays: list[np.ndarray], ref_segments: list[dict], phonon: object) -> None:
    qe_all = np.vstack(qe_arrays)
    ref_all = np.vstack([segment["frequencies"] for segment in ref_segments])
    if qe_all.shape != ref_all.shape:
        raise RuntimeError(f"Comparison shape mismatch: QE {qe_all.shape}, reference {ref_all.shape}")
    difference = qe_all - ref_all

    with BRANCH_FILE.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["branch", "mae_THz", "rmse_THz", "mean_signed_error_THz"])
        for branch in range(30):
            error = difference[:, branch]
            writer.writerow(
                [
                    branch + 1,
                    f"{np.mean(np.abs(error)):.8f}",
                    f"{np.sqrt(np.mean(error**2)):.8f}",
                    f"{np.mean(error):.8f}",
                ]
            )

    qe_gamma = qe_arrays[0][0]
    ref_gamma = ref_segments[0]["frequencies"][0]
    with GAMMA_FILE.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "mode",
                "qe_THz",
                "reference_THz",
                "difference_THz",
                "qe_cm-1",
                "reference_cm-1",
                "difference_cm-1",
            ]
        )
        for index, (qe_value, ref_value) in enumerate(zip(qe_gamma, ref_gamma), start=1):
            delta = qe_value - ref_value
            writer.writerow(
                [
                    index,
                    f"{qe_value:.8f}",
                    f"{ref_value:.8f}",
                    f"{delta:.8f}",
                    f"{qe_value * CM1_PER_THZ:.5f}",
                    f"{ref_value * CM1_PER_THZ:.5f}",
                    f"{delta * CM1_PER_THZ:.5f}",
                ]
            )

    qe_input = (ROOT / "beta_Ga2O3.scf_ideal.in").read_text().splitlines()
    cell_start = qe_input.index("CELL_PARAMETERS angstrom") + 1
    qe_cell = np.array([[float(value) for value in line.split()] for line in qe_input[cell_start : cell_start + 3]])
    qe_volume = primitive_volume(qe_cell)
    reference_volume = primitive_volume(phonon.primitive.cell)
    optical = ref_all > 0.5
    scale = float(np.sum(qe_all[optical] * ref_all[optical]) / np.sum(ref_all[optical] ** 2))

    text = f"""Ideal beta-Ga2O3 phonon comparison
=====================================

This work: Quantum ESPRESSO DFPT, PBE, experimental ideal C2/m cell,
2x2x2 q-grid, Fourier interpolation, 10-atom primitive cell.

Reference: Togo PhononDB/NIMS dataset keyed to Materials Project mp-886,
VASP finite displacements, relaxed C2/m structure, 1x3x2 supercell matrix,
including Born effective charges and dielectric tensor. This is an
MP-structure-derived PhononDB dataset, not the official MP/ABINIT dataset.

Structure
---------
QE primitive volume: {qe_volume:.6f} A^3
Reference primitive volume: {reference_volume:.6f} A^3
Relative QE-reference volume difference: {(qe_volume / reference_volume - 1) * 100:+.3f} %

Dispersion on Gamma-Y-F-L-I | I1-Z-F1 | Y-X1 | X-Gamma-N | M-Gamma
----------------------------------------------------------------------------
Compared values: {qe_all.size} frequencies ({qe_all.shape[0]} q points x 30 branches)
Pointwise sorted-branch MAE: {np.mean(np.abs(difference)):.6f} THz ({np.mean(np.abs(difference)) * CM1_PER_THZ:.3f} cm^-1)
Pointwise sorted-branch RMSE: {np.sqrt(np.mean(difference**2)):.6f} THz ({np.sqrt(np.mean(difference**2)) * CM1_PER_THZ:.3f} cm^-1)
Mean signed QE-reference difference: {np.mean(difference):+.6f} THz
Best through-origin QE/reference optical frequency scale: {scale:.6f}

QE path minimum: {qe_all.min():.6f} THz ({qe_all.min() * CM1_PER_THZ:.3f} cm^-1)
Reference path minimum: {ref_all.min():.6f} THz ({ref_all.min() * CM1_PER_THZ:.3f} cm^-1)
QE path maximum: {qe_all.max():.6f} THz ({qe_all.max() * CM1_PER_THZ:.3f} cm^-1)
Reference path maximum: {ref_all.max():.6f} THz ({ref_all.max() * CM1_PER_THZ:.3f} cm^-1)

The small QE negative values occur only in acoustic modes at Gamma and have
a maximum magnitude near 3.1 cm^-1. Branch-index errors are descriptive:
sorted branches can exchange identity at crossings. The QE result is an
initial 2x2x2-grid calculation; use the prepared 4x4x4 grid for a final
convergence-quality publication comparison.
"""
    SUMMARY_FILE.write_text(text)


def plot(ref_segments: list[dict], ref_dos: tuple[np.ndarray, np.ndarray], phonon: object) -> None:
    qe_dos_frequency, qe_dos = load_qe_dos()
    ref_dos_frequency, ref_dos_density = ref_dos

    all_reference = np.vstack([segment["frequencies"] for segment in ref_segments])
    qe_arrays = [load_qe_segment(segment["file"])[1] for segment in QE_SEGMENTS]
    all_qe = np.vstack(qe_arrays)
    frequency_max = 5.0 * np.ceil(
        max(all_reference.max(), all_qe.max(), qe_dos_frequency.max(), ref_dos_frequency.max()) / 5.0
    )

    qe_mask = qe_dos_frequency >= 0
    ref_mask = ref_dos_frequency >= 0
    dos_max = 1.08 * max(qe_dos[qe_mask].max(), ref_dos_density[ref_mask].max())

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 22,
            "axes.labelsize": 30,
            "xtick.labelsize": 25,
            "ytick.labelsize": 24,
            "axes.linewidth": 1.25,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.major.width": 1.1,
            "ytick.major.width": 1.1,
        }
    )

    fig, (ax_band, ax_dos) = plt.subplots(
        1,
        2,
        figsize=(16.0, 8.2),
        sharey=True,
        gridspec_kw={"width_ratios": [5.4, 1.45], "wspace": 0.05},
    )

    offset = 0.0
    tick_data = []
    for source, reference in zip(QE_SEGMENTS, ref_segments):
        qe_x, qe_frequency = load_qe_segment(source["file"])
        reference_frequency = np.maximum(reference["frequencies"], 0.0)
        qe_frequency_plot = np.maximum(qe_frequency, 0.0)
        reference_x = reference["distance"] - reference["distance"][0]
        reference_ticks = reference_x[source["ticks"]]
        qe_x_plot = map_qe_x(qe_x, source["ticks"], reference_ticks) + offset
        reference_x_plot = reference_x + offset

        for branch in reference_frequency.T:
            ax_band.plot(
                reference_x_plot,
                branch,
                color=REFERENCE_COLOR,
                linewidth=1.45,
                alpha=0.42,
                solid_capstyle="round",
                solid_joinstyle="round",
                zorder=1,
            )
        for branch in qe_frequency_plot.T:
            ax_band.plot(
                qe_x_plot,
                branch,
                color=QE_COLOR,
                linewidth=0.95,
                alpha=0.92,
                solid_capstyle="round",
                solid_joinstyle="round",
                zorder=2,
            )
        for position in reference_ticks + offset:
            ax_band.axvline(position, color="0.72", linewidth=0.9, zorder=0)

        tick_data.append(
            {
                "positions": [float(value) for value in reference_ticks + offset],
                "labels": list(source["labels"]),
            }
        )
        offset += reference_x[-1] + GAP

    tick_positions, tick_labels = build_ticks(tick_data)
    ax_dos.plot(
        ref_dos_density[ref_mask],
        ref_dos_frequency[ref_mask],
        color=REFERENCE_COLOR,
        linewidth=2.0,
        alpha=0.80,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=1,
    )
    ax_dos.plot(
        qe_dos[qe_mask],
        qe_dos_frequency[qe_mask],
        color=QE_COLOR,
        linewidth=1.55,
        solid_capstyle="round",
        solid_joinstyle="round",
        zorder=2,
    )

    legend_handles = [
        Line2D([0], [0], color=QE_COLOR, linewidth=1.7, label=PROJECT_LABEL),
        Line2D([0], [0], color=REFERENCE_COLOR, linewidth=2.2, label=REFERENCE_LABEL),
    ]

    ax_band.set_ylabel("Frequency (THz)")
    ax_band.set_xticks(tick_positions)
    ax_band.set_xticklabels(tick_labels)
    ax_band.set_xlim(tick_positions[0], tick_positions[-1])
    ax_band.set_ylim(0.0, frequency_max)
    ax_band.tick_params(axis="x", labelsize=27, pad=10, top=True)
    ax_band.tick_params(axis="y", labelsize=25, right=False)
    ax_band.spines["top"].set_visible(True)
    ax_band.spines["right"].set_visible(False)

    ax_dos.set_xlabel("PDOS", labelpad=10)
    ax_dos.set_xlim(0.0, dos_max)
    ax_dos.set_xticks([])
    ax_dos.tick_params(axis="x", which="both", top=False, bottom=False, labeltop=False, labelbottom=False)
    ax_dos.tick_params(axis="y", right=True, left=False, labelleft=False)
    ax_dos.spines["top"].set_visible(True)
    ax_dos.spines["left"].set_visible(False)
    ax_dos.spines["right"].set_visible(True)
    ax_band.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.985),
        frameon=False,
        handlelength=2.7,
        fontsize=22,
        borderaxespad=0.0,
    )

    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.15, top=0.975, wspace=0.05)
    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_PDF, bbox_inches="tight")
    plt.close(fig)

    write_comparison_tables(qe_arrays, ref_segments, phonon)


if __name__ == "__main__":
    reference_segments, reference_dos, phonon = calculate_reference()
    plot(reference_segments, reference_dos, phonon)
    print(f"saved {OUTPUT_PNG.name}")
    print(f"saved {OUTPUT_PDF.name}")
    print(f"saved {SUMMARY_FILE.name}")
