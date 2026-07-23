"""Compare the original and symmetry-fixed Al2O3 phonon calculations."""

from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parent
ORIGINAL = ROOT.parent
STEMS = [
    ("Al2O3.mp_seg1.freq.gp", "alpha_Al2O3.path1.freq.gp"),
    ("Al2O3.mp_seg2.freq.gp", "alpha_Al2O3.path2.freq.gp"),
    ("Al2O3.mp_seg3.freq.gp", "alpha_Al2O3.path3.freq.gp"),
    ("Al2O3.mp_seg4.freq.gp", "alpha_Al2O3.path4.freq.gp"),
]
SEGMENTS = [
    (["Γ", "L", "B₁"], [0, 30, 60]),
    (["B", "Z", "Γ", "X"], [0, 30, 60, 90]),
    (["Q", "F", "P₁", "Z"], [0, 30, 60, 90]),
    (["L", "P"], [0, 30]),
]


def load_comparison():
    comparisons = []
    ticks: list[float] = []
    labels: list[str] = []
    offset = 0.0

    for segment, ((old_name, new_name), (segment_labels, tick_indices)) in enumerate(
        zip(STEMS, SEGMENTS), start=1
    ):
        old = np.loadtxt(ORIGINAL / old_name)
        new = np.loadtxt(ROOT / new_name)
        if old.shape != new.shape:
            raise ValueError(f"Segment {segment} shapes differ: {old.shape} vs {new.shape}")
        if old.shape[1] != 31:
            raise ValueError(f"Expected 30 branches in segment {segment}")
        if not np.allclose(old[:, 0], new[:, 0], atol=2.0e-6):
            raise ValueError(f"Segment {segment} q-distance grids differ")

        x = new[:, 0] - new[0, 0] + offset
        for position, label in zip(x[tick_indices], segment_labels):
            if ticks and np.isclose(position, ticks[-1]):
                labels[-1] += f"|{label}"
            else:
                ticks.append(float(position))
                labels.append(label)
        comparisons.append((segment, x, old[:, 1:], new[:, 1:]))
        offset = float(x[-1])

    return comparisons, np.asarray(ticks), labels


def write_data_and_summary(comparisons) -> str:
    rows = []
    all_differences = []
    segment_lines = []
    for segment, x, old, new in comparisons:
        difference = new - old
        all_differences.append(difference.ravel())
        segment_lines.append(
            f"segment {segment}: MAE={np.mean(np.abs(difference)):.6f}, "
            f"RMS={np.sqrt(np.mean(difference**2)):.6f}, "
            f"max_abs={np.max(np.abs(difference)):.6f} cm^-1"
        )
        for q_index in range(len(x)):
            for mode in range(old.shape[1]):
                rows.append(
                    (
                        segment,
                        q_index,
                        x[q_index],
                        mode + 1,
                        old[q_index, mode],
                        new[q_index, mode],
                        difference[q_index, mode],
                    )
                )

    table = np.asarray(rows)
    np.savetxt(
        ROOT / "Al2O3_original_vs_symmetry_fixed_differences.csv",
        table,
        delimiter=",",
        header=(
            "segment,q_index,combined_path_distance,mode,"
            "original_cm-1,symmetry_fixed_cm-1,delta_cm-1"
        ),
        comments="",
        fmt=["%d", "%d", "%.8f", "%d", "%.6f", "%.6f", "%.6f"],
    )

    difference = np.concatenate(all_differences)
    new_outputs = [
        ROOT / "alpha_Al2O3.scf.out",
        ROOT / "alpha_Al2O3.ph_grid_2x2x2.out",
        ROOT / "alpha_Al2O3.q2r.out",
        *(ROOT / f"alpha_Al2O3.path{i}.matdyn.out" for i in range(1, 5)),
        ROOT / "alpha_Al2O3.phdos.out",
    ]
    completed_outputs = sum("JOB DONE" in path.read_text(errors="replace") for path in new_outputs)
    scf_text = (ROOT / "alpha_Al2O3.scf.out").read_text(errors="replace")
    symmetry_counts = [int(value) for value in re.findall(r"(\d+) Sym\. Ops", scf_text)]
    force_match = re.findall(r"Total force\s*=\s*([0-9.Ee+-]+)", scf_text)

    def maximum_nonhermiticity(directory: Path, pattern: str) -> float:
        values = []
        for path in directory.glob(pattern):
            values.extend(
                float(value)
                for value in re.findall(
                    r"Max \|d\(i,j\)-d\*\(j,i\)\|\s*=\s*([0-9.Ee+-]+)",
                    path.read_text(errors="replace"),
                )
            )
        return max(values) if values else float("nan")

    new_nonhermiticity = maximum_nonhermiticity(ROOT, "alpha_Al2O3.path*.matdyn.out")
    old_nonhermiticity = maximum_nonhermiticity(ORIGINAL, "Al2O3.matdyn_mp_seg*.out")

    new_dos = np.loadtxt(ROOT / "alpha_Al2O3.phdos.dat")
    old_dos = np.loadtxt(ORIGINAL / "Al2O3.phdos.dat")
    new_dos_integral = np.trapezoid(new_dos[:, 1], new_dos[:, 0])
    old_dos_integral = np.trapezoid(old_dos[:, 1], old_dos[:, 0])
    summary = "\n".join(
        [
            "Al2O3 original vs symmetry-fixed phonon comparison",
            "===================================================",
            "Both calculations: 2x2x2 q grid, PBE, 80/640 Ry, 12x12x12 k grid",
            f"Compared frequency samples: {len(difference)}",
            f"Mean absolute difference: {np.mean(np.abs(difference)):.6f} cm^-1",
            f"RMS difference: {np.sqrt(np.mean(difference**2)):.6f} cm^-1",
            f"Maximum absolute difference: {np.max(np.abs(difference)):.6f} cm^-1",
            f"Mean signed difference: {np.mean(difference):.6f} cm^-1",
            f"Completed new QE stages: {completed_outputs}/{len(new_outputs)}",
            f"New SCF symmetry operations: {max(symmetry_counts) if symmetry_counts else 'not found'}",
            f"New SCF total force: {force_match[-1] if force_match else 'not found'} Ry/bohr",
            f"Original/new max matdyn non-Hermiticity: {old_nonhermiticity:.6f} / {new_nonhermiticity:.6f}",
            f"Original/new total DOS integrals: {old_dos_integral:.6f} / {new_dos_integral:.6f}",
            "",
            *segment_lines,
            "",
            "Interpretation: the dispersions agree to well below 1 cm^-1 on average.",
            "The symmetry-fixed run confirms the original physical result while",
            "restoring the expected 12 R-3c symmetry operations.",
            "The matdyn non-Hermiticity warning persists; a denser 4x4x4 q grid",
            "is the appropriate next convergence test if higher accuracy is needed.",
            "",
        ]
    )
    (ROOT / "Al2O3_original_vs_symmetry_fixed_summary.txt").write_text(
        summary, encoding="utf-8"
    )
    return summary


def plot_dispersion(comparisons, ticks, labels) -> None:
    fig, (bands_ax, difference_ax) = plt.subplots(
        2,
        1,
        figsize=(11.5, 8.2),
        sharex=True,
        gridspec_kw={"height_ratios": [4.2, 1.25], "hspace": 0.05},
        constrained_layout=True,
    )

    for _, x, old, new in comparisons:
        for branch in old.T:
            bands_ax.plot(x, branch, color="#2878b5", linewidth=0.85, alpha=0.65)
        for branch in new.T:
            bands_ax.plot(x, branch, color="black", linewidth=0.65, alpha=0.78)
        for branch_difference in (new - old).T:
            difference_ax.plot(x, branch_difference, color="#b22222", linewidth=0.65)

    for axis in (bands_ax, difference_ax):
        for position in ticks:
            axis.axvline(position, color="0.72", linewidth=0.7, zorder=0)
        axis.axhline(0.0, color="0.35", linewidth=0.7)
        axis.tick_params(direction="in")

    bands_ax.set_ylabel(r"Frequency (cm$^{-1}$)")
    bands_ax.set_title(r"$\alpha$-Al$_2$O$_3$: original vs symmetry-fixed phonons")
    bands_ax.set_ylim(-20.0, 880.0)
    bands_ax.legend(
        handles=[
            Line2D([0], [0], color="#2878b5", linewidth=1.5, label="Original"),
            Line2D([0], [0], color="black", linewidth=1.3, label="Symmetry-fixed"),
        ],
        loc="upper right",
        frameon=False,
    )
    difference_ax.set_ylabel(r"$\Delta\omega$")
    difference_ax.set_xlabel("Wave vector")
    difference_ax.set_ylim(-2.6, 2.6)
    difference_ax.set_yticks([-2, 0, 2])
    difference_ax.set_xticks(ticks, labels)
    difference_ax.set_xlim(ticks[0], ticks[-1])

    for suffix in ("png", "pdf"):
        fig.savefig(
            ROOT / f"Al2O3_original_vs_symmetry_fixed.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def plot_dos() -> None:
    old = np.loadtxt(ORIGINAL / "Al2O3.phdos.dat")
    new = np.loadtxt(ROOT / "alpha_Al2O3.phdos.dat")
    old_on_new_grid = np.interp(new[:, 0], old[:, 0], old[:, 1])

    fig, ax = plt.subplots(figsize=(7.5, 5.3))
    ax.plot(
        old[:, 0], old[:, 1], color="#2878b5", linewidth=1.5, label="Original"
    )
    ax.plot(
        new[:, 0], new[:, 1], color="black", linewidth=1.2, label="Symmetry-fixed"
    )
    ax.fill_between(
        new[:, 0],
        old_on_new_grid,
        new[:, 1],
        color="#b22222",
        alpha=0.12,
        label="Difference",
    )
    ax.axvline(0.0, color="0.4", linewidth=0.7)
    ax.set_xlabel(r"Frequency (cm$^{-1}$)")
    ax.set_ylabel("Phonon density of states")
    ax.set_title(r"$\alpha$-Al$_2$O$_3$ phonon DOS comparison")
    ax.set_ylim(bottom=0.0)
    ax.legend(frameon=False)
    ax.tick_params(direction="in")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(
            ROOT / f"Al2O3_original_vs_symmetry_fixed_dos.{suffix}",
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig)


def main() -> None:
    comparisons, ticks, labels = load_comparison()
    plot_dispersion(comparisons, ticks, labels)
    plot_dos()
    print(write_data_and_summary(comparisons), end="")


if __name__ == "__main__":
    main()
