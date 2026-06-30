import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
EXPERIMENT_CSV = ROOT / "sapphire_20K_experiment_points.csv"
CM1_PER_THZ = 33.35640952

DIRECTIONS = [
    {
        "key": "gamma_z",
        "label": r"$\Gamma \rightarrow Z$",
        "endpoint_label": r"$Z$",
        "qe_file": ROOT / "Al2O3.exp_gz.freq.gp",
        "input_file": ROOT / "Al2O3.matdyn_exp_gz.in",
        "out_file": ROOT / "Al2O3.exp_gz.out",
    },
    {
        "key": "gamma_a",
        "label": r"$\Gamma \rightarrow A$",
        "endpoint_label": r"$A$",
        "qe_file": ROOT / "Al2O3.exp_ga.freq.gp",
        "input_file": ROOT / "Al2O3.matdyn_exp_ga.in",
        "out_file": ROOT / "Al2O3.exp_ga.out",
    },
    {
        "key": "gamma_d",
        "label": r"$\Gamma \rightarrow D$",
        "endpoint_label": r"$D$",
        "qe_file": ROOT / "Al2O3.exp_gd.freq.gp",
        "input_file": ROOT / "Al2O3.matdyn_exp_gd.in",
        "out_file": ROOT / "Al2O3.exp_gd.out",
    },
]

OUTPUT_PNG = ROOT / "Al2O3_qe_vs_experiment_20K.png"
OUTPUT_PDF = ROOT / "Al2O3_qe_vs_experiment_20K.pdf"


def load_qe_segment(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected format in {path}")

    x = data[:, 0].astype(float)
    x -= x[0]
    freqs = data[:, 1:].astype(float)
    freqs[freqs < 0.0] = 0.0
    return x, freqs


def load_experiment_points() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    grouped: dict[str, list[tuple[float, float]]] = {}

    with EXPERIMENT_CSV.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            segment = row["segment"]
            xi = float(row["xi"])
            frequency_cm1 = float(row["frequency_thz"]) * CM1_PER_THZ
            grouped.setdefault(segment, []).append((xi, frequency_cm1))

    return {
        key: (
            np.asarray([pair[0] for pair in values], dtype=float),
            np.asarray([pair[1] for pair in values], dtype=float),
        )
        for key, values in grouped.items()
    }


plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 13,
        "axes.labelsize": 17,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.major.width": 1.05,
        "ytick.major.width": 1.05,
    }
)

experiment = load_experiment_points()
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.8), sharey=True)

qe_plotted = False

for ax, direction in zip(axes, DIRECTIONS):
    exp_xi, exp_freqs = experiment[direction["key"]]
    qe_file = direction["qe_file"]

    if qe_file.exists():
        qe_x, qe_freqs = load_qe_segment(qe_file)
        x_scale = qe_x[-1] / 0.5 if qe_x[-1] > 0.0 else 1.0

        for mode in range(qe_freqs.shape[1]):
            ax.plot(
                qe_x,
                qe_freqs[:, mode],
                color="#1f5aa6",
                linewidth=1.2,
                alpha=0.9,
                zorder=1,
                label="QE" if not qe_plotted and mode == 0 else None,
            )

        qe_plotted = True
    else:
        x_scale = 1.0
        print(
            f"Missing {qe_file.name}. Run matdyn.x < {direction['input_file'].name} > {direction['out_file'].name} first."
        )

    ax.scatter(
        exp_xi * x_scale,
        exp_freqs,
        color="black",
        s=13,
        linewidths=0.0,
        zorder=3,
        label="Experiment (20 K)" if direction["key"] == DIRECTIONS[0]["key"] else None,
    )
    ax.set_title(direction["label"])
    ax.set_xlabel("Wave Vector")
    ax.set_xticks([0.0, 0.5 * x_scale])
    ax.set_xticklabels([r"$\Gamma$", direction["endpoint_label"]])
    ax.tick_params(top=True, right=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)

axes[0].set_ylabel(r"Frequency (cm$^{-1}$)")

handles, labels = axes[0].get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))

fig.subplots_adjust(left=0.085, right=0.985, bottom=0.16, top=0.84, wspace=0.10)
fig.savefig(OUTPUT_PNG, dpi=500, bbox_inches="tight")
fig.savefig(OUTPUT_PDF, bbox_inches="tight")

print(f"saved {OUTPUT_PNG.name}")
print(f"saved {OUTPUT_PDF.name}")
