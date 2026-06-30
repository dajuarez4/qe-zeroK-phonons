import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EXPERIMENT_CSV = ROOT / "sapphire_20K_experiment_points.csv"
QE_INPUT_FILE = ROOT / "Al2O3.matdyn_exp_paper_path.in"
QE_FILE = ROOT / "Al2O3.exp_paper_path.freq.gp"

THZ_TO_RADPS = 2.0 * math.pi
THZ_PER_CM1 = 33.35640952

OUTPUT_PNG = ROOT / "Al2O3_experiment_paper_style.png"
OUTPUT_PDF = ROOT / "Al2O3_experiment_paper_style.pdf"
OUTPUT_POINTS_CSV = ROOT / "Al2O3_experiment_paper_style_points.csv"


def load_qe_path(path: Path) -> tuple[list[float], list[list[float]]]:
    rows: list[list[float]] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            rows.append([float(value) for value in parts])

    if not rows:
        raise ValueError(f"No data found in {path}")

    x_values = [row[0] for row in rows]
    x0 = x_values[0]
    x_values = [value - x0 for value in x_values]

    mode_count = len(rows[0]) - 1
    freqs_by_mode: list[list[float]] = [[] for _ in range(mode_count)]
    for row in rows:
        for mode_index in range(mode_count):
            freq_cm1 = max(row[mode_index + 1], 0.0)
            freq_radps = freq_cm1 * THZ_TO_RADPS / THZ_PER_CM1
            freqs_by_mode[mode_index].append(freq_radps)

    return x_values, freqs_by_mode


def load_path_ticks(input_path: Path, freq_path: Path) -> list[float]:
    input_lines = [line.strip() for line in input_path.read_text().splitlines() if line.strip()]
    end_control = input_lines.index("/")
    point_count = int(input_lines[end_control + 1].split()[0])
    q_lines = input_lines[end_control + 2 : end_control + 2 + point_count]
    segment_counts = [int(line.split()[-1]) for line in q_lines[:-1]]

    freq_rows = []
    with freq_path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                freq_rows.append([float(value) for value in stripped.split()])

    if not freq_rows:
        raise ValueError(f"No data found in {freq_path}")

    boundary_indices = [0]
    for count in segment_counts:
        boundary_indices.append(boundary_indices[-1] + count)

    if boundary_indices[-1] >= len(freq_rows):
        raise ValueError(
            f"Path boundary index {boundary_indices[-1]} is outside {freq_path.name} with {len(freq_rows)} rows"
        )

    return [freq_rows[index][0] for index in boundary_indices]


def map_experiment_rows(ticks: list[float]) -> list[dict[str, float | str]]:
    mapped_rows: list[dict[str, float | str]] = []
    x_gamma_1, x_l, x_x, x_gamma_2, x_z = ticks

    with EXPERIMENT_CSV.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            source_segment = row["segment"]
            xi = float(row["xi"])
            frequency_thz = float(row["frequency_thz"])
            omega_radps = frequency_thz * THZ_TO_RADPS

            # The 1993 neutron tables provide Gamma-A, Gamma-D, and Gamma-Z.
            # For a plot in the later Gamma-L-X-Gamma-Z style:
            #   Gamma-L uses Gamma-A by symmetry,
            #   X-Gamma uses Gamma-D in reverse,
            #   Gamma-Z maps directly.
            # There is no direct experimental L-X scan in the table set.
            if source_segment == "gamma_a":
                mapped_segment = "gamma_l"
                x_plot = x_gamma_1 + (xi / 0.5) * (x_l - x_gamma_1)
            elif source_segment == "gamma_d":
                mapped_segment = "x_gamma"
                x_plot = x_x + ((0.5 - xi) / 0.5) * (x_gamma_2 - x_x)
            elif source_segment == "gamma_z":
                mapped_segment = "gamma_z"
                x_plot = x_gamma_2 + (xi / 0.5) * (x_z - x_gamma_2)
            else:
                continue

            mapped_rows.append(
                {
                    "source_segment": source_segment,
                    "mapped_segment": mapped_segment,
                    "xi": xi,
                    "frequency_thz": frequency_thz,
                    "omega_radps": omega_radps,
                    "x_plot": x_plot,
                }
            )

    return mapped_rows


def write_points_csv(rows: list[dict[str, float | str]], path: Path) -> None:
    fieldnames = [
        "source_segment",
        "mapped_segment",
        "xi",
        "frequency_thz",
        "omega_radps",
        "x_plot",
    ]

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def make_plot() -> None:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("Missing matplotlib. Install it to render the figure.")
        print(f"Exact remapped points were still written to {OUTPUT_POINTS_CSV.name}.")
        return

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 14,
            "axes.labelsize": 18,
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

    if not QE_FILE.exists():
        print(
            "Missing Al2O3.exp_paper_path.freq.gp. "
            "Run matdyn.x < Al2O3.matdyn_exp_paper_path.in > Al2O3.exp_paper_path.out first."
        )
        return

    ticks = load_path_ticks(QE_INPUT_FILE, QE_FILE)
    mapped_rows = map_experiment_rows(ticks)
    exp_x = [float(row["x_plot"]) for row in mapped_rows]
    exp_y = [float(row["omega_radps"]) for row in mapped_rows]

    fig, ax = plt.subplots(figsize=(6.0, 5.0))

    qe_x, qe_freqs = load_qe_path(QE_FILE)
    for mode_index, mode in enumerate(qe_freqs):
        ax.plot(
            qe_x,
            mode,
            color="#1f5aa6",
            linewidth=1.35,
            alpha=0.92,
            zorder=1,
            label="This work" if mode_index == 0 else None,
        )

    ax.scatter(
        exp_x,
        exp_y,
        color="black",
        s=12,
        linewidths=0.0,
        zorder=3,
        label="Experiment (20 K)",
    )

    for xpos in ticks[1:-1]:
        ax.axvline(xpos, color="0.65", linewidth=0.8, linestyle=":", zorder=0)

    ax.set_xlim(ticks[0], ticks[-1])
    ax.set_ylim(bottom=0.0)
    ax.set_ylabel(r"$\omega$ (rad/ps)")
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [
            r"$\Gamma$",
            "$L$\n$[\\frac{1}{2},0,0]$",
            "$X$\n$[\\frac{1}{2},0,\\frac{1}{2}]$",
            r"$\Gamma$",
            "$Z$\n$[\\frac{1}{2},\\frac{1}{2},\\frac{1}{2}]$",
        ]
    )
    ax.tick_params(top=True, right=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.legend(frameon=False, loc="upper right")

    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.18, top=0.97)
    fig.savefig(OUTPUT_PNG, dpi=500, bbox_inches="tight")
    fig.savefig(OUTPUT_PDF, bbox_inches="tight")

    print(f"saved {OUTPUT_PNG.name}")
    print(f"saved {OUTPUT_PDF.name}")


def main() -> None:
    if not QE_FILE.exists():
        print(
            "Missing Al2O3.exp_paper_path.freq.gp. "
            "Run matdyn.x < Al2O3.matdyn_exp_paper_path.in > Al2O3.exp_paper_path.out first."
        )
        return

    ticks = load_path_ticks(QE_INPUT_FILE, QE_FILE)
    mapped_rows = map_experiment_rows(ticks)
    write_points_csv(mapped_rows, OUTPUT_POINTS_CSV)
    print(f"saved {OUTPUT_POINTS_CSV.name}")
    make_plot()


if __name__ == "__main__":
    main()
