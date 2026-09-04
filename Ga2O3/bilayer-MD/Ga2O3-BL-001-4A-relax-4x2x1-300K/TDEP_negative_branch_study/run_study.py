#!/usr/bin/env python3
"""Run an independent sensitivity study of the TDEP negative phonon branch."""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


STUDY = Path(__file__).resolve().parent
ROOT = STUDY.parent
SOURCE = ROOT / "TDEP_300K"
CASES_DIR = STUDY / "cases"
TDEP_BIN = Path(os.environ.get("TDEP_BIN_DIR", "/Users/dajuarez4/Documents/Fe/tdep/build/src"))
EXTRACT = TDEP_BIN / "extract_forceconstants/extract_forceconstants"
DISPERSION = TDEP_BIN / "phonon_dispersion_relations/phonon_dispersion_relations"


def case_definitions(nframes: int) -> list[dict]:
    all_frames = np.arange(nframes)
    half = nframes // 2
    definitions = [
        {"name": "baseline_all_5p5A", "selection": all_frames, "cutoff": 5.5, "group": "baseline"},
        {"name": "cutoff_4p5A", "selection": all_frames, "cutoff": 4.5, "group": "cutoff"},
        {"name": "cutoff_5p0A", "selection": all_frames, "cutoff": 5.0, "group": "cutoff"},
        {"name": "stride_05", "selection": all_frames[::5], "cutoff": 5.5, "group": "stride"},
        {"name": "stride_10", "selection": all_frames[::10], "cutoff": 5.5, "group": "stride"},
        {"name": "stride_20", "selection": all_frames[::20], "cutoff": 5.5, "group": "stride"},
        {"name": "window_early_half", "selection": all_frames[:half], "cutoff": 5.5, "group": "window"},
        {"name": "window_late_half", "selection": all_frames[half:], "cutoff": 5.5, "group": "window"},
        {"name": "window_discard_first_20pct", "selection": all_frames[nframes // 5 :], "cutoff": 5.5, "group": "window"},
    ]
    for count in (250, 500, 750):
        if count < nframes:
            definitions.append({
                "name": f"cumulative_{count:04d}",
                "selection": all_frames[:count],
                "cutoff": 5.5,
                "group": "cumulative",
            })
    return definitions


def parse_fit_log(text: str) -> dict:
    match = re.search(
        r"^\s*second order:\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+"
        r"<-- anharmonicity measure",
        text,
        re.MULTILINE,
    )
    if not match:
        raise RuntimeError("Could not parse TDEP fit diagnostics")
    predicted_rms, residual_rms, residual_std, r2, sigma_a = map(float, match.groups())
    return {
        "predicted_force_rms_eV_A": predicted_rms,
        "residual_force_rms_eV_A": residual_rms,
        "residual_force_std_eV_A": residual_std,
        "force_fit_R2": r2,
        "anharmonicity_sigma_A": sigma_a,
    }


def prepare_and_run() -> list[dict]:
    if not EXTRACT.is_file() or not DISPERSION.is_file():
        raise FileNotFoundError(f"TDEP executables not found below {TDEP_BIN}")
    natoms, nframes, timestep_fs, _ = (SOURCE / "infile.meta").read_text().split()
    natoms, nframes, timestep_fs = int(natoms), int(nframes), float(timestep_fs)
    positions = np.loadtxt(SOURCE / "infile.positions").reshape(nframes, natoms, 3)
    forces = np.loadtxt(SOURCE / "infile.forces").reshape(nframes, natoms, 3)
    stat_lines = (SOURCE / "infile.stat").read_text().splitlines()
    temperatures = np.array([float(line.split()[5]) for line in stat_lines])
    cases = case_definitions(nframes)
    CASES_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    for number, case in enumerate(cases, 1):
        indices = case["selection"]
        folder = CASES_DIR / case["name"]
        folder.mkdir(parents=True, exist_ok=True)
        print(f"[{number}/{len(cases)}] {case['name']}: {len(indices)} frames, cutoff {case['cutoff']} A", flush=True)
        for filename in ("infile.ucposcar", "infile.ssposcar", "infile.qpoints_dispersion"):
            shutil.copy2(SOURCE / filename, folder / filename)
        np.savetxt(folder / "infile.positions", positions[indices].reshape(-1, 3), fmt="%.16e")
        np.savetxt(folder / "infile.forces", forces[indices].reshape(-1, 3), fmt="%.16e")
        (folder / "infile.stat").write_text("\n".join(stat_lines[i] for i in indices) + "\n")
        mean_temp = float(temperatures[indices].mean())
        (folder / "infile.meta").write_text(
            f"{natoms}\n{len(indices)}\n{timestep_fs:.10f}\n{mean_temp:.10f}\n"
        )
        extract_log = folder / "extract_forceconstants.log"
        with extract_log.open("w") as log:
            subprocess.run(
                [str(EXTRACT), "--secondorder_cutoff", str(case["cutoff"]),
                 "--temperature", "300", "--firstorder", "--norotational", "--nohuang"],
                cwd=folder, stdout=log, stderr=subprocess.STDOUT, check=True,
            )
        shutil.copy2(folder / "outfile.forceconstant", folder / "infile.forceconstant")
        dispersion_log = folder / "phonon_dispersion_relations.log"
        with dispersion_log.open("w") as log:
            subprocess.run(
                [str(DISPERSION), "--readpath", "--unit", "thz"],
                cwd=folder, stdout=log, stderr=subprocess.STDOUT, check=True,
            )
        data = np.loadtxt(folder / "outfile.dispersion_relations")
        frequencies = data[:, 1:]
        minimum_index = np.unravel_index(np.argmin(frequencies), frequencies.shape)
        metrics = {
            "case": case["name"], "group": case["group"],
            "cutoff_A": case["cutoff"], "n_frames": len(indices),
            "first_frame": int(indices[0] + 1), "last_frame": int(indices[-1] + 1),
            "sample_spacing_steps": int(indices[1] - indices[0]) if len(indices) > 1 else 0,
            "sample_span_fs": float((indices[-1] - indices[0]) * timestep_fs),
            "mean_temperature_K": mean_temp,
            "temperature_std_K": float(temperatures[indices].std(ddof=1)),
            **parse_fit_log(extract_log.read_text(errors="replace")),
            "minimum_frequency_THz": float(frequencies[minimum_index]),
            "minimum_q_index": int(minimum_index[0]),
            "minimum_band_1based": int(minimum_index[1] + 1),
            "negative_values": int((frequencies < -1e-6).sum()),
        }
        (folder / "case_summary.json").write_text(json.dumps(metrics, indent=2) + "\n")
        results.append(metrics)
    return results


def temperature_autocorrelation() -> dict:
    meta = (SOURCE / "infile.meta").read_text().split()
    dt = float(meta[2])
    temperature = np.array([float(line.split()[5]) for line in (SOURCE / "infile.stat").read_text().splitlines()])
    centered = temperature - temperature.mean()
    corr = np.correlate(centered, centered, mode="full")[len(centered)-1:]
    corr /= np.arange(len(centered), 0, -1)
    corr /= corr[0]
    nonpositive = np.flatnonzero(corr <= 0)
    stop = int(nonpositive[0]) if len(nonpositive) else min(len(corr), 200)
    tau_int_steps = 0.5 + float(corr[1:stop].sum())
    return {
        "lag_fs": np.arange(len(corr)) * dt,
        "acf": corr,
        "tau_int_fs": tau_int_steps * dt,
        "effective_temperature_samples": len(temperature) / max(1.0, 2 * tau_int_steps),
    }


def analyze(results: list[dict]) -> None:
    with (STUDY / "comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader(); writer.writerows(results)
    (STUDY / "comparison.json").write_text(json.dumps(results, indent=2) + "\n")

    fig, axes = plt.subplots(4, 3, figsize=(16, 18), sharex=True, sharey=True)
    for ax, row in zip(axes.flat, results):
        data = np.loadtxt(CASES_DIR / row["case"] / "outfile.dispersion_relations")
        x, freq = data[:, 0], data[:, 1:]
        ticks = x[[0, 79, 159, 239, 319]]
        for branch in freq.T: ax.plot(x, branch, color="#245b93", lw=0.65)
        for tick in ticks: ax.axvline(tick, color="#aaaaaa", lw=0.5)
        ax.axhline(0, color="black", lw=0.7)
        ax.set_title(f"{row['case']}\n{row['n_frames']} frames; min {row['minimum_frequency_THz']:.2f} THz")
        ax.set_xticks(ticks, ["Γ", "X", "S", "Y", "Γ"])
        ax.grid(axis="y", alpha=.15)
    for ax in axes[:, 0]: ax.set_ylabel("Frequency (THz)")
    fig.suptitle("TDEP negative-branch sensitivity: all fitted dispersions", fontsize=17, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, .97))
    fig.savefig(STUDY / "01_all_dispersion_comparison.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    colors = {"baseline":"black", "cutoff":"#d88916", "stride":"#2a78b8", "window":"#b23a48", "cumulative":"#478c4b"}
    for row in results:
        data = np.loadtxt(CASES_DIR / row["case"] / "outfile.dispersion_relations")
        x, freq = data[:, 0], data[:, 1:]
        soft = freq[:, np.argmin(freq.min(axis=0))]
        ax.plot(x, soft, lw=1.6, alpha=.85, color=colors[row["group"]], label=row["case"])
    ticks = x[[0,79,159,239,319]]
    for tick in ticks: ax.axvline(tick, color="#aaaaaa", lw=.6)
    ax.axhline(0, color="black", lw=.9); ax.set_ylim(-4, 3)
    ax.set_xticks(ticks, ["Γ","X","S","Y","Γ"]); ax.set_ylabel("Soft-branch frequency (THz)")
    ax.set_title("Lowest branch under sampling and cutoff changes")
    ax.legend(fontsize=8, ncol=2, frameon=False); ax.grid(axis="y", alpha=.2)
    fig.savefig(STUDY / "02_soft_branch_overlay.png", dpi=190); plt.close(fig)

    acf = temperature_autocorrelation()
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    groups = ["cutoff", "stride", "window", "cumulative"]
    for group, marker in zip(groups, ("o","s","^","D")):
        subset = [r for r in results if r["group"] == group]
        axes[0,0].scatter([r["n_frames"] for r in subset], [r["minimum_frequency_THz"] for r in subset], label=group, marker=marker, s=55)
        axes[0,1].scatter([r["n_frames"] for r in subset], [r["anharmonicity_sigma_A"] for r in subset], label=group, marker=marker, s=55)
        axes[1,0].scatter([r["mean_temperature_K"] for r in subset], [r["minimum_frequency_THz"] for r in subset], label=group, marker=marker, s=55)
    axes[0,0].axhline(0,color="black",lw=.8); axes[0,0].set(xlabel="Frames used",ylabel="Minimum frequency (THz)")
    axes[0,1].set(xlabel="Frames used",ylabel=r"Anharmonicity $\sigma_A$")
    axes[1,0].axhline(0,color="black",lw=.8); axes[1,0].set(xlabel="Mean temperature (K)",ylabel="Minimum frequency (THz)")
    axes[1,1].plot(acf["lag_fs"][:250], acf["acf"][:250], color="#6a3d9a")
    axes[1,1].axhline(0,color="black",lw=.8); axes[1,1].set(xlabel="Temperature lag (fs)",ylabel="Autocorrelation",title=f"Temperature ACF; τint ≈ {acf['tau_int_fs']:.1f} fs")
    for ax in axes.flat: ax.grid(alpha=.2)
    axes[0,0].legend(frameon=False)
    fig.suptitle("Negative-mode diagnostics", fontsize=17, weight="bold")
    fig.savefig(STUDY / "03_diagnostic_metrics.png", dpi=190); plt.close(fig)

    baseline = next(r for r in results if r["group"] == "baseline")
    late = next(r for r in results if r["case"] == "window_late_half")
    stride20 = next(r for r in results if r["case"] == "stride_20")
    cutoffs = [r for r in results if r["group"] in ("baseline", "cutoff")]
    conclusion = (
        "The negative mode is robust across the tested sampling choices."
        if all(r["minimum_frequency_THz"] < -0.1 for r in results)
        else "The negative mode is sensitive to at least one sampling or cutoff choice."
    )
    cutoff_text = ", ".join(
        f"{row['cutoff_A']:.1f} Å → {row['minimum_frequency_THz']:.3f} THz"
        for row in sorted(cutoffs, key=lambda item: item["cutoff_A"])
    )
    report = f"""# TDEP negative-branch study

This is an independent snapshot study of the {baseline['n_frames']}-frame TDEP dataset. {conclusion}

## Key results

- Baseline minimum: {baseline['minimum_frequency_THz']:.3f} THz with σA = {baseline['anharmonicity_sigma_A']:.3f}.
- Late-half minimum: {late['minimum_frequency_THz']:.3f} THz ({late['n_frames']} frames).
- Stride-20 minimum: {stride20['minimum_frequency_THz']:.3f} THz ({stride20['n_frames']} frames).
- Cutoff minima: {cutoff_text}.
- Temperature integrated autocorrelation estimate: {acf['tau_int_fs']:.1f} fs; approximate effective temperature samples: {acf['effective_temperature_samples']:.1f}.

## Interpretation

Persistence in the late window and decorrelated fits argues against the branch being caused only by duplicated adjacent MD frames. Strong cutoff sensitivity would instead implicate the finite supercell/FC2 truncation. Persistence across both tests makes a real soft interlayer/flexural mode or a structural-reference/sum-rule issue more plausible. Absolute-valued frequencies must not be used for this diagnosis because they hide the sign.

The mode eigenvector is not assigned here: branch sorting alone is insufficient near crossings. A follow-up calculation should output eigenvectors at the minimum-q point and visualize the atomic displacement pattern.
"""
    (STUDY / "REPORT.md").write_text(report)


def main() -> None:
    results = prepare_and_run()
    analyze(results)
    print(f"Study complete: {STUDY}")


if __name__ == "__main__":
    main()
