#!/usr/bin/env python3
"""Render the five final relaxed hematite structures with OVITO and build a dashboard."""

from pathlib import Path
import html
import math
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from ovito.io import import_file
from ovito.modifiers import ReplicateModifier, CreateBondsModifier
from ovito.vis import Viewport, TachyonRenderer


BOHR_TO_ANG = 0.529177210903
RY_TO_EV = 13.605693122994
VREF = 100.624948
OUTPUT_DIR = Path("relaxed_structures_ovito")


def last_float(pattern, text):
    values = re.findall(pattern, text, flags=re.I)
    if not values:
        raise ValueError(f"Missing pattern: {pattern}")
    return float(values[-1])


def qe_rhombohedral_cell(alat_bohr, cosine):
    """QE ibrav=5 primitive vectors, returned as rows in Angstrom."""
    tx = math.sqrt((1.0 - cosine) / 2.0)
    ty = math.sqrt((1.0 - cosine) / 6.0)
    tz = math.sqrt((1.0 + 2.0 * cosine) / 3.0)
    cell = np.array(
        [
            [tx, -ty, tz],
            [0.0, 2.0 * ty, tz],
            [-tx, -ty, tz],
        ]
    )
    return cell * alat_bohr * BOHR_TO_ANG


def final_positions(text):
    marker = text.lower().rfind("begin final coordinates")
    if marker < 0:
        raise ValueError("No final-coordinate block")
    tail = text[marker:]
    match = re.search(r"ATOMIC_POSITIONS\s*\(crystal\)\s*\n", tail, flags=re.I)
    if not match:
        raise ValueError("Final positions are not in crystal coordinates")
    atoms = []
    for line in tail[match.end():].splitlines()[:10]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError("Incomplete final-coordinate block")
        atoms.append((fields[0], np.array([float(value) for value in fields[1:4]])))
    return atoms


def parse_output(path):
    text = path.read_text(errors="replace")
    if "JOB DONE" not in text:
        raise ValueError(f"{path} is incomplete")
    alat = last_float(r"lattice parameter \(alat\)\s+=\s+([-+0-9.Ee]+)", text)
    cosine = last_float(r"celldm\(4\)=\s*([-+0-9.Ee]+)", text)
    volume = last_float(r"unit-cell volume\s+=\s+([-+0-9.Ee]+)", text) * BOHR_TO_ANG**3
    energy = last_float(r"!\s+total energy\s+=\s+([-+0-9.Ee]+)\s+Ry", text)
    force = last_float(r"Total force\s+=\s+([-+0-9.Ee]+)", text)
    magnetization = last_float(r"absolute magnetization\s+=\s+([-+0-9.Ee]+)", text)
    cell = qe_rhombohedral_cell(alat, cosine)
    atoms = final_positions(text)
    ratio = volume / VREF
    return {
        "label": f"V{int(round(100.0 * ratio)):03d}",
        "path": path,
        "volume": volume,
        "energy_ry": energy,
        "force": force,
        "magnetization": magnetization,
        "cell": cell,
        "atoms": atoms,
    }


def write_extxyz(record):
    path = OUTPUT_DIR / f"alpha_Fe2O3_{record['label']}_relaxed.extxyz"
    lattice = " ".join(f"{value:.10f}" for value in record["cell"].reshape(-1))
    lines = [
        str(len(record["atoms"])),
        f'Lattice="{lattice}" Properties=species:S:1:pos:R:3 pbc="T T T"',
    ]
    for species, fractional in record["atoms"]:
        cartesian = fractional @ record["cell"]
        lines.append(
            f"{species:4s} {cartesian[0]: .10f} {cartesian[1]: .10f} {cartesian[2]: .10f}"
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def render_record(record, extxyz):
    pipeline = import_file(str(extxyz))

    def style_particles(frame, data):
        types = data.particles_.particle_types_
        styling = {
            "O": ((0.86, 0.90, 0.96), 0.60),
            "Fe1": ((0.96, 0.35, 0.20), 0.82),
            "Fe2": ((0.10, 0.72, 0.78), 0.82),
        }
        for particle_type in list(types.types):
            if particle_type.name in styling:
                mutable_type = types.type_by_id_(particle_type.id)
                mutable_type.color, mutable_type.radius = styling[particle_type.name]
        data.cell_.vis.enabled = True
        data.cell_.vis.line_width = 0.18
        data.cell_.vis.rendering_color = (0.38, 0.45, 0.56)

    pipeline.modifiers.append(style_particles)
    pipeline.modifiers.append(ReplicateModifier(num_x=2, num_y=2, num_z=2, adjust_box=True))
    pipeline.modifiers.append(CreateBondsModifier(cutoff=2.30))

    def style_bonds(frame, data):
        if data.particles.bonds is not None:
            data.particles_.bonds_.vis.width = 0.11
            data.particles_.bonds_.vis.color = (0.68, 0.72, 0.78)

    pipeline.modifiers.append(style_bonds)
    pipeline.add_to_scene()
    data = pipeline.compute()
    center = np.asarray(data.cell)[:, :3].sum(axis=1) * 0.5 + np.asarray(data.cell)[:, 3]
    view_vector = np.array([1.35, -1.60, 1.05])
    view_vector /= np.linalg.norm(view_vector)
    viewport = Viewport(type=Viewport.Type.Ortho)
    viewport.camera_pos = center + 28.0 * view_vector
    viewport.camera_dir = -view_vector
    viewport.fov = 17.2
    image_path = OUTPUT_DIR / f"alpha_Fe2O3_{record['label']}_relaxed_ovito.png"
    renderer = TachyonRenderer(
        ambient_occlusion=True,
        ambient_occlusion_samples=16,
        shadows=True,
        antialiasing_samples=12,
    )
    viewport.render_image(
        filename=str(image_path), size=(900, 680),
        background=(0.025, 0.035, 0.055), renderer=renderer,
    )
    pipeline.remove_from_scene()
    return image_path


def build_dashboard(records):
    energy_min = min(record["energy_ry"] for record in records)
    colors = {"background": "#07111f", "card": "#0e1b2d", "text": "#e5edf7", "muted": "#91a3b8"}
    fig = plt.figure(figsize=(16, 10.5), facecolor=colors["background"])
    grid = fig.add_gridspec(2, 3, left=0.035, right=0.975, bottom=0.06, top=0.89, wspace=0.055, hspace=0.17)

    for index, record in enumerate(records):
        ax = fig.add_subplot(grid[index // 3, index % 3])
        ax.set_facecolor(colors["card"])
        image = plt.imread(record["image"])
        ax.imshow(image)
        delta_mev = (record["energy_ry"] - energy_min) * RY_TO_EV * 1000.0 / 2.0
        ax.set_title(
            f"{record['label']}   •   {record['volume']:.3f} Å³",
            color=colors["text"], fontsize=15, fontweight="bold", loc="left", pad=9,
        )
        ax.text(
            0.02, 0.025,
            f"ΔE {delta_mev:7.2f} meV/Fe₂O₃    Force {record['force']:.1e} Ry/bohr    |M| {record['magnetization']:.2f} μB",
            transform=ax.transAxes, color=colors["muted"], fontsize=9.2,
            bbox={"facecolor": colors["background"], "alpha": 0.86, "edgecolor": "none", "pad": 5},
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#1c3149")

    ax = fig.add_subplot(grid[1, 2])
    ax.set_facecolor(colors["card"])
    volume = np.array([record["volume"] for record in records])
    energy = np.array([record["energy_ry"] for record in records])
    relative = (energy - energy.min()) * RY_TO_EV / 2.0
    ax.plot(volume, relative, color="#f97316", lw=2.0, marker="o", ms=7)
    for record, xvalue, yvalue in zip(records, volume, relative):
        ax.annotate(record["label"], (xvalue, yvalue), xytext=(0, 8), textcoords="offset points", ha="center", color=colors["text"], fontsize=9)
    ax.set_title("Energy–volume trend", color=colors["text"], fontsize=15, fontweight="bold", loc="left", pad=12)
    ax.set_xlabel("Primitive-cell volume (Å³)", color=colors["muted"])
    ax.set_ylabel("Relative energy (eV / Fe₂O₃)", color=colors["muted"])
    ax.tick_params(colors=colors["muted"])
    ax.grid(color="#24364c", alpha=0.65, lw=0.7)
    for spine in ax.spines.values():
        spine.set_color("#1c3149")
    ax.legend(
        handles=[
            Patch(facecolor="#f55a33", label="Fe↑ (Fe1)"),
            Patch(facecolor="#1ab8c7", label="Fe↓ (Fe2)"),
            Patch(facecolor="#dbe6f5", label="O"),
        ],
        loc="upper right", frameon=False, labelcolor=colors["text"], fontsize=9,
    )

    fig.text(0.038, 0.952, "AFM α-Fe₂O₃  •  Relaxed structure dashboard", color=colors["text"], fontsize=25, fontweight="bold")
    fig.text(
        0.039, 0.914,
        "Five fixed-volume PBE+U relaxations  |  OVITO 3.10.5  |  2×2×2 periodic replicas  |  Fe–O bonds ≤ 2.30 Å",
        color=colors["muted"], fontsize=11.5,
    )
    png = OUTPUT_DIR / "alpha_Fe2O3_relaxed_structures_dashboard.png"
    pdf = OUTPUT_DIR / "alpha_Fe2O3_relaxed_structures_dashboard.pdf"
    fig.savefig(png, dpi=190, facecolor=fig.get_facecolor())
    fig.savefig(pdf, facecolor=fig.get_facecolor())
    plt.close(fig)
    return png, pdf


def build_html(records):
    energy_min = min(record["energy_ry"] for record in records)
    cards = []
    for record in records:
        delta = (record["energy_ry"] - energy_min) * RY_TO_EV * 1000.0 / 2.0
        cards.append(
            f'''<article class="card"><img src="{html.escape(record['image'].name)}" alt="{record['label']} relaxed hematite structure">
            <div class="copy"><h2>{record['label']} <span>{record['volume']:.3f} Å³</span></h2>
            <p>ΔE {delta:.2f} meV/Fe₂O₃ · Force {record['force']:.1e} Ry/bohr · |M| {record['magnetization']:.2f} μB</p></div></article>'''
        )
    page = f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>AFM α-Fe₂O₃ relaxed structures</title><style>
    :root{{--bg:#07111f;--card:#0e1b2d;--text:#e5edf7;--muted:#91a3b8;--line:#1c3149}}
    *{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--text);font-family:Inter,ui-sans-serif,system-ui,sans-serif}}
    main{{max-width:1500px;margin:auto;padding:48px 30px 64px}} h1{{margin:0;font-size:clamp(30px,4vw,54px)}} .sub{{color:var(--muted);margin:12px 0 32px}}
    .legend{{display:flex;gap:22px;margin-bottom:24px;color:var(--muted)}} .dot{{width:11px;height:11px;border-radius:50%;display:inline-block;margin-right:7px}}
    .grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(390px,1fr));gap:20px}} .card{{overflow:hidden;background:var(--card);border:1px solid var(--line);border-radius:18px;box-shadow:0 18px 45px #0006}}
    img{{display:block;width:100%;aspect-ratio:4/3;object-fit:cover}} .copy{{padding:16px 18px 18px}} h2{{margin:0 0 8px}} h2 span{{float:right;color:var(--muted);font-size:16px;font-weight:500}} p{{margin:0;color:var(--muted);font-size:14px}}
    </style></head><body><main><h1>AFM α-Fe₂O₃ relaxed structures</h1><p class="sub">Five fixed-volume PBE+U relaxations · OVITO 3.10.5 · 2×2×2 periodic replicas · Fe–O bonds ≤ 2.30 Å</p>
    <div class="legend"><span><i class="dot" style="background:#f55a33"></i>Fe↑</span><span><i class="dot" style="background:#1ab8c7"></i>Fe↓</span><span><i class="dot" style="background:#dbe6f5"></i>O</span></div>
    <section class="grid">{''.join(cards)}</section></main></body></html>'''
    path = OUTPUT_DIR / "index.html"
    path.write_text(page, encoding="utf-8")
    return path


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    records = []
    for path in sorted(Path(".").glob("*.out")):
        try:
            records.append(parse_output(path))
        except ValueError:
            pass
    records.sort(key=lambda record: record["volume"])
    if len(records) != 5:
        raise SystemExit(f"Expected five completed structures; found {len(records)}")
    for record in records:
        extxyz = write_extxyz(record)
        record["image"] = render_record(record, extxyz)
        print(f"Rendered {record['label']} from {record['path']}")
    dashboard, pdf = build_dashboard(records)
    html_path = build_html(records)
    print(f"Dashboard: {dashboard}")
    print(f"PDF: {pdf}")
    print(f"HTML: {html_path}")


if __name__ == "__main__":
    main()
