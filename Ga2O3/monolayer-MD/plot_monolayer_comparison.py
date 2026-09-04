#!/usr/bin/env python3
"""Plot a publication-style comparison of the ML-001 and ML-100 structures."""

from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
import numpy as np


ROOT = Path(__file__).resolve().parent
STRUCTURES = [
    ("ML-001", ROOT / "fdf_stru" / "Ga2O3-GGA-Bands-ml-001.fdf"),
    ("ML-100", ROOT / "fdf_stru" / "Ga2O3-GGA-Bands-ml-100.fdf"),
]
COLORS = {"Ga": "#8c62c7", "O": "#e94b4b"}
RADII = {"Ga": 0.34, "O": 0.25}


def read_fdf(path):
    text = path.read_text()

    def block(name):
        match = re.search(
            rf"%block\s+{name}\s*\n(.*?)%endblock\s+{name}", text, re.I | re.S
        )
        return [line.split() for line in match.group(1).splitlines() if line.strip()]

    cell = np.array([[float(x) for x in row[:3]] for row in block("LatticeVectors")])
    atoms = []
    for row in block("AtomicCoordinatesAndAtomicSpecies"):
        symbol = "O" if int(row[3]) == 1 else "Ga"
        frac = np.array([float(x) for x in row[:3]])
        atoms.append((symbol, frac @ cell))
    return cell, atoms


def hex_to_rgb(value):
    value = value.lstrip("#")
    return np.array([int(value[i : i + 2], 16) for i in (0, 2, 4)]) / 255


def sphere(ax, x, y, symbol, scale=1.0, zorder=5):
    radius = RADII[symbol] * scale
    base = hex_to_rgb(COLORS[symbol])
    ax.add_patch(Circle((x + 0.10 * radius, y - 0.13 * radius), 1.06 * radius,
                        color="black", alpha=0.22, linewidth=0, zorder=zorder - 2))
    for i in range(18, 0, -1):
        f = i / 18
        shade = base * (0.58 + 0.42 * (1 - f))
        ax.add_patch(Circle((x, y), radius * f, color=shade, linewidth=0,
                            zorder=zorder - 1 + (1 - f) * 0.01))
    ax.add_patch(Circle((x - 0.27 * radius, y + 0.30 * radius), 0.25 * radius,
                        color="white", alpha=0.52, linewidth=0, zorder=zorder + 1))
    ax.add_patch(Circle((x, y), radius, fill=False, edgecolor="white",
                        linewidth=0.45, alpha=0.45, zorder=zorder + 2))


def replicated_atoms(cell, atoms, repeats=(3, 2)):
    result = []
    for ia in range(repeats[0]):
        for ib in range(repeats[1]):
            shift = ia * cell[0] + ib * cell[1]
            result.extend((symbol, xyz + shift) for symbol, xyz in atoms)
    return result


def bonds(atoms, cutoff=2.25):
    pairs = []
    for i, (si, ri) in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            sj, rj = atoms[j]
            if si == sj:
                continue
            distance = np.linalg.norm(ri - rj)
            if distance < cutoff:
                pairs.append((ri, rj))
    return pairs


def style_axis(ax):
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor("#f6f7fa")


def plot_top(ax, cell, base_atoms, title):
    atoms = replicated_atoms(cell, base_atoms)
    for ri, rj in bonds(atoms):
        ax.plot([ri[0], rj[0]], [ri[1], rj[1]], color="#667085", lw=2.0,
                alpha=0.48, solid_capstyle="round", zorder=1)
    outline = np.array([[0, 0], cell[0, :2], cell[0, :2] + cell[1, :2], cell[1, :2]])
    ax.add_patch(Polygon(outline, closed=True, fill=False, edgecolor="#168aad",
                         linewidth=2.0, linestyle=(0, (4, 2)), zorder=3))
    for symbol, xyz in sorted(atoms, key=lambda item: item[1][2]):
        sphere(ax, xyz[0], xyz[1], symbol)
    xy = np.array([xyz[:2] for _, xyz in atoms])
    ax.set_xlim(xy[:, 0].min() - 0.8, xy[:, 0].max() + 0.8)
    ax.set_ylim(xy[:, 1].min() - 0.8, xy[:, 1].max() + 0.8)
    ax.set_title(f"{title} · vista superior", fontsize=15, fontweight="bold", pad=12)
    style_axis(ax)


def plot_side(ax, cell, base_atoms, title):
    atoms = replicated_atoms(cell, base_atoms, repeats=(3, 1))
    # Use distance along a as the horizontal coordinate and z relative to layer center.
    ahat = cell[0] / np.linalg.norm(cell[0])
    projected = [(symbol, np.dot(xyz, ahat), xyz[2]) for symbol, xyz in atoms]
    center_z = 0.5 * (min(z for _, _, z in projected) + max(z for _, _, z in projected))
    projected = [(s, x, z - center_z) for s, x, z in projected]
    for i, (si, xi, zi) in enumerate(projected):
        ri = atoms[i][1]
        for j in range(i + 1, len(projected)):
            sj, xj, zj = projected[j]
            if si != sj and np.linalg.norm(ri - atoms[j][1]) < 2.25:
                ax.plot([xi, xj], [zi, zj], color="#667085", lw=2.0,
                        alpha=0.48, solid_capstyle="round", zorder=1)
    for symbol, x, z in sorted(projected, key=lambda item: item[2]):
        sphere(ax, x, z, symbol, scale=0.92)
    xs = [x for _, x, _ in projected]
    zs = [z for _, _, z in projected]
    ax.set_xlim(min(xs) - 0.8, max(xs) + 0.8)
    ax.set_ylim(min(zs) - 0.8, max(zs) + 0.8)
    ax.axhline(0, color="#168aad", lw=1.0, alpha=0.28, zorder=0)
    ax.set_title(f"{title} · vista lateral", fontsize=13, fontweight="bold", pad=9)
    style_axis(ax)


def main():
    loaded = [(name, *read_fdf(path)) for name, path in STRUCTURES]
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.0), facecolor="white")
    for col, (name, cell, atoms) in enumerate(loaded):
        plot_top(axes[0, col], cell, atoms, name)
        plot_side(axes[1, col], cell, atoms, name)

    handles = []
    for symbol, label in [("Ga", "Galio (Ga)"), ("O", "Oxígeno (O)")]:
        handles.append(plt.Line2D([], [], marker="o", linestyle="", markersize=11,
                                  markerfacecolor=COLORS[symbol], markeredgecolor="white",
                                  markeredgewidth=0.8, label=label))
    handles.append(plt.Line2D([], [], color="#168aad", lw=2, linestyle=(0, (4, 2)),
                              label="Celda primitiva"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=11,
               bbox_to_anchor=(0.5, 0.018))
    fig.suptitle("Comparación estructural de monocapas de Ga₂O₃", fontsize=20,
                 fontweight="bold", y=0.985)
    fig.text(0.5, 0.935, "Geometrías relajadas de los archivos FDF originales",
             ha="center", color="#667085", fontsize=11)
    plt.subplots_adjust(left=0.035, right=0.98, top=0.875, bottom=0.09, wspace=0.08, hspace=0.20)

    png = ROOT / "Ga2O3_monolayer_structures_comparison.png"
    pdf = ROOT / "Ga2O3_monolayer_structures_comparison.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")


if __name__ == "__main__":
    main()
