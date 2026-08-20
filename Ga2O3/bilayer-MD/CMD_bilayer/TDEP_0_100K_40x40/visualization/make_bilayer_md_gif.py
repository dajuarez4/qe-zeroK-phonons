#!/usr/bin/env python3
"""Render a zoomed Ga2O3 bilayer MD trajectory as an animated GIF."""

from pathlib import Path
import argparse
import io

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def read_dump(path):
    frames = []
    with open(path) as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(handle.readline())
            handle.readline()  # NUMBER OF ATOMS
            natoms = int(handle.readline())
            handle.readline()  # BOX BOUNDS
            bounds = [list(map(float, handle.readline().split())) for _ in range(3)]
            columns = handle.readline().split()[2:]
            col = {name: i for i, name in enumerate(columns)}
            atoms = np.array(
                [list(map(float, handle.readline().split())) for _ in range(natoms)]
            )
            order = np.argsort(atoms[:, col["id"]].astype(int))
            atoms = atoms[order]
            pos = atoms[:, [col["xu"], col["yu"], col["zu"]]]
            types = atoms[:, col["type"]].astype(int)
            frames.append((step, bounds, types, pos))
    return frames


def project(pos):
    # Oblique side view: x is horizontal, z separates the layers, y adds depth.
    return np.column_stack((pos[:, 0] + 0.28 * pos[:, 1], pos[:, 2] + 0.08 * pos[:, 1]))


def bonds(pos, types, cutoff=2.35):
    # The crop is small enough for a direct pair search; only Ga-O bonds are drawn.
    delta = pos[:, None, :] - pos[None, :, :]
    dist2 = np.einsum("ijk,ijk->ij", delta, delta)
    mask = (types[:, None] != types[None, :]) & (dist2 < cutoff**2) & (dist2 > 0.3**2)
    return np.argwhere(np.triu(mask, 1))


def render(frames, output, amplify=5.0, window=24.0, fps=6):
    ref = frames[0][3]
    xmid, ymid = np.median(ref[:, 0]), np.median(ref[:, 1])
    keep = (
        (np.abs(ref[:, 0] - xmid) <= window / 2)
        & (np.abs(ref[:, 1] - ymid) <= window / 2)
    )
    ref = ref[keep]
    types = frames[0][2][keep]
    pair_list = bonds(ref, types)

    rendered = []
    limits = project(ref)
    xmin, ymin = limits.min(axis=0) - (2.0, 1.8)
    xmax, ymax = limits.max(axis=0) + (2.0, 1.8)
    for index, (step, _bounds, _types, current_all) in enumerate(frames):
        current = current_all[keep]
        shown = ref + amplify * (current - ref)
        xy = project(shown)
        fig, ax = plt.subplots(figsize=(9.5, 5.6), dpi=110, facecolor="#0b1020")
        ax.set_facecolor("#0b1020")
        for i, j in pair_list:
            ax.plot(xy[[i, j], 0], xy[[i, j], 1], color="#aab3c5", lw=0.7, alpha=0.55, zorder=1)
        ga = types == 1
        oxygen = ~ga
        ax.scatter(xy[oxygen, 0], xy[oxygen, 1], s=42, c="#ef4444", edgecolors="#fecaca", linewidths=0.35, label="O", zorder=2)
        ax.scatter(xy[ga, 0], xy[ga, 1], s=82, c="#38bdf8", edgecolors="#e0f2fe", linewidths=0.45, label="Ga", zorder=3)
        ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), xlabel="in-plane direction", ylabel="out-of-plane direction")
        ax.set_aspect("equal", adjustable="box")
        ax.tick_params(colors="#dbeafe")
        for spine in ax.spines.values():
            spine.set_color("#475569")
        ax.xaxis.label.set_color("#dbeafe")
        ax.yaxis.label.set_color("#dbeafe")
        ax.set_title("Ga$_2$O$_3$ bilayer MD — 50 K", color="white", fontsize=15, pad=12)
        ax.text(0.015, 0.975, f"frame {index + 1}/{len(frames)}   timestep {step}\ndisplacements ×{amplify:g}", transform=ax.transAxes, va="top", color="white", fontsize=10, bbox=dict(boxstyle="round", fc="#111827", ec="#64748b", alpha=0.9))
        legend = ax.legend(loc="upper right", frameon=True, facecolor="#111827", edgecolor="#64748b")
        for label in legend.get_texts():
            label.set_color("white")
        fig.tight_layout()
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", facecolor=fig.get_facecolor())
        plt.close(fig)
        buffer.seek(0)
        rendered.append(Image.open(buffer).convert("P", palette=Image.Palette.ADAPTIVE))

    duration = round(1000 / fps)
    rendered[0].save(output, save_all=True, append_images=rendered[1:], duration=duration, loop=0, optimize=True, disposal=2)
    return keep.sum(), len(frames)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--amplify", type=float, default=5.0)
    parser.add_argument("--window", type=float, default=24.0)
    parser.add_argument("--fps", type=int, default=6)
    args = parser.parse_args()
    frames = read_dump(args.trajectory)
    if not frames:
        raise SystemExit("No frames found")
    atoms, count = render(frames, args.output, args.amplify, args.window, args.fps)
    print(f"Created {args.output} from {count} frames; showing {atoms} atoms")


if __name__ == "__main__":
    main()
