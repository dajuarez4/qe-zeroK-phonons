#!/usr/bin/env python3
"""Create a side-view GIF from Quantum ESPRESSO MD output."""

from pathlib import Path
import argparse
import io
import re

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

BOHR_TO_ANG = 0.529177210903


def read_qe(path):
    lines = path.read_text(errors="replace").splitlines()
    alat = nat = None
    cell = None
    for i, line in enumerate(lines):
        if alat is None and "lattice parameter (alat)" in line:
            alat = float(line.split("=")[1].split()[0]) * BOHR_TO_ANG
        if nat is None and "number of atoms/cell" in line:
            nat = int(line.split("=")[1])
        if cell is None and "crystal axes:" in line:
            rows = []
            for row in lines[i + 1:i + 4]:
                match = re.search(r"\(\s*([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*\)", row)
                rows.append([float(match.group(j)) for j in range(1, 4)])
            cell = np.asarray(rows)
    if alat is None or nat is None or cell is None:
        raise ValueError("Could not read cell or atom count")
    cell *= alat

    steps, species, fractional = [], None, []
    current_step = 0
    for i, line in enumerate(lines):
        match = re.search(r"Entering Dynamics:\s+iteration\s*=\s*(\d+)", line)
        if match:
            current_step = int(match.group(1))
        if line.strip().upper() == "ATOMIC_POSITIONS (CRYSTAL)":
            symbols, coords = [], []
            for row in lines[i + 1:i + 1 + nat]:
                fields = row.split()
                if len(fields) < 4 or fields[0] not in ("Ga", "O"):
                    break
                symbols.append(fields[0])
                coords.append([float(v) for v in fields[1:4]])
            if len(coords) == nat:
                if species is None:
                    species = np.asarray(symbols)
                fractional.append(np.asarray(coords))
                steps.append(current_step)
    if not fractional:
        raise ValueError("No complete QE MD frames found")
    return cell, species, np.asarray(fractional), np.asarray(steps)


def unwrap_and_align(frac, cell, species):
    unwrapped = np.empty_like(frac)
    unwrapped[0] = frac[0]
    for frame in range(1, len(frac)):
        delta = frac[frame] - frac[frame - 1]
        delta -= np.round(delta)
        unwrapped[frame] = unwrapped[frame - 1] + delta
    cart = unwrapped @ cell
    masses = np.where(species == "Ga", 69.723, 15.999)
    com = (cart * masses[None, :, None]).sum(axis=1) / masses.sum()
    cart -= (com - com[0])[:, None, :]
    return cart


def project(pos):
    return np.column_stack((pos[:, 0] + 0.25 * pos[:, 1], pos[:, 2] + 0.06 * pos[:, 1]))


def bond_pairs(pos, species, cutoff=2.35):
    delta = pos[:, None, :] - pos[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", delta, delta)
    mask = (species[:, None] != species[None, :]) & (d2 > 0.25) & (d2 < cutoff**2)
    return np.argwhere(np.triu(mask, 1))


def render(cell, species, positions, steps, output, max_frames=60, fps=8):
    indices = np.unique(np.linspace(0, len(positions) - 1, min(max_frames, len(positions))).astype(int))
    positions = positions[indices]
    steps = steps[indices]
    projected = [project(pos) for pos in positions]
    extent = np.concatenate(projected)
    span = np.ptp(extent, axis=0)
    margin = np.maximum(span * 0.10, [1.8, 1.8])
    low, high = extent.min(axis=0) - margin, extent.max(axis=0) + margin

    images = []
    for number, (pos, xy, step) in enumerate(zip(positions, projected, steps), 1):
        pairs = bond_pairs(pos, species)
        fig, ax = plt.subplots(figsize=(9.5, 5.6), dpi=110, facecolor="#0b1020")
        ax.set_facecolor("#0b1020")
        for i, j in pairs:
            ax.plot(xy[[i, j], 0], xy[[i, j], 1], color="#aab3c5", lw=1.0, alpha=0.65, zorder=1)
        ga = species == "Ga"
        ax.scatter(xy[~ga, 0], xy[~ga, 1], s=58, c="#ef4444", edgecolors="#fecaca", linewidths=0.4, label="O", zorder=2)
        ax.scatter(xy[ga, 0], xy[ga, 1], s=105, c="#38bdf8", edgecolors="#e0f2fe", linewidths=0.5, label="Ga", zorder=3)
        ax.set_xlim(low[0], high[0])
        ax.set_ylim(low[1], high[1])
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("in-plane direction", color="#dbeafe")
        ax.set_ylabel("out-of-plane direction", color="#dbeafe")
        ax.tick_params(colors="#dbeafe")
        for spine in ax.spines.values():
            spine.set_color("#475569")
        ax.set_title("Ga$_2$O$_3$ bilayer — QE molecular dynamics at 300 K", color="white", fontsize=15, pad=12)
        ax.text(0.015, 0.975, f"frame {number}/{len(positions)}   MD step {step}\nactual displacement ×1", transform=ax.transAxes, va="top", color="white", fontsize=10, bbox=dict(boxstyle="round", fc="#111827", ec="#64748b", alpha=0.9))
        legend = ax.legend(loc="upper right", facecolor="#111827", edgecolor="#64748b")
        for label in legend.get_texts():
            label.set_color("white")
        fig.tight_layout()
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", facecolor=fig.get_facecolor())
        plt.close(fig)
        buffer.seek(0)
        images.append(Image.open(buffer).convert("P", palette=Image.Palette.ADAPTIVE))
    images[0].save(output, save_all=True, append_images=images[1:], duration=round(1000 / fps), loop=0, optimize=True, disposal=2)
    return len(indices)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("qe_output", type=Path)
    parser.add_argument("gif", type=Path)
    parser.add_argument("--max-frames", type=int, default=60)
    parser.add_argument("--fps", type=int, default=8)
    args = parser.parse_args()
    cell, species, frac, steps = read_qe(args.qe_output)
    positions = unwrap_and_align(frac, cell, species)
    used = render(cell, species, positions, steps, args.gif, args.max_frames, args.fps)
    print(f"Created {args.gif}: {used}/{len(frac)} frames, {len(species)} atoms, displacement x1")


if __name__ == "__main__":
    main()
