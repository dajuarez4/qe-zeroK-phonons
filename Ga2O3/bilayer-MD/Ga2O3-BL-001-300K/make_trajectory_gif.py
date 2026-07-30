#!/usr/bin/env python3
"""Render complete QE MD frames as a dependency-free animated GIF."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import tempfile
from pathlib import Path

WIDTH = 900
HEIGHT = 650
MARGIN = 65
BACKGROUND = (248, 249, 251)
CELL_COLOR = (100, 108, 118)
BOND_COLOR = (145, 150, 158)
ATOM_STYLE = {
    "Ga": ((65, 105, 190), 15),
    "O": ((220, 70, 62), 10),
}


def parse_input(path: Path):
    lines = path.read_text().splitlines()

    cell_start = next(i for i, line in enumerate(lines) if line.startswith("CELL_PARAMETERS")) + 1
    cell = [
        [float(value) for value in lines[cell_start + offset].split()[:3]]
        for offset in range(3)
    ]

    pos_start = next(i for i, line in enumerate(lines) if line.startswith("ATOMIC_POSITIONS")) + 1
    symbols = []
    for line in lines[pos_start:]:
        fields = line.split()
        if len(fields) < 4 or fields[0] in {"K_POINTS", "CELL_PARAMETERS"}:
            break
        symbols.append(fields[0])
    return cell, symbols


def parse_complete_frames(path: Path, nat: int):
    text = path.read_text(errors="replace")
    pattern = re.compile(
        r"^ATOMIC_POSITIONS\s+\(crystal\)\s*$((?:\n\s*[A-Za-z]+\s+"
        r"[-+0-9.Ee]+\s+[-+0-9.Ee]+\s+[-+0-9.Ee]+[^\n]*){" + str(nat) + r"})",
        re.MULTILINE,
    )
    frames = []
    for match in pattern.finditer(text):
        frame = []
        for line in match.group(1).strip().splitlines():
            fields = line.split()
            frame.append([float(value) for value in fields[1:4]])
        frames.append(frame)

    force_blocks = len(re.findall(r"^\s*Forces acting on atoms .*:$", text, re.MULTILINE))
    # The first force block is the pre-MD force evaluation.
    complete = min(len(frames), max(0, force_blocks - 1))
    if complete == 0:
        raise ValueError("No complete QE MD position/force pairs were found")
    return frames[:complete], len(frames) - complete


def cartesian(frac, cell):
    return [
        sum(frac[k] * cell[k][axis] for k in range(3))
        for axis in range(3)
    ]


def projection_basis(point, zcenter):
    angle = math.radians(25.0)
    x, y, z = point
    u = x * math.cos(angle) + y * math.sin(angle)
    depth = -x * math.sin(angle) + y * math.cos(angle)
    v = 1.05 * (z - zcenter) + 0.36 * depth
    return u, v, depth


def put_pixel(image, x, y, color):
    if 0 <= x < WIDTH and 0 <= y < HEIGHT:
        index = 3 * (y * WIDTH + x)
        image[index : index + 3] = bytes(color)


def circle(image, cx, cy, radius, color, outline=None):
    outer = radius + (2 if outline else 0)
    for dy in range(-outer, outer + 1):
        span = int(math.sqrt(max(0, outer * outer - dy * dy)))
        for dx in range(-span, span + 1):
            distance2 = dx * dx + dy * dy
            pixel_color = outline if outline and distance2 > radius * radius else color
            put_pixel(image, cx + dx, cy + dy, pixel_color)


def line(image, x1, y1, x2, y2, color, width=2):
    steps = max(abs(x2 - x1), abs(y2 - y1), 1)
    for step in range(steps + 1):
        fraction = step / steps
        x = round(x1 + fraction * (x2 - x1))
        y = round(y1 + fraction * (y2 - y1))
        circle(image, x, y, width, color)


def write_ppm(path: Path, image):
    with path.open("wb") as handle:
        handle.write(f"P6\n{WIDTH} {HEIGHT}\n255\n".encode())
        handle.write(image)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("Ga2O3-BL-001.md.in"))
    parser.add_argument("--output", type=Path, default=Path("Ga2O3-BL-001.md.out"))
    parser.add_argument("--gif", type=Path, default=Path("Ga2O3-BL-001-trajectory.gif"))
    parser.add_argument("--fps", type=int, default=4)
    args = parser.parse_args()

    cell, symbols = parse_input(args.input)
    frames, excluded = parse_complete_frames(args.output, len(symbols))
    cart_frames = [[cartesian(position, cell) for position in frame] for frame in frames]

    all_points = [point for frame in cart_frames for point in frame]
    zcenter = sum(point[2] for point in all_points) / len(all_points)
    cell_corners = [
        [0.0, 0.0, zcenter],
        [*cell[0][:2], zcenter],
        [cell[0][0] + cell[1][0], cell[0][1] + cell[1][1], zcenter],
        [*cell[1][:2], zcenter],
    ]
    projected_all = [
        projection_basis(point, zcenter)
        for point in all_points + cell_corners
    ]
    umin = min(item[0] for item in projected_all)
    umax = max(item[0] for item in projected_all)
    vmin = min(item[1] for item in projected_all)
    vmax = max(item[1] for item in projected_all)
    scale = min(
        (WIDTH - 2 * MARGIN) / max(umax - umin, 1e-9),
        (HEIGHT - 2 * MARGIN) / max(vmax - vmin, 1e-9),
    )

    def screen(point):
        u, v, depth = projection_basis(point, zcenter)
        x = round(MARGIN + (u - umin) * scale)
        y = round(HEIGHT - MARGIN - (v - vmin) * scale)
        return x, y, depth

    cell_screen = [screen(point) for point in cell_corners]

    with tempfile.TemporaryDirectory(prefix="ga2o3_trajectory_") as temporary:
        temporary_path = Path(temporary)
        for frame_number, points in enumerate(cart_frames):
            image = bytearray(BACKGROUND * (WIDTH * HEIGHT))

            for first, second in zip(cell_screen, cell_screen[1:] + cell_screen[:1]):
                line(image, first[0], first[1], second[0], second[1], CELL_COLOR, width=1)

            projected = [screen(point) for point in points]
            bonds = []
            for i in range(len(points)):
                for j in range(i + 1, len(points)):
                    if symbols[i] == symbols[j]:
                        continue
                    distance = math.sqrt(
                        sum((points[i][axis] - points[j][axis]) ** 2 for axis in range(3))
                    )
                    if distance <= 2.25:
                        bonds.append((i, j, (projected[i][2] + projected[j][2]) / 2))
            for i, j, _ in sorted(bonds, key=lambda item: item[2], reverse=True):
                line(
                    image,
                    projected[i][0],
                    projected[i][1],
                    projected[j][0],
                    projected[j][1],
                    BOND_COLOR,
                    width=2,
                )

            atom_order = sorted(range(len(points)), key=lambda i: projected[i][2], reverse=True)
            for index in atom_order:
                color, radius = ATOM_STYLE[symbols[index]]
                circle(
                    image,
                    projected[index][0],
                    projected[index][1],
                    radius,
                    color,
                    outline=(35, 39, 46),
                )

            write_ppm(temporary_path / f"frame_{frame_number:04d}.ppm", image)

        subprocess.run(
            [
                "ffmpeg",
                "-y",
                "-loglevel",
                "error",
                "-framerate",
                str(args.fps),
                "-i",
                str(temporary_path / "frame_%04d.ppm"),
                "-filter_complex",
                "[0:v]split[a][b];"
                "[a]palettegen=max_colors=128:stats_mode=diff[p];"
                "[b][p]paletteuse=dither=sierra2_4a",
                "-loop",
                "0",
                str(args.gif),
            ],
            check=True,
        )

    print(f"Wrote {args.gif} from {len(frames)} complete MD frames")
    print(f"Unmatched position frames excluded: {excluded}")


if __name__ == "__main__":
    main()
