#!/usr/bin/env python3
"""Run the local No-Vito GIF renderer with slab-focused axis limits."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


NOVITO_RENDERER = Path(
    "/Users/dajuarez4/Documents/Fe/HELD_fcc/No-Vito/utils/qe_npz_to_gif.py"
)


def slab_focused_limits(ax, cell_ang: np.ndarray) -> None:
    """Focus on the bilayer and remove most of the vacuum from the view."""
    a, b, c = np.asarray(cell_ang, dtype=float)
    vertices = np.asarray(
        [
            np.zeros(3),
            a,
            b,
            c,
            a + b,
            a + c,
            b + c,
            a + b + c,
        ]
    )
    minima = vertices.min(axis=0)
    maxima = vertices.max(axis=0)

    xy_padding = 0.04 * np.max(maxima[:2] - minima[:2])
    x_limits = (minima[0] - xy_padding, maxima[0] + xy_padding)
    y_limits = (minima[1] - xy_padding, maxima[1] + xy_padding)

    # The bilayer is centered in the 30-Angstrom out-of-plane cell.
    z_center = 0.5 * (minima[2] + maxima[2])
    z_half_span = 5.3
    z_limits = (z_center - z_half_span, z_center + z_half_span)

    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_zlim(*z_limits)
    ax.set_box_aspect(
        (
            x_limits[1] - x_limits[0],
            y_limits[1] - y_limits[0],
            z_limits[1] - z_limits[0],
        )
    )


def main() -> None:
    spec = importlib.util.spec_from_file_location("novito_qe_renderer", NOVITO_RENDERER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load No-Vito renderer: {NOVITO_RENDERER}")
    renderer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(renderer)
    renderer.set_axes_limits = slab_focused_limits
    renderer.main()


if __name__ == "__main__":
    main()
