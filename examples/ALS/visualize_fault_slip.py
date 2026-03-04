#!/usr/bin/env python3
"""Visualize the location, magnitude, and direction of slip from a PyLith SpatialDB file.

Usage:
    python visualize_fault_slip.py [path_to_spatialdb]

If no path is given, defaults to the fault_1_slip.spatialdb in the same example directory.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from pathlib import Path


def parse_spatialdb(filepath):
    """Parse a PyLith SimpleDB spatial database file.

    Returns
    -------
    coords : ndarray, shape (N, 3)
        (x, y, z) coordinates of each slip location.
    values : dict
        Mapping of value-name -> 1-D array of length N.
    metadata : dict
        Header metadata (value-names, value-units, crs-string, etc.).
    """
    with open(filepath) as f:
        lines = f.readlines()

    # ---- parse header ----
    metadata = {}
    header_end = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("num-values"):
            metadata["num-values"] = int(stripped.split("=")[1])
        elif stripped.startswith("value-names"):
            metadata["value-names"] = stripped.split("=")[1].split()
        elif stripped.startswith("value-units"):
            metadata["value-units"] = stripped.split("=")[1].split()
        elif stripped.startswith("num-locs"):
            metadata["num-locs"] = int(stripped.split("=")[1])
        elif stripped.startswith("data-dim"):
            metadata["data-dim"] = int(stripped.split("=")[1])
        elif stripped.startswith("space-dim"):
            metadata["space-dim"] = int(stripped.split("=")[1])
        elif stripped.startswith("crs-string"):
            metadata["crs-string"] = stripped.split("=")[1].strip()
        if stripped == "}":
            # Check if this closes the top-level SimpleDB block
            # (the last closing brace before data begins)
            header_end = i + 1

    # ---- parse data ----
    num_values = metadata["num-values"]
    space_dim = metadata["space-dim"]
    value_names = metadata["value-names"]

    data_lines = [l.strip() for l in lines[header_end:] if l.strip()]
    raw = np.array([[float(v) for v in l.split()] for l in data_lines])

    coords = raw[:, :space_dim]
    value_arrays = {name: raw[:, space_dim + j] for j, name in enumerate(value_names)}

    return coords, value_arrays, metadata


def main(filepath=None):
    if filepath is None:
        filepath = (
            Path(__file__).resolve().parent
            / "20260206_190857_pylith_inputs"
            / "faults"
            / "fault_1_slip.spatialdb"
        )
    filepath = Path(filepath)
    print(f"Reading: {filepath}")

    coords, values, meta = parse_spatialdb(filepath)
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]

    ll = values.get("final_slip_left_lateral", np.zeros(len(x)))
    rev = values.get("final_slip_reverse", np.zeros(len(x)))
    opening = values.get("final_slip_opening", np.zeros(len(x)))

    # Total slip magnitude (Euclidean norm of the three components)
    magnitude = np.sqrt(ll**2 + rev**2 + opening**2)

    # ---- Estimate local fault-strike direction from neighbours ----
    # Sort points by distance along the fault trace (approximate with a
    # nearest-neighbour chain starting from one end).
    n = len(x)
    visited = np.zeros(n, dtype=bool)
    order = [0]
    visited[0] = True
    for _ in range(n - 1):
        last = order[-1]
        dists = (x - x[last])**2 + (y - y[last])**2
        dists[visited] = np.inf
        nearest = np.argmin(dists)
        order.append(nearest)
        visited[nearest] = True
    order = np.array(order)

    # Reorder all arrays by the chain order
    x_o, y_o = x[order], y[order]
    ll_o, rev_o, opening_o = ll[order], rev[order], opening[order]
    mag_o = magnitude[order]

    # Local tangent (strike) direction via finite differences
    dx = np.gradient(x_o)
    dy = np.gradient(y_o)
    length = np.sqrt(dx**2 + dy**2)
    length[length == 0] = 1.0
    tx, ty = dx / length, dy / length  # unit tangent along strike

    # Normal direction (90° CCW rotation of tangent) – approximates dip-direction
    # projected onto the surface.  "reverse" slip acts in the normal direction.
    nx, ny = -ty, tx

    # Build slip vector in map coordinates:
    #   left-lateral  ➜ along-strike (tangent direction)
    #   reverse       ➜ along fault-normal (up-dip direction projected)
    #   opening       ➜ along fault-normal (tensile)
    slip_x = ll_o * tx + (rev_o + opening_o) * nx
    slip_y = ll_o * ty + (rev_o + opening_o) * ny

    # ---- Figure ----
    fig, ax = plt.subplots(figsize=(12, 10))

    # Compute arrow scale: make arrows a visible fraction of the map extent
    extent = max(x_o.max() - x_o.min(), y_o.max() - y_o.min())
    max_slip_vec = np.max(np.sqrt(slip_x**2 + slip_y**2))
    if max_slip_vec > 0:
        arrow_scale = 0.025 * extent / max_slip_vec
    else:
        arrow_scale = 1.0

    # Color each arrow by slip magnitude
    norm = Normalize(vmin=mag_o.min(), vmax=max(mag_o.max(), 1e-12))
    cmap = plt.colormaps["plasma"]
    colors = cmap(norm(mag_o))

    # Plot each point as a dot + individual arrow
    for i in range(len(x_o)):
        dx = slip_x[i] * arrow_scale
        dy = slip_y[i] * arrow_scale
        c = colors[i]
        # Point
        ax.plot(x_o[i], y_o[i], "o", color=c, markersize=4, markeredgecolor="k",
                markeredgewidth=0.3, zorder=3)
        # Arrow showing slip direction & relative magnitude
        if mag_o[i] > 0:
            ax.annotate(
                "", xy=(x_o[i] + dx, y_o[i] + dy), xytext=(x_o[i], y_o[i]),
                arrowprops=dict(arrowstyle="-|>", color=c, lw=1.2),
                zorder=4,
            )

    # Colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Slip magnitude (m)")

    ax.set_xlabel(f"Easting (m)  [{meta.get('crs-string', '')}]")
    ax.set_ylabel("Northing (m)")
    ax.set_title(f"Fault 1 Slip — {len(x)} locations\n"
                 f"(left-lateral / reverse / opening)", fontsize=13, fontweight="bold")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    out_png = filepath.with_suffix(".png")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Saved figure to: {out_png}")
    plt.show()


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else None
    main(path)
