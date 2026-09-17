"""Matplotlib helpers with a consistent, publication-oriented style."""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from .lattice import Box

STYLE = {
    "figure.dpi": 130,
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "lines.linewidth": 1.4,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "figure.facecolor": "white",
}

#: Qualitative palette (Okabe-Ito, colour-blind safe).
COLORS = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#000000",
]


def use_style() -> None:
    mpl.rcParams.update(STYLE)


def plot_configuration(
    ax,
    positions: np.ndarray,
    box: Box,
    *,
    color=None,
    cmap="viridis",
    vmin=None,
    vmax=None,
    size=18.0,
    title=None,
    cbar_label=None,
    show_box=True,
):
    """Scatter a 2-D configuration with the periodic cell drawn."""
    sc = ax.scatter(
        positions[:, 0],
        positions[:, 1],
        c=color if color is not None else COLORS[0],
        cmap=cmap if color is not None else None,
        vmin=vmin,
        vmax=vmax,
        s=size,
        linewidths=0.3,
        edgecolors="0.25",
    )
    if show_box:
        ax.add_patch(
            plt.Rectangle((0, 0), box.lx, box.ly, fill=False, ec="0.4", lw=0.8, ls="--")
        )
    ax.set_xlim(-0.03 * box.lx, 1.03 * box.lx)
    ax.set_ylim(-0.03 * box.ly, 1.03 * box.ly)
    ax.set_aspect("equal")
    ax.set_xlabel("x (A)")
    ax.set_ylabel("y (A)")
    ax.grid(False)
    if title:
        ax.set_title(title)
    if color is not None and cbar_label:
        cb = ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(cbar_label)
    return sc


def annotate_shells(ax, radii, counts, y=None, max_shells=6):
    """Mark ideal-lattice coordination shells on a g(r) plot."""
    y = ax.get_ylim()[1] * 0.95 if y is None else y
    for r, c in list(zip(radii, counts))[:max_shells]:
        ax.axvline(r, color="0.7", lw=0.6, ls=":", zorder=0)
        ax.text(r, y, str(int(c)), ha="center", va="top", fontsize=6.5, color="0.45")


__all__ = ["use_style", "STYLE", "COLORS", "plot_configuration", "annotate_shells"]
