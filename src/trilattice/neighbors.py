"""Neighbour finding: linked cell lists plus a Verlet list with a skin.

The original code built eight periodic replicas of the whole system and then a
``9N x 9N`` distance matrix -- 81 N^2 distances per step, of which ~0.03% were
inside the cutoff.  At N = 50 that is merely wasteful; at N = 1000 it is 81
million entries per step and the run never finishes.

The standard fix is two-layered:

1. **Cell list.**  Partition the box into cells of edge >= r_list.  A particle
   can only interact with particles in its own cell and the eight neighbours, so
   the pair search is O(N) with a prefactor set by the occupancy.
2. **Verlet list.**  Cache the pairs within ``r_list = r_cut + skin`` and reuse
   them until some particle has moved more than ``skin/2`` (at which point a pair
   that was outside ``r_list`` could have entered ``r_cut``).  This amortises the
   search over tens of steps.

Only the "half" stencil is enumerated, so every pair appears exactly once and
Newton's third law can be used to halve the force work.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .lattice import Box

#: Half stencil: self-cell plus four of the eight neighbours.  Combined with the
#: ``i < j`` filter inside the self-cell this visits each pair exactly once.
_HALF_STENCIL = ((0, 0), (1, 0), (-1, 1), (0, 1), (1, 1))


@dataclass
class NeighborList:
    """Verlet pair list for an orthorhombic periodic box."""

    box: Box
    cutoff: float
    skin: float = 1.0

    def __post_init__(self) -> None:
        if self.cutoff <= 0.0:
            raise ValueError("cutoff must be positive")
        if self.skin < 0.0:
            raise ValueError("skin must be non-negative")
        if self.r_list > self.box.max_cutoff:
            raise ValueError(
                f"r_list = {self.r_list:.3f} A exceeds half the smallest box edge "
                f"({self.box.max_cutoff:.3f} A); the minimum-image convention would "
                "be ambiguous.  Use a bigger supercell or a smaller cutoff."
            )
        self._i = np.zeros(0, dtype=np.int64)
        self._j = np.zeros(0, dtype=np.int64)
        self._reference = np.zeros((0, 2), dtype=np.float64)
        self.n_builds = 0
        self.n_queries = 0

    @property
    def r_list(self) -> float:
        return self.cutoff + self.skin

    @property
    def pairs(self) -> tuple[np.ndarray, np.ndarray]:
        return self._i, self._j

    @property
    def n_pairs(self) -> int:
        return self._i.size

    # ------------------------------------------------------------------ #
    def needs_rebuild(self, positions: np.ndarray) -> bool:
        """True once any particle has drifted more than half the skin."""
        if self._reference.shape != positions.shape:
            return True
        d = self.box.minimum_image(positions - self._reference)
        return bool(np.max(d[:, 0] ** 2 + d[:, 1] ** 2) > (0.5 * self.skin) ** 2)

    def update(self, positions: np.ndarray, force: bool = False) -> bool:
        """Rebuild the list if necessary.  Returns True if a rebuild happened."""
        self.n_queries += 1
        if force or self.needs_rebuild(positions):
            self.build(positions)
            return True
        return False

    # ------------------------------------------------------------------ #
    def build(self, positions: np.ndarray) -> None:
        wrapped = self.box.wrap(positions)
        i, j = self._cell_list_pairs(wrapped)

        if i.size:
            d = self.box.minimum_image(wrapped[i] - wrapped[j])
            r2 = d[:, 0] ** 2 + d[:, 1] ** 2
            inside = r2 < self.r_list**2
            i, j = i[inside], j[inside]

        self._i = np.ascontiguousarray(i)
        self._j = np.ascontiguousarray(j)
        self._reference = positions.copy()
        self.n_builds += 1

    # ------------------------------------------------------------------ #
    def _cell_list_pairs(self, wrapped: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Candidate pairs from a padded, vectorised linked-cell decomposition.

        Particles are bucketed into a ``(ncx, ncy)`` grid stored as a dense
        ``(ncells, max_occupancy)`` index table padded with -1.  The pair
        enumeration is then a broadcast outer product per stencil offset: no
        Python-level loop over particles anywhere.
        """
        n = wrapped.shape[0]
        if n < 2:
            return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)

        ncx = max(1, int(self.box.lx / self.r_list))
        ncy = max(1, int(self.box.ly / self.r_list))
        # With fewer than 3 cells per direction the half stencil would wrap onto
        # itself and double-count; fall back to the direct O(N^2) enumeration.
        if ncx < 3 or ncy < 3:
            return self._all_pairs(n)

        cx = np.minimum((wrapped[:, 0] / self.box.lx * ncx).astype(np.int64), ncx - 1)
        cy = np.minimum((wrapped[:, 1] / self.box.ly * ncy).astype(np.int64), ncy - 1)
        cell = cx * ncy + cy
        ncells = ncx * ncy

        counts = np.bincount(cell, minlength=ncells)
        max_occ = int(counts.max())
        order = np.argsort(cell, kind="stable")
        slot = np.arange(n) - np.repeat(np.concatenate(([0], np.cumsum(counts)[:-1])), counts)

        table = np.full((ncells, max_occ), -1, dtype=np.int64)
        table[cell[order], slot] = order

        cx_grid, cy_grid = np.divmod(np.arange(ncells), ncy)

        chunks_i: list[np.ndarray] = []
        chunks_j: list[np.ndarray] = []
        for dx, dy in _HALF_STENCIL:
            nbr = ((cx_grid + dx) % ncx) * ncy + (cy_grid + dy) % ncy
            left = table[:, :, None]                      # (ncells, occ, 1)
            right = table[nbr][:, None, :]                # (ncells, 1, occ)
            a = np.broadcast_to(left, (ncells, max_occ, max_occ))
            b = np.broadcast_to(right, (ncells, max_occ, max_occ))
            if dx == 0 and dy == 0:
                keep = (a >= 0) & (b >= 0) & (a < b)
            else:
                keep = (a >= 0) & (b >= 0)
            chunks_i.append(a[keep])
            chunks_j.append(b[keep])

        return np.concatenate(chunks_i), np.concatenate(chunks_j)

    @staticmethod
    def _all_pairs(n: int) -> tuple[np.ndarray, np.ndarray]:
        i, j = np.triu_indices(n, k=1)
        return i.astype(np.int64), j.astype(np.int64)

    def stats(self) -> dict:
        return {
            "pairs": self.n_pairs,
            "builds": self.n_builds,
            "queries": self.n_queries,
            "rebuild_fraction": self.n_builds / max(1, self.n_queries),
            "r_list": self.r_list,
        }


def all_pairs_within(box: Box, positions: np.ndarray, cutoff: float):
    """Reference O(N^2) pair search -- used only to validate the cell list."""
    n = positions.shape[0]
    i, j = np.triu_indices(n, k=1)
    d = box.minimum_image(positions[i] - positions[j])
    r2 = d[:, 0] ** 2 + d[:, 1] ** 2
    m = r2 < cutoff**2
    return i[m], j[m]


__all__ = ["NeighborList", "all_pairs_within"]
