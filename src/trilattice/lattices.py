"""The lattice types the code can build, and what each one implies.

Every structure here is generated inside an orthorhombic cell whose edges are
integer multiples of the structure's own periods, so the crystal is never
strained by its own periodic images.

===========  =====  =========  =========  ====================================
kind         z      psi order  rho a^2    remark
===========  =====  =========  =========  ====================================
triangular   6      6          1.1547     close packed
square       4      4          1.0000     metastable
honeycomb    3      6          0.7698     graphene topology
kagome       4      6          0.8660     corner-sharing triangles
===========  =====  =========  =========  ====================================

**Only the triangular lattice is a mechanically stable ground state of an
isotropic pair potential.**  Close packing wins in 2-D for the same reason it
does in 3-D: with a central force there is nothing to pay for bond angles, so
the energy is minimised by maximising the number of neighbours at the potential
minimum.  Square, honeycomb and kagome are all *metastable at best* under
Lennard-Jones -- they survive at low temperature because there is a barrier to
rearranging, and collapse towards triangular once there is enough thermal energy
to cross it.  Measured with the default parameters (eps = 0.27 eV, sigma = 2.88 A,
a = 3.3567 A), the highest temperature at which each survives 5 ps is

===========  ==============  ===========================================
kind         survives to     what happens above it
===========  ==============  ===========================================
triangular   --              melts near 2400 K
kagome       ~100 K          z goes 4 -> ~5.4, energy drops by 0.2 eV/atom
honeycomb    ~25 K           z goes 3 -> ~5.2
square       never           shears at once, even at 10 K
===========  ==============  ===========================================

The square lattice has no barrier at all: its shear modulus is negative, so any
perturbation is downhill.  The others have a small barrier that thermal energy
crosses quickly.  In every case the potential energy *falls* on collapse, which
is the point -- close packing is simply lower.

A trap worth knowing about: on collapse, kagome turns into a triangular lattice,
and ``|psi_6|`` is 1 for *both*.  A run that watched only the order parameter
would report kagome as intact at 200 K, where it has in fact already collapsed.
The coordination number catches it (4 -> 5.4), which is why both are tracked.

That is not a limitation to work around; it is the physical content.  Real
honeycomb and kagome films (graphene, kagome metals, molecular networks) are
held open by *directional* bonding -- sp2 hybridisation, coordination chemistry,
a three-body term -- which an isotropic pair potential does not have.  Watching
them fall in is a direct demonstration of why those materials need more than a
pair potential to exist.

Each entry carries the symmetry metadata the rest of the code needs, because the
diagnostics are not transferable between lattices:

* the bond-orientational order parameter must match the local symmetry --
  ``|psi_6|`` is 1 on a triangular lattice but **0** on a square one, so using it
  everywhere would report a perfect square crystal as molten.  Honeycomb is a
  subtler case: its three bonds sit 120 degrees apart, which suggests ``n = 3``,
  but the two sublattices are rotated 60 degrees from each other and so carry
  ``psi_3`` of opposite phase.  The global average cancels exactly, and a perfect
  honeycomb crystal would read as zero.  ``n = 6`` is the sublattice-independent
  choice and is what the spec uses; per-site ``|psi_3|`` from
  :func:`~trilattice.observables.psi_n` is still 1 and remains useful locally;
* the first-shell cutoff must sit between the first and second neighbour shells,
  which are at different ratios for each structure;
* the ideal coordination differs, so "defect" means a different number;
* Delaunay coordination is the right topological tool only for the triangular
  lattice, where the triangulation is unique.  On a square lattice every
  plaquette is a degenerate quadrilateral that the triangulation splits
  arbitrarily, so a distance cutoff is used instead.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .lattice import SQRT3, Box, Configuration

__all__ = [
    "LatticeSpec",
    "LATTICES",
    "LATTICE_NAMES",
    "build_lattice",
    "lattice_spec",
    "cells_for",
    "triangular",
    "square",
    "honeycomb",
    "kagome",
    "emptiest_points",
]


# --------------------------------------------------------------------------- #
# Builders.  ``a`` is always the nearest-neighbour distance.
# --------------------------------------------------------------------------- #


def _pack(positions: np.ndarray, lx: float, ly: float, mass: float) -> Configuration:
    box = Box(lx, ly)
    pos = box.wrap(np.asarray(positions, dtype=np.float64))
    return Configuration(pos, box, np.full(pos.shape[0], float(mass)))


def _tile(basis: np.ndarray, nx: int, ny: int, cx: float, cy: float) -> np.ndarray:
    """Replicate a cell basis over an ``nx`` by ``ny`` grid of cells."""
    ix, iy = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    shifts = np.column_stack((ix.ravel() * cx, iy.ravel() * cy))
    return (shifts[:, None, :] + basis[None, :, :]).reshape(-1, 2)


def triangular(nx: int, ny: int, a: float, mass: float = 24.305) -> Configuration:
    """Close-packed lattice: 6 neighbours at ``a``.  2 atoms per cell."""
    basis = np.array([[0.0, 0.0], [0.5 * a, 0.5 * a * SQRT3]])
    return _pack(_tile(basis, nx, ny, a, a * SQRT3), nx * a, ny * a * SQRT3, mass)


def square(nx: int, ny: int, a: float, mass: float = 24.305) -> Configuration:
    """Simple square lattice: 4 neighbours at ``a``, 4 more at ``a*sqrt(2)``."""
    basis = np.array([[0.0, 0.0]])
    return _pack(_tile(basis, nx, ny, a, a), nx * a, ny * a, mass)


def honeycomb(nx: int, ny: int, a: float, mass: float = 24.305) -> Configuration:
    """Graphene topology: 3 neighbours at the bond length ``a``.

    Built on the standard rectangular cell of width ``sqrt(3) a`` and height
    ``3 a`` containing four sites, which tiles the plane without strain.
    """
    w, h = SQRT3 * a, 3.0 * a
    basis = np.array([
        [0.0, 0.0],
        [0.0, a],
        [0.5 * w, 1.5 * a],
        [0.5 * w, 2.5 * a],
    ])
    return _pack(_tile(basis, nx, ny, w, h), nx * w, ny * h, mass)


def kagome(nx: int, ny: int, a: float, mass: float = 24.305) -> Configuration:
    """Corner-sharing triangles: 4 neighbours at ``a``.

    A triangular Bravais lattice of constant ``2a`` with a three-site basis;
    the same centred-rectangular supercell as :func:`triangular` is used, so a
    cell holds six atoms.
    """
    lat = 2.0 * a
    basis = np.array([
        [0.0, 0.0], [a, 0.0], [0.5 * a, 0.5 * a * SQRT3],
        [0.5 * lat, 0.5 * lat * SQRT3],
        [0.5 * lat + a, 0.5 * lat * SQRT3],
        [0.5 * lat + 0.5 * a, 0.5 * lat * SQRT3 + 0.5 * a * SQRT3],
    ])
    return _pack(_tile(basis, nx, ny, lat, lat * SQRT3), nx * lat, ny * lat * SQRT3, mass)


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class LatticeSpec:
    """Everything the diagnostics need to know about a structure."""

    name: str
    builder: Callable[..., Configuration]
    coordination: int          #: nearest neighbours in the ideal structure
    psi_order: int             #: n of the n-fold bond-orientational order parameter
    angular_order: int         #: n of the three-body penalty [1-cos(n theta)]/2
    cutoff_factor: float       #: first-shell cutoff in units of a
    second_shell: float        #: second-neighbour distance in units of a
    atoms_per_cell: int
    cell_factors: tuple[float, float]   #: cell width and height in units of a
    #: |G1| * a -- the length of the shortest reciprocal-lattice vector, which
    #: sets how far an S(k) window has to reach to show the first Bragg ring
    bragg_factor: float
    density_factor: float      #: rho * a^2 for the ideal structure
    lj_ground_state: bool
    #: highest temperature (K) at which the structure was observed to survive
    #: 5 ps with the default parameters; 0 means it collapses at any temperature
    metastable_below: float
    note: str

    def cutoff(self, a: float) -> float:
        return self.cutoff_factor * a

    def density(self, a: float) -> float:
        return self.density_factor / a**2

    def first_bragg(self, a: float) -> float:
        """Length of the shortest reciprocal-lattice vector, in 1/A."""
        return self.bragg_factor / a

    def k_window(self, a: float, rings: float = 3.0) -> float:
        """Half-width of an S(k) window that shows ``rings`` Bragg shells."""
        return rings * self.first_bragg(a)

    def cells_for(self, n_target: int) -> tuple[int, int]:
        """Cell counts giving roughly ``n_target`` atoms in a roughly square box.

        Without this, switching lattice at fixed ``(nx, ny)`` changes the particle
        count by more than an order of magnitude -- a kagome cell holds six atoms
        against a square cell's one -- and the dashboard would stall the moment
        anyone picked kagome.
        """
        fx, fy = self.cell_factors
        cells = max(1.0, n_target / self.atoms_per_cell)
        nx = max(2, int(round(math.sqrt(cells * fy / fx))))
        ny = max(2, int(round(cells / nx)))
        return nx, ny

    def use_delaunay(self) -> bool:
        """Delaunay coordination is only unambiguous on the triangular lattice."""
        return self.name == "triangular"


LATTICES: dict[str, LatticeSpec] = {
    "triangular": LatticeSpec(
        "triangular", triangular, 6, 6, 6, 1.35, SQRT3, 2, (1.0, SQRT3), 4.0 * math.pi / SQRT3, 2.0 / SQRT3, True, float("inf"),
        "close packed; the ground state of any isotropic pair potential in 2-D",
    ),
    "square": LatticeSpec(
        "square", square, 4, 4, 4, 1.20, math.sqrt(2.0), 1, (1.0, 1.0), 2.0 * math.pi, 1.0, False, 0.0,
        "mechanically unstable under Lennard-Jones: it shears towards triangular "
        "immediately, at any temperature",
    ),
    "honeycomb": LatticeSpec(
        "honeycomb", honeycomb, 3, 6, 3, 1.30, SQRT3, 4, (SQRT3, 3.0), 4.0 * math.pi / 3.0, 4.0 / (3.0 * SQRT3), False, 25.0,
        "graphene topology; the most open of the four -- needs directional "
        "bonding to stay open, so it collapses under a pair potential",
    ),
    "kagome": LatticeSpec(
        "kagome", kagome, 4, 6, 6, 1.20, SQRT3, 6, (2.0, 2.0 * SQRT3), 2.0 * math.pi / SQRT3, 0.5 * SQRT3, False, 100.0,
        "corner-sharing triangles; four neighbours but only 75% of close-packed "
        "density, so it collapses readily",
    ),
}

LATTICE_NAMES = tuple(LATTICES)


def lattice_spec(kind: str) -> LatticeSpec:
    try:
        return LATTICES[kind]
    except KeyError:
        raise ValueError(
            f"unknown lattice {kind!r}; choose from {list(LATTICES)}"
        ) from None


def cells_for(kind: str, n_target: int) -> tuple[int, int]:
    """Convenience wrapper around :meth:`LatticeSpec.cells_for`."""
    return lattice_spec(kind).cells_for(n_target)


def build_lattice(kind: str, nx: int, ny: int, a: float,
                  mass: float = 24.305) -> Configuration:
    """Build any registered lattice.  ``a`` is the nearest-neighbour distance."""
    return lattice_spec(kind).builder(nx, ny, a, mass)


# --------------------------------------------------------------------------- #
# Lattice-agnostic interstitial placement
# --------------------------------------------------------------------------- #


def emptiest_points(positions: np.ndarray, box: Box, n: int = 1,
                    grid: int = 160) -> np.ndarray:
    """The ``n`` points of the cell furthest from any particle.

    A coarse grid search rather than a lattice-specific construction, so it works
    on a square lattice, a honeycomb, a half-melted configuration or a glass.
    Picked points are excluded from each other's neighbourhood so that several
    interstitials do not all land in the same hole.
    """
    gx = (np.arange(grid) + 0.5) * box.lx / grid
    gy = (np.arange(grid) + 0.5) * box.ly / grid
    xx, yy = np.meshgrid(gx, gy, indexing="ij")
    pts = np.column_stack((xx.ravel(), yy.ravel()))

    d = box.minimum_image(pts[:, None, :] - positions[None, :, :])
    dmin = np.min(np.hypot(d[:, :, 0], d[:, :, 1]), axis=1)

    chosen: list[np.ndarray] = []
    for _ in range(n):
        k = int(np.argmax(dmin))
        p = pts[k]
        chosen.append(p)
        dp = box.minimum_image(pts - p)
        dmin = np.minimum(dmin, np.hypot(dp[:, 0], dp[:, 1]))
    return np.asarray(chosen)
