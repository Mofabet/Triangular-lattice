"""Geometry: the periodic box and the triangular lattice that fills it.

The two-dimensional triangular (hexagonal) lattice is the ground state of almost
every isotropic pair potential in 2-D, which is why it is the natural starting
configuration here.  Its primitive vectors are

    a1 = a (1, 0),      a2 = a (1/2, sqrt(3)/2)

which do **not** span a rectangle, so they cannot be used directly with the
orthorhombic minimum-image convention.  Instead we use the *centred rectangular*
conventional cell

    A1 = a (1, 0),      A2 = a (0, sqrt(3)),     basis = {(0, 0), (a/2, a*sqrt(3)/2)}

which tiles the plane with exactly the same point set, contains two atoms, and is
rectangular.  A supercell of ``nx x ny`` such cells holds ``N = 2 nx ny`` atoms in
a box of ``Lx = nx a`` by ``Ly = ny a sqrt(3)``.

Commensurability matters: if the box edges are not integer multiples of the
lattice periods the crystal is strained by the periodic images and melts at the
wrong temperature.  :func:`triangular_lattice` guarantees commensurability by
construction; that is the entire reason for preferring it over "place atoms then
pick a box".
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

SQRT3 = math.sqrt(3.0)


# --------------------------------------------------------------------------- #
# Periodic box
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Box:
    """A 2-D orthorhombic periodic cell of size ``lx x ly``."""

    lx: float
    ly: float

    def __post_init__(self) -> None:
        if not (self.lx > 0.0 and self.ly > 0.0):
            raise ValueError(f"box edges must be positive, got ({self.lx}, {self.ly})")

    # -- basic properties --------------------------------------------------- #
    @property
    def lengths(self) -> np.ndarray:
        return np.array([self.lx, self.ly], dtype=np.float64)

    @property
    def area(self) -> float:
        return self.lx * self.ly

    @property
    def max_cutoff(self) -> float:
        """Largest cutoff for which the minimum-image convention is unambiguous."""
        return 0.5 * min(self.lx, self.ly)

    # -- the two operations that every PBC code needs ----------------------- #
    def wrap(self, positions: np.ndarray) -> np.ndarray:
        """Fold positions back into ``[0, L)`` (returns a new array)."""
        return np.mod(positions, self.lengths)

    def wrap_inplace(self, positions: np.ndarray) -> None:
        np.mod(positions, self.lengths, out=positions)

    def minimum_image(self, delta: np.ndarray) -> np.ndarray:
        """Map separation vectors onto the nearest periodic image.

        Uses rounding rather than ``if`` chains so that it stays correct for
        displacements larger than one box length (which do occur in a Verlet-list
        rebuild after a violent step, and which the ``if dx > L: dx -= L`` idiom
        silently gets wrong).
        """
        lengths = self.lengths
        return delta - lengths * np.round(delta / lengths)

    def scaled(self, factor: float) -> "Box":
        return Box(self.lx * factor, self.ly * factor)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"Box(lx={self.lx:.4f}, ly={self.ly:.4f}, area={self.area:.2f})"


# --------------------------------------------------------------------------- #
# Configuration container
# --------------------------------------------------------------------------- #


@dataclass
class Configuration:
    """Positions + per-particle metadata + the box they live in."""

    positions: np.ndarray  #: (N, 2) float64, angstrom
    box: Box
    masses: np.ndarray  #: (N,) float64, u
    types: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int32))
    #: index of the perfect-lattice site each particle started from, or -1
    site_index: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int32))

    def __post_init__(self) -> None:
        self.positions = np.ascontiguousarray(self.positions, dtype=np.float64)
        if self.positions.ndim != 2 or self.positions.shape[1] != 2:
            raise ValueError("positions must have shape (N, 2)")
        n = self.positions.shape[0]
        self.masses = np.ascontiguousarray(np.broadcast_to(self.masses, (n,)), dtype=np.float64)
        if self.types.size == 0:
            self.types = np.zeros(n, dtype=np.int32)
        if self.site_index.size == 0:
            self.site_index = np.arange(n, dtype=np.int32)

    @property
    def n_particles(self) -> int:
        return self.positions.shape[0]

    @property
    def n_types(self) -> int:
        return int(self.types.max()) + 1 if self.types.size else 0

    @property
    def density(self) -> float:
        """Number density in A^-2."""
        return self.n_particles / self.box.area

    def copy(self) -> "Configuration":
        return Configuration(
            self.positions.copy(),
            self.box,
            self.masses.copy(),
            self.types.copy(),
            self.site_index.copy(),
        )


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #


def triangular_lattice(
    nx: int,
    ny: int,
    a: float,
    mass: float = 24.305,
    *,
    strain: tuple[float, float] = (0.0, 0.0),
) -> Configuration:
    """Build a commensurate ``2*nx*ny``-atom triangular lattice.

    Parameters
    ----------
    nx, ny
        Number of centred-rectangular cells along x and y.
    a
        Nearest-neighbour distance (the lattice constant) in angstrom.
    mass
        Particle mass in u (24.305 = magnesium, the original project's choice).
    strain
        Optional ``(exx, eyy)`` engineering strain applied to both the box and
        the atoms, for computing elastic constants by finite differences.

    Notes
    -----
    Sublattice A sits at ``(i a, j a sqrt(3))`` and sublattice B at
    ``(i a + a/2, j a sqrt(3) + a sqrt(3)/2)``.  Each atom then has exactly six
    neighbours at distance ``a`` -- the defining property of the triangular
    lattice -- and the construction is periodic in both directions by
    construction.
    """
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be >= 1")
    if a <= 0.0:
        raise ValueError("lattice constant must be positive")

    ix, iy = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    ix = ix.ravel().astype(np.float64)
    iy = iy.ravel().astype(np.float64)

    sub_a = np.column_stack((ix * a, iy * a * SQRT3))
    sub_b = np.column_stack((ix * a + 0.5 * a, iy * a * SQRT3 + 0.5 * a * SQRT3))

    positions = np.empty((2 * nx * ny, 2), dtype=np.float64)
    positions[0::2] = sub_a
    positions[1::2] = sub_b

    lx, ly = nx * a, ny * a * SQRT3
    exx, eyy = strain
    if exx or eyy:
        positions[:, 0] *= 1.0 + exx
        positions[:, 1] *= 1.0 + eyy
        lx *= 1.0 + exx
        ly *= 1.0 + eyy

    box = Box(lx, ly)
    return Configuration(positions, box, np.full(positions.shape[0], float(mass)))


def interstitial_sites(config: Configuration, a: float) -> np.ndarray:
    """The two honeycomb hollow sites per lattice cell.

    A triangular lattice has two inequivalent three-fold hollows per atom; these
    are the physically meaningful interstitial positions.  Dropping an extra atom
    at a *uniformly random* point -- as the original code did -- puts it on top of
    a host atom with probability ~sigma^2 rho, which at this density is ~10% per
    insertion and makes the integrator explode on step one.
    """
    base = config.positions
    offsets = np.array(
        [[0.5 * a, a * SQRT3 / 6.0], [0.0, a * SQRT3 / 3.0]],
        dtype=np.float64,
    )
    sites = (base[:, None, :] + offsets[None, :, :]).reshape(-1, 2)
    return config.box.wrap(sites)


def add_defects(
    config: Configuration,
    *,
    n_vacancies: int = 0,
    n_interstitials: int = 0,
    a: float | None = None,
    interstitial_mass: float | None = None,
    interstitial_type: int = 0,
    rng: np.random.Generator | None = None,
    min_separation: float = 0.0,
) -> Configuration:
    """Remove ``n_vacancies`` atoms and insert ``n_interstitials`` new ones.

    Interstitials are placed on true hollow sites (see :func:`interstitial_sites`)
    and rejected if they land closer than ``min_separation`` to an existing atom.

    Note the geometric ceiling: a three-fold hollow of a triangular lattice sits
    at ``a/sqrt(3) ~ 0.577 a`` from each of its three neighbours, so any
    ``min_separation`` above that rejects every candidate.  ``0.5 a`` is a safe
    default.
    """
    rng = np.random.default_rng() if rng is None else rng
    cfg = config.copy()

    if n_vacancies:
        if n_vacancies >= cfg.n_particles:
            raise ValueError("cannot remove every particle")
        keep = rng.permutation(cfg.n_particles)[n_vacancies:]
        keep.sort()
        cfg = Configuration(
            cfg.positions[keep],
            cfg.box,
            cfg.masses[keep],
            cfg.types[keep],
            cfg.site_index[keep],
        )

    if n_interstitials:
        if a is None:
            raise ValueError("lattice constant `a` is required to place interstitials")
        candidates = interstitial_sites(cfg, a)
        rng.shuffle(candidates)
        accepted: list[np.ndarray] = []
        occupied = cfg.positions
        for cand in candidates:
            if len(accepted) >= n_interstitials:
                break
            reference = np.vstack([occupied, *accepted]) if accepted else occupied
            d = cfg.box.minimum_image(reference - cand)
            if np.min(np.hypot(d[:, 0], d[:, 1])) >= min_separation:
                accepted.append(cand)
        if len(accepted) < n_interstitials:
            raise RuntimeError(
                f"only placed {len(accepted)}/{n_interstitials} interstitials with "
                f"min_separation={min_separation:.3f} A"
            )
        extra = np.asarray(accepted, dtype=np.float64)
        mass = cfg.masses[0] if interstitial_mass is None else float(interstitial_mass)
        cfg = Configuration(
            np.vstack([cfg.positions, extra]),
            cfg.box,
            np.concatenate([cfg.masses, np.full(len(extra), mass)]),
            np.concatenate([cfg.types, np.full(len(extra), interstitial_type, dtype=np.int32)]),
            np.concatenate([cfg.site_index, np.full(len(extra), -1, dtype=np.int32)]),
        )

    return cfg


def make_binary_mixture(
    config: Configuration,
    fraction_b: float,
    *,
    mass_b: float | None = None,
    rng: np.random.Generator | None = None,
) -> Configuration:
    """Randomly relabel a fraction of the atoms as type 1.

    Combined with a size-asymmetric pair table this turns the crystal into a
    Kob-Andersen-like glass former: the lattice can no longer accommodate both
    species, so on cooling the system arrests in an amorphous state instead of
    recrystallising.  It is the cheapest interesting extension of the original
    single-component model.
    """
    rng = np.random.default_rng() if rng is None else rng
    cfg = config.copy()
    n_b = int(round(fraction_b * cfg.n_particles))
    idx = rng.permutation(cfg.n_particles)[:n_b]
    cfg.types = cfg.types.copy()
    cfg.types[idx] = 1
    if mass_b is not None:
        cfg.masses = cfg.masses.copy()
        cfg.masses[idx] = mass_b
    return cfg


def replicate(config: Configuration, nx: int, ny: int) -> Configuration:
    """Tile a configuration ``nx`` by ``ny`` times, enlarging the box to match.

    Needed when a cell is too small for its own cutoff: the minimum-image
    convention requires ``r_cut < L/2`` in every direction, and a smaller cell
    makes an atom interact with a neighbour *and* with that neighbour's own
    periodic image at the same time -- silently double counting.
    """
    if nx < 1 or ny < 1:
        raise ValueError("replication factors must be >= 1")
    shifts = np.array(
        [[i * config.box.lx, j * config.box.ly] for i in range(nx) for j in range(ny)],
        dtype=np.float64,
    )
    positions = (config.positions[None, :, :] + shifts[:, None, :]).reshape(-1, 2)
    reps = nx * ny
    return Configuration(
        positions,
        Box(config.box.lx * nx, config.box.ly * ny),
        np.tile(config.masses, reps),
        np.tile(config.types, reps),
        np.tile(config.site_index, reps),
    )


def required_replication(box: Box, r_list: float) -> tuple[int, int]:
    """Smallest ``(nx, ny)`` tiling for which ``r_list < L/2`` holds."""
    return (
        max(1, int(math.ceil(2.0 * r_list / box.lx))),
        max(1, int(math.ceil(2.0 * r_list / box.ly))),
    )


def perfect_lattice_neighbour_shells(a: float, r_max: float) -> tuple[np.ndarray, np.ndarray]:
    """Radii and multiplicities of the coordination shells of the ideal lattice.

    Returned as ``(radii, counts)`` with ``radii <= r_max``.  Shell radii are
    ``a*sqrt(m^2 + m n + n^2)`` for integer ``m, n`` -- the norm form of the
    triangular lattice.  Used to annotate g(r) plots and to check lattice sums.
    """
    k = int(math.ceil(r_max / a)) + 2
    m, n = np.meshgrid(np.arange(-k, k + 1), np.arange(-k, k + 1), indexing="ij")
    norm2 = (m.astype(np.float64) ** 2 + m * n + n.astype(np.float64) ** 2)
    norm2 = norm2[norm2 > 0]
    r = a * np.sqrt(norm2)
    r = r[r <= r_max]
    radii, counts = np.unique(np.round(r, 9), return_counts=True)
    return radii, counts


__all__ = [
    "SQRT3",
    "Box",
    "Configuration",
    "triangular_lattice",
    "interstitial_sites",
    "add_defects",
    "replicate",
    "required_replication",
    "make_binary_mixture",
    "perfect_lattice_neighbour_shells",
]
