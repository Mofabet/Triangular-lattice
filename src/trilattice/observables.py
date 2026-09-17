"""Structural and dynamical observables for a 2-D system.

Dimensionality is not cosmetic here.  The radial distribution normalisation is
``2 pi r dr`` and not ``4 pi r^2 dr``; the diffusion constant is ``MSD/(4t)`` and
not ``MSD/(6t)``; the virial pressure divides by ``2A`` and not ``3V``; and the
ordinary Lindemann criterion does not exist at all, because in two dimensions
Mermin-Wagner fluctuations make ``<u^2>`` diverge logarithmically with system
size.  Each of those is a place where a 3-D formula silently produces a
plausible-looking wrong number.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .lattice import Box
from .units import KB, MVV2E


# --------------------------------------------------------------------------- #
# Thermodynamics
# --------------------------------------------------------------------------- #


def pressure_2d(temperature: float, virial: float, n_particles: int, area: float) -> float:
    r"""2-D virial pressure in eV/A^2.

    .. math::
        P = \frac{1}{d A}\left(d N k_B T + \sum_{i<j}\mathbf{r}_{ij}\cdot
        \mathbf{f}_{ij}\right),\qquad d = 2
    """
    return (2.0 * n_particles * KB * temperature + virial) / (2.0 * area)


def heat_capacity_nvt(energies: np.ndarray, temperature: float, n_particles: int) -> float:
    r"""Heat capacity per particle in units of k_B, from energy fluctuations.

    .. math::
        \frac{C_v}{N k_B} = \frac{\langle E^2\rangle - \langle E\rangle^2}
                                  {N k_B^2 T^2}

    Note the *square* of the Boltzmann constant: the variance of an energy has
    units of energy squared.  A harmonic 2-D crystal gives 2 k_B per atom
    (k_B kinetic + k_B potential), which is the Dulong-Petit value in two
    dimensions and a useful check on the low-temperature end of any scan.
    """
    var = float(np.var(np.asarray(energies, dtype=np.float64)))
    return var / (KB**2 * temperature**2 * n_particles)


# --------------------------------------------------------------------------- #
# Pair structure
# --------------------------------------------------------------------------- #


@dataclass
class RDFAccumulator:
    """Running histogram for the radial distribution function g(r)."""

    box: Box
    r_max: float
    n_bins: int = 300

    def __post_init__(self) -> None:
        self.r_max = min(self.r_max, self.box.max_cutoff)
        self.edges = np.linspace(0.0, self.r_max, self.n_bins + 1)
        self.centers = 0.5 * (self.edges[1:] + self.edges[:-1])
        self.hist = np.zeros(self.n_bins, dtype=np.float64)
        self.n_frames = 0
        self.n_particles = 0

    def accumulate(self, positions: np.ndarray) -> None:
        n = positions.shape[0]
        i, j = np.triu_indices(n, k=1)
        d = self.box.minimum_image(positions[i] - positions[j])
        r = np.hypot(d[:, 0], d[:, 1])
        self.hist += np.histogram(r, bins=self.edges)[0]
        self.n_frames += 1
        self.n_particles = n

    def result(self) -> tuple[np.ndarray, np.ndarray]:
        if self.n_frames == 0:
            return self.centers, np.zeros_like(self.centers)
        rho = self.n_particles / self.box.area
        shell_area = np.pi * (self.edges[1:] ** 2 - self.edges[:-1] ** 2)
        ideal = 0.5 * self.n_particles * rho * shell_area  # pairs counted once
        return self.centers, self.hist / (self.n_frames * ideal)


def coordination_number(g_r: np.ndarray, r: np.ndarray, density: float, r_cut: float) -> float:
    """First-shell coordination ``2 pi rho int_0^{r_cut} g(r) r dr``."""
    m = r <= r_cut
    return float(2.0 * np.pi * density * np.trapezoid(g_r[m] * r[m], r[m]))


def structure_factor(
    positions: np.ndarray, box: Box, n_max: int = 24
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """S(k) on the discrete reciprocal lattice of the periodic box.

    Returns ``(kx, ky, S)`` with ``S`` of shape ``(2*n_max+1, 2*n_max+1)``.
    Only wavevectors commensurate with the box are physically meaningful; using
    a continuous k-grid instead produces spurious ringing.
    """
    n = positions.shape[0]
    nx = np.arange(-n_max, n_max + 1)
    ny = np.arange(-n_max, n_max + 1)
    kx = 2.0 * np.pi * nx / box.lx
    ky = 2.0 * np.pi * ny / box.ly
    phase_x = np.exp(-1j * np.outer(kx, positions[:, 0]))  # (nkx, N)
    phase_y = np.exp(-1j * np.outer(ky, positions[:, 1]))  # (nky, N)
    rho_k = np.einsum("an,bn->ab", phase_x, phase_y)
    s = (np.abs(rho_k) ** 2) / n
    return kx, ky, s


# --------------------------------------------------------------------------- #
# Bond-orientational order
# --------------------------------------------------------------------------- #


def neighbour_pairs_within(box: Box, positions: np.ndarray, cutoff: float):
    n = positions.shape[0]
    i, j = np.triu_indices(n, k=1)
    d = box.minimum_image(positions[i] - positions[j])
    r2 = d[:, 0] ** 2 + d[:, 1] ** 2
    m = r2 < cutoff**2
    return i[m], j[m], d[m]


def psi6(positions: np.ndarray, box: Box, cutoff: float) -> np.ndarray:
    r"""Per-particle hexatic order parameter.

    .. math::
        \psi_6(j) = \frac{1}{n_j}\sum_{k \in \mathrm{nn}(j)} e^{6 i \theta_{jk}}

    ``|psi_6| = 1`` for a perfect triangular environment and averages to ~0 in an
    isotropic liquid.  ``cutoff`` should sit in the first minimum of g(r); the
    canonical choice for a triangular lattice of spacing ``a`` is ``1.35 a``.
    """
    n = positions.shape[0]
    i, j, d = neighbour_pairs_within(box, positions, cutoff)
    theta = np.arctan2(d[:, 1], d[:, 0])
    phase = np.exp(6j * theta)
    acc = np.zeros(n, dtype=np.complex128)
    cnt = np.zeros(n, dtype=np.float64)
    # theta_ji = theta_ij + pi, and exp(6 i pi) = 1, so both ends take the
    # same phase -- the six-fold symmetry is exactly what makes psi_6 bond-
    # direction agnostic.
    np.add.at(acc, i, phase)
    np.add.at(acc, j, phase)
    np.add.at(cnt, i, 1.0)
    np.add.at(cnt, j, 1.0)
    cnt[cnt == 0] = 1.0
    return acc / cnt


def global_psi6(positions: np.ndarray, box: Box, cutoff: float) -> float:
    """``|<psi_6>|`` averaged over particles: the hexatic order parameter."""
    return float(np.abs(np.mean(psi6(positions, box, cutoff))))


def psi6_correlation(
    positions: np.ndarray, box: Box, cutoff: float, r_max: float, n_bins: int = 60
):
    r"""Orientational correlation ``g_6(r) = <psi_6^*(0) psi_6(r)> / g(r)``.

    The decay law distinguishes the three KTHNY phases: ``g_6 -> const`` in the
    solid, ``g_6 ~ r^{-eta_6}`` with ``eta_6 <= 1/4`` in the hexatic, and
    exponential in the liquid.
    """
    p = psi6(positions, box, cutoff)
    n = positions.shape[0]
    i, j = np.triu_indices(n, k=1)
    d = box.minimum_image(positions[i] - positions[j])
    r = np.hypot(d[:, 0], d[:, 1])
    m = r < r_max
    r, i, j = r[m], i[m], j[m]
    w = np.real(np.conj(p[i]) * p[j])
    edges = np.linspace(0.0, r_max, n_bins + 1)
    num = np.histogram(r, bins=edges, weights=w)[0]
    den = np.histogram(r, bins=edges)[0]
    centers = 0.5 * (edges[1:] + edges[:-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        g6 = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    return centers, g6


# --------------------------------------------------------------------------- #
# Topological defects
# --------------------------------------------------------------------------- #


def coordination_by_delaunay(positions: np.ndarray, box: Box) -> np.ndarray:
    """Number of Delaunay (= Voronoi-face) neighbours of every particle.

    A 3x3 tiling is triangulated so that the result is properly periodic.  In the
    2-D melting picture this is *the* defect diagnostic: the crystal is a lattice
    of six-fold sites, dislocations are bound 5-7 pairs, and the hexatic-liquid
    transition is the unbinding of those pairs into free disclinations.
    """
    from scipy.spatial import Delaunay

    n = positions.shape[0]
    shifts = np.array(
        [[dx * box.lx, dy * box.ly] for dx in (-1, 0, 1) for dy in (-1, 0, 1)],
        dtype=np.float64,
    )
    tiled = (positions[None, :, :] + shifts[:, None, :]).reshape(-1, 2)
    tri = Delaunay(tiled)
    indptr, indices = tri.vertex_neighbor_vertices

    coord = np.zeros(n, dtype=np.int32)
    center_offset = 4 * n  # the (0, 0) shift is the 5th of the nine
    for k in range(n):
        v = center_offset + k
        coord[k] = indptr[v + 1] - indptr[v]
    return coord


def defect_fraction(positions: np.ndarray, box: Box) -> float:
    """Fraction of particles whose coordination differs from six."""
    c = coordination_by_delaunay(positions, box)
    return float(np.mean(c != 6))


# --------------------------------------------------------------------------- #
# Dynamics
# --------------------------------------------------------------------------- #


def mean_squared_displacement(unwrapped: np.ndarray, origin: np.ndarray | None = None) -> np.ndarray:
    """MSD(t) from a trajectory of unwrapped coordinates, shape ``(T, N, 2)``."""
    traj = np.asarray(unwrapped, dtype=np.float64)
    ref = traj[0] if origin is None else origin
    d = traj - ref[None, :, :]
    return np.mean(np.sum(d**2, axis=2), axis=1)


def diffusion_coefficient(times: np.ndarray, msd: np.ndarray, fit_from: float = 0.5) -> float:
    r"""Einstein relation in 2-D: ``D = MSD/(4t)``, fitted on the late-time half."""
    times = np.asarray(times)
    msd = np.asarray(msd)
    m = times >= fit_from * times[-1]
    if m.sum() < 3:
        return float("nan")
    slope = np.polyfit(times[m], msd[m], 1)[0]
    return float(slope / 4.0)


def velocity_autocorrelation(velocities: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Normalised VACF from a ``(T, N, 2)`` velocity trajectory."""
    v = np.asarray(velocities, dtype=np.float64)
    t = v.shape[0]
    max_lag = t // 2 if max_lag is None else min(max_lag, t - 1)
    c = np.empty(max_lag + 1)
    for lag in range(max_lag + 1):
        c[lag] = np.mean(np.sum(v[: t - lag] * v[lag:], axis=2))
    return c / c[0]


def lindemann_2d(unwrapped: np.ndarray, box: Box, a: float) -> float:
    r"""Modified (Bedanov-Gadiyak-Lozovik) Lindemann parameter.

    .. math::
        \gamma_m = \frac{\langle |\mathbf{u}_i - \mathbf{u}_j|^2 \rangle_{\langle ij\rangle}}{a^2}

    Built from *relative* displacements of neighbouring pairs, which cancels the
    long-wavelength phonons responsible for the Mermin-Wagner divergence of the
    ordinary Lindemann ratio.  Melting is signalled near ``gamma_m ~ 0.1``.
    """
    traj = np.asarray(unwrapped, dtype=np.float64)
    ref = traj[0]
    i, j, _ = neighbour_pairs_within(box, box.wrap(ref), 1.35 * a)
    if i.size == 0:
        return float("nan")
    u = traj - ref[None, :, :]
    du = u[:, i, :] - u[:, j, :]
    return float(np.mean(np.sum(du**2, axis=2)) / a**2)


__all__ = [
    "pressure_2d",
    "heat_capacity_nvt",
    "RDFAccumulator",
    "coordination_number",
    "structure_factor",
    "psi6",
    "global_psi6",
    "psi6_correlation",
    "coordination_by_delaunay",
    "defect_fraction",
    "mean_squared_displacement",
    "diffusion_coefficient",
    "velocity_autocorrelation",
    "lindemann_2d",
    "neighbour_pairs_within",
]
