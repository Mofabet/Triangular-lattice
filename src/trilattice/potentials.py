"""The Lennard-Jones 12-6 pair potential, truncated three different ways.

.. math::

    u(r) = 4\\varepsilon\\left[(\\sigma/r)^{12} - (\\sigma/r)^{6}\\right],
    \\qquad
    f(r) = -\\frac{\\mathrm{d}u}{\\mathrm{d}r}
         = \\frac{24\\varepsilon}{r}\\left[2(\\sigma/r)^{12} - (\\sigma/r)^{6}\\right]

Two details decide whether an MD code conserves energy or not.

**The factor of 2.**  The repulsive term differentiates to ``-12 sigma^12/r^13``
and the attractive one to ``+6 sigma^6/r^7``; the ratio is 2, not 1.  Dropping it
(as the original project did) produces a potential whose minimum sits at
``r = sigma`` instead of ``2^{1/6} sigma`` and whose force does not integrate to
the energy, so no integrator can conserve anything.

**The sign.**  ``f(r)`` is negative for ``r > 2^{1/6} sigma`` -- that is the
cohesion that holds the crystal together.  Taking ``abs()`` of the force (again,
as the original did) makes every interaction repulsive and the lattice sublimates
on the first step.

Truncation
----------
``cut``
    Bare truncation.  ``u`` jumps by ``u(r_c)`` at the cutoff, so every crossing
    injects energy; usable only with a thermostat, and then only with the
    analytic tail corrections in :func:`tail_energy` / :func:`tail_pressure`.
``shifted``
    ``u(r) - u(r_c)``.  Energy continuous, force still discontinuous.  This is
    the textbook default and what most published LJ data uses.
``shifted-force``
    Stoddard-Ford: ``u(r) - u(r_c) - (r - r_c) u'(r_c)``.  Both ``u`` and ``f``
    vanish smoothly at ``r_c``; energy conservation in NVE improves by two to
    three orders of magnitude.  This is the default here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal, Sequence

import numpy as np

TruncationMode = Literal["cut", "shifted", "shifted-force"]
_MODES = ("cut", "shifted", "shifted-force")

#: Position of the minimum of the bare LJ potential, in units of sigma.
R_MIN_OVER_SIGMA = 2.0 ** (1.0 / 6.0)


@dataclass
class LennardJones:
    """Lennard-Jones pair style, optionally multi-component.

    Parameters
    ----------
    epsilon, sigma
        Scalars for a one-component system, or sequences of length ``n_types``
        for a mixture (combined with Lorentz-Berthelot rules).
    cutoff
        Global cutoff radius in angstrom.
    mode
        One of ``"cut"``, ``"shifted"``, ``"shifted-force"``.
    """

    epsilon: float | Sequence[float]
    sigma: float | Sequence[float]
    cutoff: float
    mode: TruncationMode = "shifted-force"

    # derived tables, filled in __post_init__
    eps_table: np.ndarray = field(init=False, repr=False)
    sig_table: np.ndarray = field(init=False, repr=False)
    sig2_table: np.ndarray = field(init=False, repr=False)
    energy_shift: np.ndarray = field(init=False, repr=False)
    force_at_cutoff: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.mode not in _MODES:
            raise ValueError(f"mode must be one of {_MODES}, got {self.mode!r}")
        if self.cutoff <= 0.0:
            raise ValueError("cutoff must be positive")

        eps = np.atleast_1d(np.asarray(self.epsilon, dtype=np.float64))
        sig = np.atleast_1d(np.asarray(self.sigma, dtype=np.float64))
        if eps.size != sig.size:
            raise ValueError("epsilon and sigma must have the same number of types")
        if np.any(eps <= 0) or np.any(sig <= 0):
            raise ValueError("epsilon and sigma must be positive")

        # Lorentz-Berthelot mixing
        self.eps_table = np.sqrt(np.outer(eps, eps))
        self.sig_table = 0.5 * (sig[:, None] + sig[None, :])
        self.sig2_table = self.sig_table**2

        rc = self.cutoff
        src2 = self.sig2_table / rc**2
        src6 = src2**3
        src12 = src6**2
        u_rc = 4.0 * self.eps_table * (src12 - src6)
        f_rc = 24.0 * self.eps_table * (2.0 * src12 - src6) / rc

        zeros = np.zeros_like(u_rc)
        if self.mode == "cut":
            self.energy_shift = zeros
            self.force_at_cutoff = zeros
        elif self.mode == "shifted":
            self.energy_shift = -u_rc
            self.force_at_cutoff = zeros
        else:  # shifted-force
            self.energy_shift = -u_rc - rc * f_rc
            self.force_at_cutoff = f_rc

    # ------------------------------------------------------------------ #
    @property
    def n_types(self) -> int:
        return self.eps_table.shape[0]

    @property
    def r_min(self) -> np.ndarray:
        """Position of the potential minimum for each type pair."""
        return R_MIN_OVER_SIGMA * self.sig_table

    def well_depth(self) -> np.ndarray:
        return self.eps_table.copy()

    # ------------------------------------------------------------------ #
    # Scalar / array evaluation (for plots, tests and lattice sums)
    # ------------------------------------------------------------------ #
    def energy(self, r: np.ndarray | float, ti: int = 0, tj: int = 0) -> np.ndarray:
        """Pair energy, already truncated according to ``mode``."""
        r = np.asarray(r, dtype=np.float64)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            sr2 = self.sig2_table[ti, tj] / r**2
            sr6 = sr2**3
            u = 4.0 * self.eps_table[ti, tj] * (sr6 * sr6 - sr6)
            u = u + self.energy_shift[ti, tj] + r * self.force_at_cutoff[ti, tj]
        return np.where(r < self.cutoff, u, 0.0)

    def force(self, r: np.ndarray | float, ti: int = 0, tj: int = 0) -> np.ndarray:
        """Scalar pair force ``-du/dr`` (positive = repulsive)."""
        r = np.asarray(r, dtype=np.float64)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            sr2 = self.sig2_table[ti, tj] / r**2
            sr6 = sr2**3
            f = 24.0 * self.eps_table[ti, tj] * (2.0 * sr6 * sr6 - sr6) / r
            f = f - self.force_at_cutoff[ti, tj]
        return np.where(r < self.cutoff, f, 0.0)

    def bare_energy(self, r: np.ndarray | float, ti: int = 0, tj: int = 0) -> np.ndarray:
        """Untruncated ``4 eps [(sig/r)^12 - (sig/r)^6]``."""
        r = np.asarray(r, dtype=np.float64)
        sr6 = (self.sig2_table[ti, tj] / r**2) ** 3
        return 4.0 * self.eps_table[ti, tj] * (sr6 * sr6 - sr6)

    # ------------------------------------------------------------------ #
    # 2-D long-range corrections (only meaningful for mode="cut")
    # ------------------------------------------------------------------ #
    def tail_energy(self, density: float, n_particles: int) -> float:
        r"""Mean-field estimate of the energy beyond the cutoff.

        In :math:`d = 2`, assuming :math:`g(r) \simeq 1` past the cutoff,

        .. math::
            U_\mathrm{tail} = \tfrac12 N \rho \int_{r_c}^{\infty} u(r)\, 2\pi r\,
            \mathrm{d}r = 4\pi N \rho \varepsilon
            \left[\frac{\sigma^{12}}{10 r_c^{10}} - \frac{\sigma^{6}}{4 r_c^{4}}\right].
        """
        if self.mode != "cut":
            return 0.0
        eps = float(self.eps_table[0, 0])
        sig = float(self.sig_table[0, 0])
        rc = self.cutoff
        return (
            4.0
            * math.pi
            * n_particles
            * density
            * eps
            * (sig**12 / (10.0 * rc**10) - sig**6 / (4.0 * rc**4))
        )

    def tail_pressure(self, density: float) -> float:
        r"""Mean-field virial correction to the 2-D pressure [eV/A^2]."""
        if self.mode != "cut":
            return 0.0
        eps = float(self.eps_table[0, 0])
        sig = float(self.sig_table[0, 0])
        rc = self.cutoff
        return (
            12.0
            * math.pi
            * density**2
            * eps
            * (sig**12 / (5.0 * rc**10) - sig**6 / (4.0 * rc**4))
        )

    # ------------------------------------------------------------------ #
    def lattice_energy(self, a: float, n_shells: int = 40) -> float:
        """Static energy per atom of a perfect triangular lattice, direct sum.

        Sums the *untruncated* potential over all sites of the infinite lattice
        out to ``n_shells * a``; used to locate the T = 0 lattice constant and to
        check that the simulation's potential energy at low temperature is right.
        """
        k = n_shells
        m, n = np.meshgrid(np.arange(-k, k + 1), np.arange(-k, k + 1), indexing="ij")
        norm2 = m.astype(np.float64) ** 2 + m * n + n.astype(np.float64) ** 2
        norm2 = norm2[norm2 > 0]
        r = a * np.sqrt(norm2)
        r = r[r <= k * a]
        return 0.5 * float(np.sum(self.bare_energy(r)))

    def summary(self) -> str:  # pragma: no cover - cosmetic
        e, s = self.eps_table[0, 0], self.sig_table[0, 0]
        return (
            f"LJ(eps={e:.4f} eV, sigma={s:.4f} A, rc={self.cutoff:.3f} A"
            f" = {self.cutoff / s:.2f} sigma, mode={self.mode}, ntypes={self.n_types})"
        )


def kob_andersen_like(epsilon: float, sigma: float, cutoff: float, **kw) -> LennardJones:
    """A size-asymmetric two-component table that frustrates crystallisation.

    Uses ``sigma_BB = 0.88 sigma``, ``eps_BB = 0.5 eps`` -- the same asymmetry as
    the Kob-Andersen mixture, which in 2-D is a well-studied glass former.
    """
    return LennardJones(
        epsilon=[epsilon, 0.5 * epsilon],
        sigma=[sigma, 0.88 * sigma],
        cutoff=cutoff,
        **kw,
    )


__all__ = ["LennardJones", "TruncationMode", "R_MIN_OVER_SIGMA", "kob_andersen_like"]
