"""A three-body angular term, so that open lattices stop collapsing.

An isotropic pair potential has exactly one ground state in two dimensions: the
triangular lattice.  Nothing about the starting configuration changes that,
because with a central force there is no energy cost attached to bond *angles*,
so the minimum is whatever packs the most neighbours at the potential minimum.
Square, honeycomb and kagome are therefore metastable at best under Lennard-Jones
alone (see :mod:`trilattice.lattices`).

Real open two-dimensional films -- graphene, silicene, kagome metals, molecular
networks -- are held open by *directional* bonding.  The cheapest way to put that
into a classical model is the Stillinger-Weber construction: add a term that
costs energy whenever a bond angle departs from the ones the target structure
wants.

.. math::

    U_3 = \\lambda \\sum_i \\sum_{j<k \\in \\mathrm{nn}(i)}
          f(r_{ij})\\, f(r_{ik})\\, P(\\cos\\theta_{jik})

The angular penalty used here is

.. math::

    P(\\theta) = \\tfrac12\\left[1 - \\cos(n\\theta)\\right]

which is zero exactly at every multiple of :math:`2\\pi/n` and rises to 1 halfway
between.  Choosing ``n`` to match the structure makes *all* of its bond angles
cost nothing:

===========  ===  =================================  ================
lattice      n    bond angles at a site              U_3 of the ideal
===========  ===  =================================  ================
honeycomb    3    120 deg                            exactly 0
square       4    90, 180 deg                        exactly 0
triangular   6    60, 120, 180 deg                   exactly 0
kagome       6    60, 120, 180 deg                   exactly 0
===========  ===  =================================  ================

Two consequences worth being explicit about.

**The ideal lattice keeps its Lennard-Jones energy.**  Because the penalty
vanishes at every angle the structure contains, ``U_3 = 0`` for the perfect
crystal.  The term adds no binding; it only raises the energy of *departures*
from the wanted geometry.  So the T = 0 lattice constant and cohesive energy are
unchanged, and the term can be switched on and off without re-fitting anything.

**This is a per-structure parametrisation, not a universal potential.**  ``n = 3``
makes honeycomb stable and simultaneously penalises the triangular lattice
heavily -- its 60 degree angles sit at the maximum of the ``n = 3`` penalty.
That is not a defect of the implementation; it is what a Stillinger-Weber
potential *is*.  SW for silicon is fitted to the tetrahedral angle and does not
describe close-packed silicon either.  Pick ``n`` for the structure you mean to
study, and say so when you report the result.

Note also that kagome shares ``n = 6`` with triangular, since both have bond
angles at multiples of 60 degrees.  The angular term alone therefore does not
distinguish them; what keeps kagome open is its coordination of four, held by
the radial cutoff of the three-body term.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .lattice import Box

try:  # pragma: no cover - environment dependent
    from numba import njit

    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*a, **k):  # type: ignore
        def deco(f):
            return f

        return deco if not a else a[0]


__all__ = ["ThreeBodyAngular", "angular_penalty", "MAX_NEIGHBOURS"]

#: Upper bound on neighbours inside the three-body cutoff.  The cutoff is set
#: between the first and second shells, so 6 is the physical maximum and this
#: leaves room for a compressed defect.
MAX_NEIGHBOURS = 12


# --------------------------------------------------------------------------- #
# Angular penalty
# --------------------------------------------------------------------------- #


def angular_penalty(cos_theta, n: int):
    """``P = (1 - cos(n*theta))/2`` and ``dP/d(cos theta)``, as arrays.

    Evaluated through the Chebyshev recurrences ``cos(n t) = T_n(c)`` and
    ``dT_n/dc = n U_{n-1}(c)``, so no inverse trigonometry is needed anywhere and
    the result is a polynomial -- cheap and free of branch cuts at 0 and pi.
    """
    c = np.clip(np.asarray(cos_theta, dtype=np.float64), -1.0, 1.0)
    t_prev, t_cur = np.ones_like(c), c.copy()
    u_prev, u_cur = np.ones_like(c), 2.0 * c
    for _ in range(n - 1):
        t_prev, t_cur = t_cur, 2.0 * c * t_cur - t_prev
        u_prev, u_cur = u_cur, 2.0 * c * u_cur - u_prev
    if n == 0:
        return np.zeros_like(c), np.zeros_like(c)
    p = 0.5 * (1.0 - t_cur)
    dp = -0.5 * n * u_prev
    return p, dp


@njit(cache=True, inline="always")  # pragma: no cover - compiled
def _penalty_scalar(c, n):
    if c > 1.0:
        c = 1.0
    elif c < -1.0:
        c = -1.0
    t_prev = 1.0
    t_cur = c
    u_prev = 1.0
    u_cur = 2.0 * c
    for _ in range(n - 1):
        t_new = 2.0 * c * t_cur - t_prev
        t_prev = t_cur
        t_cur = t_new
        u_new = 2.0 * c * u_cur - u_prev
        u_prev = u_cur
        u_cur = u_new
    return 0.5 * (1.0 - t_cur), -0.5 * n * u_prev


@njit(cache=True, inline="always")  # pragma: no cover - compiled
def _switch(r, r1, r2):
    """C^1 radial envelope: 1 below ``r1``, 0 above ``r2``.  Returns ``(f, df/dr)``."""
    if r <= r1:
        return 1.0, 0.0
    if r >= r2:
        return 0.0, 0.0
    x = math.pi * (r - r1) / (r2 - r1)
    return 0.5 * (1.0 + math.cos(x)), -0.5 * math.pi / (r2 - r1) * math.sin(x)


def _switch_np(r, r1, r2):
    r = np.asarray(r, dtype=np.float64)
    x = np.pi * (r - r1) / (r2 - r1)
    f = np.where(r <= r1, 1.0, np.where(r >= r2, 0.0, 0.5 * (1.0 + np.cos(x))))
    df = np.where((r > r1) & (r < r2), -0.5 * np.pi / (r2 - r1) * np.sin(x), 0.0)
    return f, df


# --------------------------------------------------------------------------- #
# Kernel
# --------------------------------------------------------------------------- #


@njit(cache=True)  # pragma: no cover - compiled
def _kernel(pos, lx, ly, pair_i, pair_j, lam, n_order, r1, r2):
    n = pos.shape[0]
    nb = np.full((n, MAX_NEIGHBOURS), -1, dtype=np.int64)
    cnt = np.zeros(n, dtype=np.int64)

    r2sq = r2 * r2
    for k in range(pair_i.shape[0]):
        i = pair_i[k]
        j = pair_j[k]
        dx = pos[i, 0] - pos[j, 0]
        dy = pos[i, 1] - pos[j, 1]
        dx -= lx * np.round(dx / lx)
        dy -= ly * np.round(dy / ly)
        if dx * dx + dy * dy < r2sq:
            if cnt[i] < MAX_NEIGHBOURS:
                nb[i, cnt[i]] = j
                cnt[i] += 1
            if cnt[j] < MAX_NEIGHBOURS:
                nb[j, cnt[j]] = i
                cnt[j] += 1

    forces = np.zeros((n, 2))
    per_atom = np.zeros(n)
    energy = 0.0
    virial = 0.0

    for i in range(n):
        for a in range(cnt[i]):
            j = nb[i, a]
            ux = pos[j, 0] - pos[i, 0]
            uy = pos[j, 1] - pos[i, 1]
            ux -= lx * np.round(ux / lx)
            uy -= ly * np.round(uy / ly)
            ru = math.sqrt(ux * ux + uy * uy)
            if ru <= 0.0:
                continue
            fu, dfu = _switch(ru, r1, r2)
            if fu == 0.0:
                continue

            for b in range(a + 1, cnt[i]):
                k = nb[i, b]
                vx = pos[k, 0] - pos[i, 0]
                vy = pos[k, 1] - pos[i, 1]
                vx -= lx * np.round(vx / lx)
                vy -= ly * np.round(vy / ly)
                rv = math.sqrt(vx * vx + vy * vy)
                if rv <= 0.0:
                    continue
                fv, dfv = _switch(rv, r1, r2)
                if fv == 0.0:
                    continue

                c = (ux * vx + uy * vy) / (ru * rv)
                p, dp = _penalty_scalar(c, n_order)

                e = lam * fu * fv * p
                energy += e
                per_atom[i] += e / 3.0
                per_atom[j] += e / 3.0
                per_atom[k] += e / 3.0

                # dc/du = v/(ru rv) - c u/ru^2 ;  dc/dv = u/(ru rv) - c v/rv^2
                dcdux = vx / (ru * rv) - c * ux / (ru * ru)
                dcduy = vy / (ru * rv) - c * uy / (ru * ru)
                dcdvx = ux / (ru * rv) - c * vx / (rv * rv)
                dcdvy = uy / (ru * rv) - c * vy / (rv * rv)

                # dE/du and dE/dv
                gu = lam * (dfu * fv * p * ux / ru + fu * fv * dp * dcdux)
                guy = lam * (dfu * fv * p * uy / ru + fu * fv * dp * dcduy)
                gv = lam * (dfv * fu * p * vx / rv + fu * fv * dp * dcdvx)
                gvy = lam * (dfv * fu * p * vy / rv + fu * fv * dp * dcdvy)

                # force on j is -dE/du, on k is -dE/dv, on i the balance
                forces[j, 0] -= gu
                forces[j, 1] -= guy
                forces[k, 0] -= gv
                forces[k, 1] -= gvy
                forces[i, 0] += gu + gv
                forces[i, 1] += guy + gvy

                # W = sum over the independent vectors of r . f
                virial += -(ux * gu + uy * guy) - (vx * gv + vy * gvy)

    return forces, energy, virial, per_atom


def _kernel_numpy(pos, box: Box, pair_i, pair_j, lam, n_order, r1, r2):
    """Readable reference implementation; loops atoms in Python."""
    n = pos.shape[0]
    d = box.minimum_image(pos[pair_i] - pos[pair_j])
    inside = (d[:, 0] ** 2 + d[:, 1] ** 2) < r2 * r2
    ii, jj = pair_i[inside], pair_j[inside]
    adjacency = [[] for _ in range(n)]
    for i, j in zip(ii, jj):
        adjacency[i].append(j)
        adjacency[j].append(i)

    forces = np.zeros((n, 2))
    per_atom = np.zeros(n)
    energy = 0.0
    virial = 0.0
    for i in range(n):
        nbs = adjacency[i]
        for a in range(len(nbs)):
            for b in range(a + 1, len(nbs)):
                j, k = nbs[a], nbs[b]
                u = box.minimum_image((pos[j] - pos[i])[None, :])[0]
                v = box.minimum_image((pos[k] - pos[i])[None, :])[0]
                ru, rv = np.hypot(*u), np.hypot(*v)
                fu, dfu = (float(x) for x in _switch_np(ru, r1, r2))
                fv, dfv = (float(x) for x in _switch_np(rv, r1, r2))
                if fu == 0.0 or fv == 0.0:
                    continue
                c = float(np.dot(u, v) / (ru * rv))
                p, dp = angular_penalty(c, n_order)
                p, dp = float(p), float(dp)

                e = lam * fu * fv * p
                energy += e
                per_atom[[i, j, k]] += e / 3.0

                dcdu = v / (ru * rv) - c * u / ru**2
                dcdv = u / (ru * rv) - c * v / rv**2
                gu = lam * (dfu * fv * p * u / ru + fu * fv * dp * dcdu)
                gv = lam * (dfv * fu * p * v / rv + fu * fv * dp * dcdv)

                forces[j] -= gu
                forces[k] -= gv
                forces[i] += gu + gv
                virial += -float(np.dot(u, gu)) - float(np.dot(v, gv))
    return forces, energy, virial, per_atom


# --------------------------------------------------------------------------- #


@dataclass
class ThreeBodyAngular:
    """Stillinger-Weber-style angular term.

    Parameters
    ----------
    strength
        ``lambda`` in eV.  This is the energy cost of putting a bond at the worst
        possible angle, with both radial envelopes at 1.  Values comparable to
        the pair well depth are what it takes to hold an open lattice together.
    order
        ``n`` in the penalty ``[1 - cos(n theta)]/2``.  Must match the target
        structure; :class:`~trilattice.lattices.LatticeSpec` carries the right
        one as ``angular_order``.
    r_inner, r_outer
        The radial envelope goes smoothly from 1 to 0 between these radii.
        ``r_outer`` must sit below the second neighbour shell, or second
        neighbours enter the triplet sum and the penalty stops describing the
        structure it was chosen for.
    """

    strength: float
    order: int
    r_inner: float
    r_outer: float

    def __post_init__(self) -> None:
        if self.order < 1:
            raise ValueError("angular order must be >= 1")
        if not (0.0 < self.r_inner < self.r_outer):
            raise ValueError("need 0 < r_inner < r_outer")

    @property
    def cutoff(self) -> float:
        return self.r_outer

    def __call__(self, positions, box: Box, pair_i, pair_j, use_numba: bool = True):
        args = (np.ascontiguousarray(positions), box.lx, box.ly,
                np.ascontiguousarray(pair_i), np.ascontiguousarray(pair_j),
                float(self.strength), int(self.order),
                float(self.r_inner), float(self.r_outer))
        if use_numba and HAVE_NUMBA:
            f, e, w, pa = _kernel(*args)
            return f, float(e), float(w), pa
        return _kernel_numpy(positions, box, pair_i, pair_j, float(self.strength),
                             int(self.order), float(self.r_inner), float(self.r_outer))

    @classmethod
    def for_lattice(cls, spec, a: float, strength: float = 1.0) -> "ThreeBodyAngular":
        """Build the term that stabilises ``spec``'s structure at spacing ``a``."""
        r_outer = spec.cutoff_factor * a
        return cls(strength=strength, order=spec.angular_order,
                   r_inner=0.85 * r_outer, r_outer=r_outer)

    def summary(self) -> str:  # pragma: no cover - cosmetic
        return (f"3-body angular(lambda={self.strength:.3f} eV, n={self.order}, "
                f"r={self.r_inner:.2f}-{self.r_outer:.2f} A)")
