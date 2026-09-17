"""Force, energy and virial evaluation over a Verlet pair list.

Everything is accumulated with :func:`numpy.bincount`, which is the fastest way
in pure NumPy to scatter-add per-pair contributions onto per-particle arrays
(``np.add.at`` is correct but roughly an order of magnitude slower).

If Numba is installed an equivalent compiled kernel is used instead; it is
typically 5-10x faster because it avoids materialising the per-pair temporaries.
The two paths are checked against each other in the test suite.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .lattice import Box
from .potentials import LennardJones

try:  # pragma: no cover - environment dependent
    from numba import njit

    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*a, **k):  # type: ignore
        def deco(f):
            return f

        return deco if not a else a[0]


@dataclass
class ForceResult:
    """Everything one force evaluation produces."""

    forces: np.ndarray  #: (N, 2) eV/A
    energy: float  #: total potential energy, eV
    virial: float  #: sum_{i<j} r_ij . f_ij, eV
    per_atom_energy: np.ndarray  #: (N,) eV


# --------------------------------------------------------------------------- #
# NumPy kernel
# --------------------------------------------------------------------------- #


def _evaluate_numpy(
    positions: np.ndarray,
    box: Box,
    pair_i: np.ndarray,
    pair_j: np.ndarray,
    types: np.ndarray,
    eps_t: np.ndarray,
    sig2_t: np.ndarray,
    eshift_t: np.ndarray,
    fcut_t: np.ndarray,
    cutoff: float,
) -> ForceResult:
    n = positions.shape[0]
    forces = np.zeros((n, 2), dtype=np.float64)
    per_atom = np.zeros(n, dtype=np.float64)
    if pair_i.size == 0:
        return ForceResult(forces, 0.0, 0.0, per_atom)

    d = box.minimum_image(positions[pair_i] - positions[pair_j])
    r2 = d[:, 0] ** 2 + d[:, 1] ** 2
    inside = r2 < cutoff * cutoff
    if not np.any(inside):
        return ForceResult(forces, 0.0, 0.0, per_atom)

    d = d[inside]
    r2 = r2[inside]
    ii = pair_i[inside]
    jj = pair_j[inside]

    if eps_t.shape[0] == 1:
        eps = eps_t[0, 0]
        sig2 = sig2_t[0, 0]
        eshift = eshift_t[0, 0]
        fcut = fcut_t[0, 0]
    else:
        ti, tj = types[ii], types[jj]
        eps = eps_t[ti, tj]
        sig2 = sig2_t[ti, tj]
        eshift = eshift_t[ti, tj]
        fcut = fcut_t[ti, tj]

    inv_r2 = 1.0 / r2
    sr2 = sig2 * inv_r2
    sr6 = sr2 * sr2 * sr2
    sr12 = sr6 * sr6

    r = np.sqrt(r2)
    # scalar force f(r) = -du/dr, then fpair = f(r)/r so that F_vec = fpair * d
    fpair = 24.0 * eps * (2.0 * sr12 - sr6) * inv_r2
    energy = 4.0 * eps * (sr12 - sr6) + eshift
    if np.any(fcut):
        fpair = fpair - fcut / r
        energy = energy + r * fcut

    fx = fpair * d[:, 0]
    fy = fpair * d[:, 1]

    forces[:, 0] = np.bincount(ii, weights=fx, minlength=n) - np.bincount(jj, weights=fx, minlength=n)
    forces[:, 1] = np.bincount(ii, weights=fy, minlength=n) - np.bincount(jj, weights=fy, minlength=n)

    half = 0.5 * energy
    per_atom = np.bincount(ii, weights=half, minlength=n) + np.bincount(jj, weights=half, minlength=n)

    total_energy = float(np.sum(energy))
    virial = float(np.sum(fpair * r2))  # r_ij . f_ij summed over pairs
    return ForceResult(forces, total_energy, virial, per_atom)


# --------------------------------------------------------------------------- #
# Numba kernel
# --------------------------------------------------------------------------- #


@njit(cache=True, fastmath=True)
def _kernel_numba(pos, lx, ly, pi, pj, types, eps_t, sig2_t, eshift_t, fcut_t, rc2, multi):  # pragma: no cover
    n = pos.shape[0]
    forces = np.zeros((n, 2))
    per_atom = np.zeros(n)
    etot = 0.0
    virial = 0.0
    for k in range(pi.shape[0]):
        i = pi[k]
        j = pj[k]
        dx = pos[i, 0] - pos[j, 0]
        dy = pos[i, 1] - pos[j, 1]
        dx -= lx * np.round(dx / lx)
        dy -= ly * np.round(dy / ly)
        r2 = dx * dx + dy * dy
        if r2 >= rc2:
            continue
        if multi:
            ti = types[i]
            tj = types[j]
        else:
            ti = 0
            tj = 0
        eps = eps_t[ti, tj]
        sig2 = sig2_t[ti, tj]
        inv_r2 = 1.0 / r2
        sr2 = sig2 * inv_r2
        sr6 = sr2 * sr2 * sr2
        sr12 = sr6 * sr6
        fpair = 24.0 * eps * (2.0 * sr12 - sr6) * inv_r2
        e = 4.0 * eps * (sr12 - sr6) + eshift_t[ti, tj]
        fc = fcut_t[ti, tj]
        if fc != 0.0:
            r = np.sqrt(r2)
            fpair -= fc / r
            e += r * fc
        fx = fpair * dx
        fy = fpair * dy
        forces[i, 0] += fx
        forces[i, 1] += fy
        forces[j, 0] -= fx
        forces[j, 1] -= fy
        etot += e
        per_atom[i] += 0.5 * e
        per_atom[j] += 0.5 * e
        virial += fpair * r2
    return forces, etot, virial, per_atom


# --------------------------------------------------------------------------- #


class ForceField:
    """Binds a potential to a box and evaluates it on a pair list."""

    def __init__(self, potential: LennardJones, box: Box, *, use_numba: bool | None = None):
        self.potential = potential
        self.box = box
        self.use_numba = HAVE_NUMBA if use_numba is None else (use_numba and HAVE_NUMBA)
        self.n_evaluations = 0

    def __call__(self, positions, pair_i, pair_j, types) -> ForceResult:
        p = self.potential
        self.n_evaluations += 1
        if self.use_numba:
            f, e, w, pa = _kernel_numba(
                np.ascontiguousarray(positions),
                self.box.lx,
                self.box.ly,
                pair_i,
                pair_j,
                types.astype(np.int64),
                p.eps_table,
                p.sig2_table,
                p.energy_shift,
                p.force_at_cutoff,
                p.cutoff**2,
                p.n_types > 1,
            )
            return ForceResult(f, float(e), float(w), pa)
        return _evaluate_numpy(
            positions,
            self.box,
            pair_i,
            pair_j,
            types,
            p.eps_table,
            p.sig2_table,
            p.energy_shift,
            p.force_at_cutoff,
            p.cutoff,
        )


def numerical_forces(potential, box, positions, types, h: float = 1e-6) -> np.ndarray:
    """Central-difference forces -- the ground truth for the analytic kernel."""
    from .neighbors import all_pairs_within

    ff = ForceField(potential, box, use_numba=False)

    def energy_of(p):
        i, j = all_pairs_within(box, box.wrap(p), potential.cutoff)
        return ff(p, i, j, types).energy

    f = np.zeros_like(positions)
    for a in range(positions.shape[0]):
        for c in range(2):
            plus = positions.copy()
            minus = positions.copy()
            plus[a, c] += h
            minus[a, c] -= h
            f[a, c] = -(energy_of(plus) - energy_of(minus)) / (2.0 * h)
    return f


__all__ = ["ForceField", "ForceResult", "numerical_forces", "HAVE_NUMBA"]
