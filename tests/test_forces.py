import numpy as np
import pytest

from trilattice import LennardJones, triangular_lattice
from trilattice.forces import ForceField, HAVE_NUMBA, numerical_forces
from trilattice.neighbors import all_pairs_within


def _system(seed=0, jitter=0.25):
    rng = np.random.default_rng(seed)
    cfg = triangular_lattice(6, 4, 3.2)
    cfg.positions += jitter * rng.standard_normal(cfg.positions.shape)
    return cfg


@pytest.mark.parametrize("mode", ["cut", "shifted", "shifted-force"])
def test_analytic_forces_match_central_differences(mode):
    cfg = _system()
    pot = LennardJones(0.27, 2.88, 8.0, mode=mode)
    i, j = all_pairs_within(cfg.box, cfg.box.wrap(cfg.positions), pot.cutoff)
    f = ForceField(pot, cfg.box, use_numba=False)(cfg.positions, i, j, cfg.types).forces
    assert np.allclose(f, numerical_forces(pot, cfg.box, cfg.positions, cfg.types), atol=1e-6)


def test_newtons_third_law():
    cfg = _system(1)
    pot = LennardJones(0.27, 2.88, 8.0)
    i, j = all_pairs_within(cfg.box, cfg.box.wrap(cfg.positions), pot.cutoff)
    f = ForceField(pot, cfg.box, use_numba=False)(cfg.positions, i, j, cfg.types).forces
    assert np.allclose(f.sum(axis=0), 0.0, atol=1e-10)


@pytest.mark.skipif(not HAVE_NUMBA, reason="numba not installed")
@pytest.mark.parametrize("mode", ["cut", "shifted", "shifted-force"])
def test_numba_and_numpy_kernels_agree(mode):
    cfg = _system(2)
    pot = LennardJones(0.27, 2.88, 8.0, mode=mode)
    i, j = all_pairs_within(cfg.box, cfg.box.wrap(cfg.positions), pot.cutoff)
    a = ForceField(pot, cfg.box, use_numba=False)(cfg.positions, i, j, cfg.types)
    b = ForceField(pot, cfg.box, use_numba=True)(cfg.positions, i, j, cfg.types)
    assert np.allclose(a.forces, b.forces, atol=1e-12)
    assert a.energy == pytest.approx(b.energy, abs=1e-12)
    assert a.virial == pytest.approx(b.virial, abs=1e-10)


def test_forces_vanish_on_the_perfect_lattice_at_a0():
    # At the lattice sum minimum every atom sits at a symmetry point, so the
    # force must vanish identically -- a strong end-to-end check of PBC, the
    # neighbour list and the kernel at once.
    cfg = triangular_lattice(8, 5, 3.2010)
    pot = LennardJones(0.27, 2.88, 8.0, mode="cut")
    i, j = all_pairs_within(cfg.box, cfg.positions, pot.cutoff)
    f = ForceField(pot, cfg.box, use_numba=False)(cfg.positions, i, j, cfg.types).forces
    assert np.max(np.abs(f)) < 1e-12


def test_per_atom_energies_sum_to_the_total():
    cfg = _system(3)
    pot = LennardJones(0.27, 2.88, 8.0)
    i, j = all_pairs_within(cfg.box, cfg.box.wrap(cfg.positions), pot.cutoff)
    r = ForceField(pot, cfg.box, use_numba=False)(cfg.positions, i, j, cfg.types)
    assert r.per_atom_energy.sum() == pytest.approx(r.energy)


def test_energy_is_translation_and_pbc_invariant():
    cfg = _system(4)
    pot = LennardJones(0.27, 2.88, 8.0)

    def energy(p):
        w = cfg.box.wrap(p)
        i, j = all_pairs_within(cfg.box, w, pot.cutoff)
        return ForceField(pot, cfg.box, use_numba=False)(w, i, j, cfg.types).energy

    shifted = cfg.positions + np.array([17.3, -41.9])
    assert energy(cfg.positions) == pytest.approx(energy(shifted), rel=1e-12)
