import numpy as np
import pytest

from trilattice import NeighborList, triangular_lattice
from trilattice.neighbors import all_pairs_within


def _as_set(i, j):
    return set(map(tuple, np.sort(np.stack([i, j], axis=1), axis=1)))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_cell_list_finds_exactly_the_brute_force_pairs(seed):
    rng = np.random.default_rng(seed)
    cfg = triangular_lattice(14, 8, 3.2)
    cfg.positions += 0.4 * rng.standard_normal(cfg.positions.shape)
    nl = NeighborList(cfg.box, cutoff=8.0, skin=1.0)
    nl.build(cfg.positions)
    ref = all_pairs_within(cfg.box, cfg.box.wrap(cfg.positions), nl.r_list)
    assert _as_set(*nl.pairs) == _as_set(*ref)


def test_no_pair_is_listed_twice():
    cfg = triangular_lattice(12, 7, 3.2)
    nl = NeighborList(cfg.box, 8.0, skin=1.0)
    nl.build(cfg.positions)
    i, j = nl.pairs
    assert len(_as_set(i, j)) == i.size


def test_small_box_falls_back_to_brute_force_without_double_counting():
    cfg = triangular_lattice(3, 2, 3.2)          # box smaller than 3 r_list
    nl = NeighborList(cfg.box, cutoff=4.0, skin=0.5)
    nl.build(cfg.positions)
    ref = all_pairs_within(cfg.box, cfg.box.wrap(cfg.positions), nl.r_list)
    assert _as_set(*nl.pairs) == _as_set(*ref)


def test_cutoff_larger_than_half_the_box_is_rejected():
    cfg = triangular_lattice(3, 2, 3.2)
    with pytest.raises(ValueError, match="minimum-image"):
        NeighborList(cfg.box, cutoff=cfg.box.ly, skin=1.0)


def test_rebuild_is_triggered_only_after_half_the_skin():
    cfg = triangular_lattice(10, 6, 3.2)
    nl = NeighborList(cfg.box, 8.0, skin=1.0)
    nl.build(cfg.positions)
    p = cfg.positions.copy()
    p[0, 0] += 0.4
    assert not nl.needs_rebuild(p)
    p[0, 0] += 0.4
    assert nl.needs_rebuild(p)


def test_pair_count_is_linear_in_n():
    counts = []
    for nx, ny in [(10, 6), (20, 12), (30, 18)]:
        cfg = triangular_lattice(nx, ny, 3.2)
        nl = NeighborList(cfg.box, 8.0, skin=1.0)
        nl.build(cfg.positions)
        counts.append(nl.n_pairs / cfg.n_particles)
    assert max(counts) - min(counts) < 0.5     # pairs per atom is intensive
