import numpy as np
import pytest

from trilattice import Box, triangular_lattice, add_defects, perfect_lattice_neighbour_shells
from trilattice.observables import neighbour_pairs_within


def test_every_site_has_six_neighbours_at_exactly_a():
    a = 3.2
    cfg = triangular_lattice(6, 4, a)
    i, j, d = neighbour_pairs_within(cfg.box, cfg.positions, 1.01 * a)
    counts = np.bincount(np.concatenate([i, j]), minlength=cfg.n_particles)
    assert np.all(counts == 6)
    assert np.allclose(np.hypot(d[:, 0], d[:, 1]), a)


def test_box_is_commensurate_with_the_lattice():
    a, nx, ny = 3.2, 5, 3
    cfg = triangular_lattice(nx, ny, a)
    assert cfg.box.lx == pytest.approx(nx * a)
    assert cfg.box.ly == pytest.approx(ny * a * np.sqrt(3))
    assert cfg.n_particles == 2 * nx * ny


def test_density_matches_the_analytic_value():
    a = 3.2
    cfg = triangular_lattice(7, 5, a)
    assert cfg.density == pytest.approx(2.0 / (a**2 * np.sqrt(3)))


def test_minimum_image_handles_displacements_larger_than_the_box():
    box = Box(10.0, 12.0)
    d = np.array([[26.0, -31.0], [4.9, 5.9], [5.1, 6.1]])
    m = box.minimum_image(d)
    assert np.all(np.abs(m[:, 0]) <= 5.0 + 1e-12)
    assert np.all(np.abs(m[:, 1]) <= 6.0 + 1e-12)
    # the naive `if dx > L: dx -= L` idiom fails on the first row; rounding does not
    assert m[0] == pytest.approx([-4.0, 5.0])


def test_wrap_is_idempotent_and_in_range():
    box = Box(10.0, 12.0)
    rng = np.random.default_rng(0)
    p = rng.uniform(-40, 40, size=(200, 2))
    w = box.wrap(p)
    assert np.all((w >= 0) & (w < box.lengths))
    assert np.allclose(w, box.wrap(w))


def test_vacancies_and_interstitials_change_the_count_correctly():
    cfg = triangular_lattice(6, 4, 3.2)
    n = cfg.n_particles
    rng = np.random.default_rng(1)
    out = add_defects(cfg, n_vacancies=3, n_interstitials=2, a=3.2, rng=rng, min_separation=1.5)
    assert out.n_particles == n - 3 + 2


def test_min_separation_above_the_hollow_site_radius_is_refused():
    """A three-fold hollow sits at a/sqrt(3); nothing can be placed further out."""
    cfg = triangular_lattice(6, 4, 3.2)
    with pytest.raises(RuntimeError, match="only placed"):
        add_defects(cfg, n_interstitials=1, a=3.2,
                    rng=np.random.default_rng(0), min_separation=3.2 / np.sqrt(3) + 0.01)


def test_interstitials_never_overlap_a_host_atom():
    cfg = triangular_lattice(8, 5, 3.2)
    rng = np.random.default_rng(2)
    out = add_defects(cfg, n_interstitials=5, a=3.2, rng=rng, min_separation=1.5)
    i, j = np.triu_indices(out.n_particles, 1)
    d = out.box.minimum_image(out.positions[i] - out.positions[j])
    assert np.min(np.hypot(d[:, 0], d[:, 1])) >= 1.5


def test_shell_radii_follow_the_triangular_norm_form():
    radii, counts = perfect_lattice_neighbour_shells(1.0, 2.1)
    assert radii[:3] == pytest.approx([1.0, np.sqrt(3), 2.0])
    assert counts[0] == 6 and counts[1] == 6 and counts[2] == 6


def test_replicate_preserves_density_and_lattice_perfection():
    from trilattice import replicate

    cfg = triangular_lattice(4, 3, 3.2)
    big = replicate(cfg, 3, 2)
    assert big.n_particles == 6 * cfg.n_particles
    assert big.density == pytest.approx(cfg.density)
    i, j, d = neighbour_pairs_within(big.box, big.positions, 1.01 * 3.2)
    counts = np.bincount(np.concatenate([i, j]), minlength=big.n_particles)
    assert np.all(counts == 6)


def test_required_replication_satisfies_minimum_image():
    from trilattice import replicate, required_replication

    cfg = triangular_lattice(5, 5, 3.2)          # the original 50-atom deck
    r_list = 9.0
    assert cfg.box.max_cutoff < r_list           # the shipped deck is too small
    nx, ny = required_replication(cfg.box, r_list)
    assert replicate(cfg, nx, ny).box.max_cutoff >= r_list
