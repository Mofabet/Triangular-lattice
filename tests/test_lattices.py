"""Each lattice must be built correctly and diagnosed with its own symmetry."""
import numpy as np
import pytest

import trilattice as tl
from trilattice import observables as obs
from trilattice.lattices import (
    LATTICES,
    LATTICE_NAMES,
    build_lattice,
    cells_for,
    emptiest_points,
    lattice_spec,
)

A = 3.2


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_every_site_has_the_ideal_coordination_at_exactly_a(kind):
    spec = lattice_spec(kind)
    cfg = build_lattice(kind, 5, 4, A)
    i, j, d = obs.neighbour_pairs_within(cfg.box, cfg.positions, spec.cutoff(A))
    z = np.bincount(np.concatenate([i, j]), minlength=cfg.n_particles)
    assert np.all(z == spec.coordination)
    assert np.allclose(np.hypot(d[:, 0], d[:, 1]), A)


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_density_matches_the_tabulated_factor(kind):
    spec = lattice_spec(kind)
    cfg = build_lattice(kind, 5, 4, A)
    assert cfg.density * A**2 == pytest.approx(spec.density_factor, rel=1e-12)
    assert cfg.n_particles == 5 * 4 * spec.atoms_per_cell


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_the_box_matches_the_declared_cell_factors(kind):
    spec = lattice_spec(kind)
    cfg = build_lattice(kind, 5, 4, A)
    fx, fy = spec.cell_factors
    assert cfg.box.lx == pytest.approx(5 * fx * A)
    assert cfg.box.ly == pytest.approx(4 * fy * A)


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_the_specs_order_parameter_is_one_on_its_own_lattice(kind):
    spec = lattice_spec(kind)
    cfg = build_lattice(kind, 5, 4, A)
    value = obs.global_psi_n(cfg.positions, cfg.box, spec.cutoff(A), spec.psi_order)
    assert value == pytest.approx(1.0, abs=1e-9)


def test_the_hexatic_parameter_is_zero_on_a_square_lattice():
    """Why the order has to follow the lattice: psi_6 would call it molten."""
    cfg = build_lattice("square", 6, 6, A)
    cut = lattice_spec("square").cutoff(A)
    assert obs.global_psi_n(cfg.positions, cfg.box, cut, 6) == pytest.approx(0.0, abs=1e-9)
    assert obs.global_psi_n(cfg.positions, cfg.box, cut, 4) == pytest.approx(1.0, abs=1e-9)


def test_psi3_cancels_on_a_triangular_lattice_but_not_locally_on_honeycomb():
    """The sublattice subtlety that forces honeycomb onto n = 6 globally."""
    tri = build_lattice("triangular", 5, 4, A)
    hon = build_lattice("honeycomb", 4, 3, A)
    tri_cut = lattice_spec("triangular").cutoff(A)
    hon_cut = lattice_spec("honeycomb").cutoff(A)
    # triangular: the six bonds cancel in threes
    assert obs.global_psi_n(tri.positions, tri.box, tri_cut, 3) == pytest.approx(0.0, abs=1e-9)
    # honeycomb: every site has |psi_3| = 1 ...
    local = np.abs(obs.psi_n(hon.positions, hon.box, hon_cut, 3))
    assert np.allclose(local, 1.0, atol=1e-9)
    # ... but the two sublattices are 180 degrees apart in phase, so it averages away
    assert obs.global_psi_n(hon.positions, hon.box, hon_cut, 3) == pytest.approx(0.0, abs=1e-9)
    assert lattice_spec("honeycomb").psi_order == 6


def test_kagome_and_triangular_are_indistinguishable_by_psi6_alone():
    """The false positive that makes the coordination check necessary."""
    kag = build_lattice("kagome", 4, 3, A)
    tri = build_lattice("triangular", 5, 4, A)
    k_cut = lattice_spec("kagome").cutoff(A)
    t_cut = lattice_spec("triangular").cutoff(A)
    assert obs.global_psi_n(kag.positions, kag.box, k_cut, 6) == pytest.approx(1.0, abs=1e-9)
    assert obs.global_psi_n(tri.positions, tri.box, t_cut, 6) == pytest.approx(1.0, abs=1e-9)
    # only the coordination separates them
    assert obs.coordination_by_cutoff(kag.positions, kag.box, k_cut).mean() == 4
    assert obs.coordination_by_cutoff(tri.positions, tri.box, t_cut).mean() == 6


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_cells_for_keeps_the_particle_count_near_the_target(kind):
    nx, ny = cells_for(kind, 300)
    n = build_lattice(kind, nx, ny, A).n_particles
    assert 0.75 * 300 <= n <= 1.35 * 300


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_cells_for_keeps_the_box_roughly_square(kind):
    nx, ny = cells_for(kind, 300)
    box = build_lattice(kind, nx, ny, A).box
    assert 0.7 < box.lx / box.ly < 1.45


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_defect_fraction_is_zero_on_the_ideal_lattice(kind):
    spec = lattice_spec(kind)
    cfg = build_lattice(kind, 5, 4, A)
    cut = None if spec.use_delaunay() else spec.cutoff(A)
    assert obs.defect_fraction(cfg.positions, cfg.box, spec.coordination, cut) == 0.0


def test_delaunay_is_only_claimed_for_the_triangular_lattice():
    assert lattice_spec("triangular").use_delaunay()
    for kind in ("square", "honeycomb", "kagome"):
        assert not lattice_spec(kind).use_delaunay()


def test_only_triangular_is_declared_a_ground_state():
    assert lattice_spec("triangular").lj_ground_state
    assert sum(s.lj_ground_state for s in LATTICES.values()) == 1


def test_triangular_is_the_lowest_energy_of_the_four():
    """The claim the whole module rests on, checked rather than asserted."""
    pot = tl.LennardJones(0.27, 2.88, 8.0, mode="cut")
    energies = {}
    for kind in LATTICE_NAMES:
        nx, ny = cells_for(kind, 300)
        cfg = build_lattice(kind, nx, ny, 3.3567)
        sim = tl.Simulation(cfg, pot, timestep=0.002, seed=1)
        energies[kind] = sim.e_pot / cfg.n_particles
    assert energies["triangular"] == min(energies.values())


def test_emptiest_points_finds_the_hollow_of_a_triangular_lattice():
    cfg = build_lattice("triangular", 6, 4, A)
    p = emptiest_points(cfg.positions, cfg.box, n=1, grid=240)[0]
    d = cfg.box.minimum_image(cfg.positions - p)
    assert np.min(np.hypot(d[:, 0], d[:, 1])) == pytest.approx(A / np.sqrt(3), rel=0.03)


def test_emptiest_points_spreads_several_picks_apart():
    cfg = build_lattice("square", 8, 8, A)
    pts = emptiest_points(cfg.positions, cfg.box, n=4, grid=160)
    d = cfg.box.minimum_image(pts[:, None, :] - pts[None, :, :])
    r = np.hypot(d[:, :, 0], d[:, :, 1])
    assert np.min(r[~np.eye(4, dtype=bool)]) > A


def test_unknown_lattice_is_rejected():
    with pytest.raises(ValueError, match="unknown lattice"):
        build_lattice("hexagonal-ish", 4, 4, A)
