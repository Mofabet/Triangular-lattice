"""The angular term: correctness first, then whether it does its job."""
import numpy as np
import pytest

import trilattice as tl
from trilattice import observables as obs
from trilattice.forces import ForceField
from trilattice.lattices import LATTICE_NAMES, build_lattice, cells_for, lattice_spec
from trilattice.neighbors import all_pairs_within
from trilattice.threebody import HAVE_NUMBA, ThreeBodyAngular, angular_penalty

A = 3.3567
POT = tl.LennardJones(0.27, 2.88, 8.0, mode="shifted-force")


# --- the penalty ----------------------------------------------------------- #
@pytest.mark.parametrize("n,angles", [(3, [120]), (4, [90, 180]), (6, [60, 120, 180])])
def test_penalty_vanishes_at_the_lattice_angles(n, angles):
    p, _ = angular_penalty(np.cos(np.radians(angles)), n)
    assert np.allclose(p, 0.0, atol=1e-12)


@pytest.mark.parametrize("n", [3, 4, 6])
def test_penalty_is_bounded_and_peaks_between_the_zeros(n):
    theta = np.linspace(0, np.pi, 2001)
    p, _ = angular_penalty(np.cos(theta), n)
    assert p.min() >= -1e-12 and p.max() <= 1.0 + 1e-12
    assert p.max() == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize("n", [3, 4, 6])
def test_penalty_derivative_matches_finite_differences(n):
    c = np.linspace(-0.98, 0.98, 97)
    h = 1e-6
    _, dp = angular_penalty(c, n)
    numeric = (angular_penalty(c + h, n)[0] - angular_penalty(c - h, n)[0]) / (2 * h)
    assert np.allclose(dp, numeric, atol=1e-5)


# --- the term on ideal lattices -------------------------------------------- #
@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_the_ideal_lattice_costs_exactly_nothing(kind):
    """The whole design rests on this: no penalty for the wanted geometry, so
    the term shifts no equilibrium property and can be toggled freely."""
    spec = lattice_spec(kind)
    cfg = build_lattice(kind, 4, 3, A)
    term = ThreeBodyAngular.for_lattice(spec, A, strength=1.0)
    i, j = all_pairs_within(cfg.box, cfg.positions, term.cutoff)
    forces, energy, _, _ = term(cfg.positions, cfg.box, i, j)
    assert energy == pytest.approx(0.0, abs=1e-12)
    assert np.max(np.abs(forces)) < 1e-10


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_the_term_does_not_change_the_pair_energy_of_the_ideal_lattice(kind):
    spec = lattice_spec(kind)
    nx, ny = cells_for(kind, 200)      # big enough for the 8 A pair cutoff
    cfg = build_lattice(kind, nx, ny, A)
    plain = tl.Simulation(cfg, POT, timestep=0.002, seed=1)
    with_term = tl.Simulation(cfg, POT, timestep=0.002, seed=1,
                              three_body=ThreeBodyAngular.for_lattice(spec, A, 1.0))
    assert with_term.e_pot == pytest.approx(plain.e_pot, abs=1e-10)


# --- forces ---------------------------------------------------------------- #
def _perturbed(kind, seed=0, amp=0.12):
    rng = np.random.default_rng(seed)
    cfg = build_lattice(kind, 4, 3, A)
    cfg.positions = cfg.box.wrap(cfg.positions + amp * rng.standard_normal(cfg.positions.shape))
    return cfg


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_three_body_forces_match_central_differences(kind):
    spec = lattice_spec(kind)
    cfg = _perturbed(kind)
    term = ThreeBodyAngular.for_lattice(spec, A, strength=1.0)

    def energy(q):
        w = cfg.box.wrap(q)
        i, j = all_pairs_within(cfg.box, w, term.cutoff)
        return term(w, cfg.box, i, j)[1]

    i, j = all_pairs_within(cfg.box, cfg.positions, term.cutoff)
    analytic = term(cfg.positions, cfg.box, i, j)[0]

    h = 1e-6
    numeric = np.zeros_like(analytic)
    for atom in range(cfg.n_particles):
        for comp in range(2):
            plus, minus = cfg.positions.copy(), cfg.positions.copy()
            plus[atom, comp] += h
            minus[atom, comp] -= h
            numeric[atom, comp] = -(energy(plus) - energy(minus)) / (2 * h)
    assert np.allclose(analytic, numeric, atol=1e-6)


@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_three_body_forces_sum_to_zero(kind):
    cfg = _perturbed(kind, seed=1)
    term = ThreeBodyAngular.for_lattice(lattice_spec(kind), A, 1.0)
    i, j = all_pairs_within(cfg.box, cfg.positions, term.cutoff)
    assert np.allclose(term(cfg.positions, cfg.box, i, j)[0].sum(axis=0), 0.0, atol=1e-10)


@pytest.mark.skipif(not HAVE_NUMBA, reason="numba not installed")
@pytest.mark.parametrize("kind", LATTICE_NAMES)
def test_numba_and_numpy_three_body_kernels_agree(kind):
    cfg = _perturbed(kind, seed=2)
    term = ThreeBodyAngular.for_lattice(lattice_spec(kind), A, 1.0)
    i, j = all_pairs_within(cfg.box, cfg.positions, term.cutoff)
    fa, ea, wa, pa = term(cfg.positions, cfg.box, i, j, use_numba=True)
    fb, eb, wb, pb = term(cfg.positions, cfg.box, i, j, use_numba=False)
    assert np.allclose(fa, fb, atol=1e-10)
    assert ea == pytest.approx(eb, abs=1e-12)
    assert wa == pytest.approx(wb, abs=1e-10)
    assert np.allclose(pa, pb, atol=1e-12)


def test_per_atom_energies_sum_to_the_total():
    cfg = _perturbed("honeycomb", seed=3)
    term = ThreeBodyAngular.for_lattice(lattice_spec("honeycomb"), A, 1.0)
    i, j = all_pairs_within(cfg.box, cfg.positions, term.cutoff)
    _, e, _, per_atom = term(cfg.positions, cfg.box, i, j)
    assert per_atom.sum() == pytest.approx(e, rel=1e-12)


def test_the_forcefield_sums_pair_and_three_body_contributions():
    cfg = _perturbed("square", seed=4)
    term = ThreeBodyAngular.for_lattice(lattice_spec("square"), A, 1.0)
    i, j = all_pairs_within(cfg.box, cfg.positions, POT.cutoff)
    pair = ForceField(POT, cfg.box)(cfg.positions, i, j, cfg.types)
    both = ForceField(POT, cfg.box, three_body=term)(cfg.positions, i, j, cfg.types)
    only3 = term(cfg.positions, cfg.box, i, j)
    assert both.energy == pytest.approx(pair.energy + only3[1], rel=1e-12)
    assert np.allclose(both.forces, pair.forces + only3[0], atol=1e-12)


def test_energy_is_conserved_in_nve_with_the_angular_term():
    """A force that is not the gradient of the energy shows up here first."""
    spec = lattice_spec("honeycomb")
    nx, ny = cells_for("honeycomb", 200)
    cfg = build_lattice("honeycomb", nx, ny, A)
    sim = tl.Simulation(cfg, POT, timestep=0.001, seed=5,
                        three_body=ThreeBodyAngular.for_lattice(spec, A, 0.5))
    sim.set_temperature(200.0)
    sim.run(2000, log_every=0)
    sim.clear_log()
    log = sim.run(8000, log_every=20)
    drift = np.std(np.array(log.e_tot)) / sim.n_particles
    assert drift < 1e-4


# --- does it actually hold the lattices? ----------------------------------- #
@pytest.mark.parametrize("kind", ["square", "honeycomb", "kagome"])
def test_the_open_lattices_survive_with_the_term_and_collapse_without(kind):
    spec = lattice_spec(kind)
    nx, ny = cells_for(kind, 200)
    cut = spec.cutoff(A)
    results = {}
    for label, term in (("pair", None),
                        ("3body", ThreeBodyAngular.for_lattice(spec, A, 0.5))):
        sim = tl.Simulation(build_lattice(kind, nx, ny, A), POT, timestep=0.002,
                            three_body=term,
                            thermostat=tl.Bussi(temperature=300.0, tau=0.2), seed=3)
        sim.set_temperature(300.0)
        sim.run(2500, log_every=0)
        results[label] = (
            obs.global_psi_n(sim.positions, sim.box, cut, spec.psi_order),
            obs.coordination_by_cutoff(sim.positions, sim.box, cut).mean(),
        )
    psi_on, z_on = results["3body"]
    psi_off, z_off = results["pair"]
    assert psi_on > 0.9, f"{kind}: angular term failed to hold the structure"
    assert z_on == pytest.approx(spec.coordination, abs=0.1)
    assert z_off > spec.coordination + 0.5, f"{kind}: expected collapse without the term"


def test_triangular_needs_no_help():
    spec = lattice_spec("triangular")
    nx, ny = cells_for("triangular", 200)
    sim = tl.Simulation(build_lattice("triangular", nx, ny, A), POT, timestep=0.002,
                        thermostat=tl.Bussi(temperature=300.0, tau=0.2), seed=3)
    sim.set_temperature(300.0)
    sim.run(2500, log_every=0)
    assert obs.global_psi_n(sim.positions, sim.box, spec.cutoff(A), 6) > 0.9


def test_rejects_a_bad_parametrisation():
    with pytest.raises(ValueError):
        ThreeBodyAngular(strength=1.0, order=0, r_inner=1.0, r_outer=2.0)
    with pytest.raises(ValueError):
        ThreeBodyAngular(strength=1.0, order=3, r_inner=2.0, r_outer=1.0)
