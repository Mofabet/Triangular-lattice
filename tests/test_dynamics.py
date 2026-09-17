import numpy as np
import pytest

from trilattice import Bussi, LennardJones, Simulation, triangular_lattice
from trilattice.thermostats import NoseHooverChain, Berendsen, Langevin
from trilattice.units import KB

POT = LennardJones(0.27, 2.88, 8.0, mode="shifted-force")


def _sim(**kw):
    cfg = triangular_lattice(8, 5, 3.2010)
    return Simulation(cfg, POT, timestep=0.002, seed=5, **kw)


def test_nve_conserves_energy():
    sim = _sim()
    sim.set_temperature(400.0)
    sim.run(2000, log_every=0)
    sim.clear_log()
    log = sim.run(10000, log_every=25)
    e = np.array(log.e_tot) / sim.n_particles
    assert np.std(e) < 1e-5            # eV/atom
    assert abs(e[-1] - e[0]) < 1e-4


def test_energy_error_scales_as_dt_squared():
    rms = []
    for dt in (0.001, 0.002, 0.004):
        cfg = triangular_lattice(8, 5, 3.2010)
        sim = Simulation(cfg, POT, timestep=dt, seed=5)
        sim.set_temperature(400.0)
        sim.run(2000, log_every=0)
        sim.clear_log()
        rms.append(np.std(sim.run(6000, log_every=10).e_tot))
    ratios = [rms[1] / rms[0], rms[2] / rms[1]]
    assert all(2.8 < r < 5.2 for r in ratios)


def test_momentum_stays_at_zero():
    sim = _sim()
    sim.set_temperature(500.0)
    log = sim.run(5000, log_every=50)
    assert max(log.momentum) < 1e-8


def test_time_reversibility():
    sim = _sim()
    sim.set_temperature(300.0)
    sim.run(500, log_every=0)
    r0, v0 = sim.positions.copy(), sim.velocities.copy()
    sim.run(500, log_every=0)
    sim.velocities *= -1.0
    sim.run(500, log_every=0)
    d = sim.box.minimum_image(sim.positions - r0)
    assert np.max(np.abs(d)) < 1e-6
    assert np.allclose(sim.velocities, -v0, atol=1e-6)


@pytest.mark.parametrize("name,kw", [
    ("bussi", dict(tau=0.1)),
    ("nose-hoover", dict(tau=0.1)),
    ("langevin", dict(friction=10.0)),
    ("berendsen", dict(tau=0.1)),
])
def test_thermostats_reach_the_target_temperature(name, kw):
    from trilattice import make_thermostat

    cfg = triangular_lattice(8, 5, 3.2010)
    sim = Simulation(cfg, POT, timestep=0.002,
                     thermostat=make_thermostat(name, 500.0, **kw), seed=9)
    sim.set_temperature(50.0)
    sim.run(15000, log_every=0)
    sim.clear_log()
    log = sim.run(30000, log_every=20)
    assert np.mean(log.temperature) == pytest.approx(500.0, rel=0.05)


@pytest.mark.parametrize("name,kw", [("bussi", dict(tau=0.1)), ("nose-hoover", dict(tau=0.1))])
def test_canonical_thermostats_have_the_right_temperature_variance(name, kw):
    """var(T)/T^2 = 2/N_f in the canonical ensemble.  Berendsen fails this."""
    from trilattice import make_thermostat

    cfg = triangular_lattice(8, 5, 3.2010)
    sim = Simulation(cfg, POT, timestep=0.002,
                     thermostat=make_thermostat(name, 600.0, **kw), seed=13)
    sim.set_temperature(600.0)
    sim.run(20000, log_every=0)
    sim.clear_log()
    log = sim.run(120000, log_every=10)
    t = np.array(log.temperature)
    assert np.var(t) / t.mean() ** 2 == pytest.approx(2.0 / sim.dof, rel=0.25)


def test_berendsen_underestimates_the_fluctuations():
    from trilattice import make_thermostat

    cfg = triangular_lattice(8, 5, 3.2010)
    sim = Simulation(cfg, POT, timestep=0.002,
                     thermostat=make_thermostat("berendsen", 600.0, tau=0.1), seed=13)
    sim.set_temperature(600.0)
    sim.run(20000, log_every=0)
    sim.clear_log()
    t = np.array(sim.run(120000, log_every=10).temperature)
    assert np.var(t) / t.mean() ** 2 < 0.7 * (2.0 / sim.dof)


def test_minimize_finds_the_perfect_lattice():
    cfg = triangular_lattice(8, 5, 3.2010)
    rng = np.random.default_rng(0)
    cfg.positions += 0.08 * rng.standard_normal(cfg.positions.shape)
    sim = Simulation(cfg, POT, timestep=0.002, seed=1)
    e_before = sim.e_pot
    info = sim.minimize(max_steps=3000, f_tol=1e-8)
    assert sim.e_pot < e_before
    assert info["f_max"] < 1e-6
