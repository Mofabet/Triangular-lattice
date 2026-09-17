"""The viewer is exercised headlessly: no window, but every code path runs."""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

import trilattice as tl
from trilattice.animate import COLOR_MODES, LiveViewer, ViewerConfig, is_headless

A = 3.2010
POT = tl.LennardJones(0.27, 2.88, 8.0, mode="shifted-force")


def _viewer(**kw):
    cfg = tl.triangular_lattice(8, 5, A)
    sim = tl.Simulation(cfg, POT, timestep=0.002,
                        thermostat=tl.Bussi(temperature=400.0, tau=0.2), seed=4)
    sim.set_temperature(400.0)
    return LiveViewer(sim, A, ViewerConfig(steps_per_frame=5, **kw))


def test_headless_detection():
    assert is_headless()


@pytest.mark.parametrize("mode", COLOR_MODES)
def test_every_colour_mode_produces_finite_per_atom_values(mode):
    v = _viewer(color_mode=mode)
    v.update()
    colors, _, _ = v._colors()
    assert colors.shape == (v.sim.n_particles,)
    assert np.all(np.isfinite(colors))


def test_update_advances_time_and_fills_the_traces():
    v = _viewer()
    for _ in range(5):
        v.update()
    assert v.sim.step_count == 25
    assert len(v.psi6_hist) == 5
    assert v.sim.elapsed == pytest.approx(25 * 0.002)
    assert 0.0 <= v.psi6_hist[-1] <= 1.0


def test_pause_freezes_the_dynamics():
    v = _viewer()
    v.update()
    n = v.sim.step_count
    v.paused = True
    v.update()
    assert v.sim.step_count == n


def test_history_is_bounded_by_the_config():
    v = _viewer(history=4)
    for _ in range(12):
        v.update()
    assert len(v.t_hist) == 4 and len(v.psi6_hist) == 4


def test_temperature_keys_drive_the_thermostat():
    v = _viewer()
    t0 = v.target_temperature

    class Key:
        key = "up"

    v._on_key(Key())
    assert v.target_temperature == pytest.approx(t0 + 50.0)
    Key.key = "shift+down"
    v._on_key(Key())
    assert v.target_temperature == pytest.approx(t0 - 200.0)
    assert v.sim.thermostat.temperature == pytest.approx(t0 - 200.0)


def test_colour_cycling_key():
    v = _viewer(color_mode="psi6")

    class Key:
        key = "c"

    v._on_key(Key())
    assert v.cfg.color_mode == COLOR_MODES[1]


def test_reset_restores_the_initial_configuration():
    v = _viewer()
    start = v.sim.positions.copy()
    for _ in range(10):
        v.update()
    assert not np.allclose(v.sim.positions, start)

    class Key:
        key = "r"

    v._on_key(Key())
    assert np.allclose(v.sim.positions, start)
    assert v.sim.elapsed == 0.0


def test_save_writes_a_gif(tmp_path):
    v = _viewer()
    out = tmp_path / "a.gif"
    v.save(str(out), frames=3, fps=4, dpi=40)
    assert out.exists() and out.stat().st_size > 0
