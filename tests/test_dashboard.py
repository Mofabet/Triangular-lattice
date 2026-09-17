"""The dashboard is driven headlessly: callbacks are invoked directly."""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

import matplotlib.pyplot as plt

import trilattice as tl
from trilattice.dashboard import COLOR_MODES, PANELS, RATES, Dashboard, DashboardConfig, hollow_sites

A = 3.2010
POT = tl.LennardJones(0.27, 2.88, 8.0, mode="shifted-force")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _dash(**kw):
    cfg = tl.triangular_lattice(9, 5, A)
    sim = tl.Simulation(cfg, POT, timestep=0.002,
                        thermostat=tl.Bussi(temperature=500.0, tau=0.2), seed=6)
    sim.set_temperature(500.0)
    kw.setdefault("rate_index", RATES.index(5.0))
    kw.setdefault("heavy_every", 2)
    return Dashboard(sim, A, DashboardConfig(**kw))


# --- rate control --------------------------------------------------------- #
def test_rate_zero_freezes_the_dynamics():
    d = _dash(rate_index=0)
    for _ in range(6):
        d.update()
    assert d.sim.step_count == 0


def test_slow_motion_advances_less_than_one_step_per_frame():
    """The point of slow motion: fewer steps per frame, never a smaller dt."""
    d = _dash(rate_index=RATES.index(0.25))
    dt_before = d.sim.dt
    for _ in range(8):
        d.update()
    assert d.sim.step_count == 2           # 8 frames * 1/4
    assert d.sim.dt == dt_before


def test_rate_slider_sets_the_rate():
    d = _dash()
    d.s_rate.set_val(RATES.index(40.0))
    assert d.rate == 40.0
    assert "40" in d._rate_label()


# --- temperature ---------------------------------------------------------- #
def test_temperature_slider_reaches_the_thermostat():
    d = _dash()
    d.s_temp.set_val(1234.0)
    assert d.sim.thermostat.temperature == pytest.approx(1234.0)


def test_thermostat_toggle_switches_to_nve_and_back():
    d = _dash()
    assert d.thermostat_on
    d.checks.set_active(3)
    assert not d.thermostat_on
    t = d.target_temperature          # still readable while off
    assert t > 0
    d.checks.set_active(3)
    assert d.thermostat_on
    assert d.sim.thermostat.temperature == pytest.approx(t)


# --- time travel ---------------------------------------------------------- #
def test_scrubbing_shows_a_past_frame_without_advancing():
    d = _dash()
    for _ in range(12):
        d.update()
    n = d.sim.step_count
    d.s_scrub.set_val(0.3)
    d.update()
    assert d.scrubbing
    assert d.sim.step_count == n               # frozen while scrubbing
    shown = d.scatter.get_offsets()
    assert np.allclose(shown, d.frames[d.scrub_index].positions)


def test_resume_here_rewinds_the_simulation_and_drops_the_future():
    d = _dash()
    for _ in range(16):
        d.update()
    d.s_scrub.set_val(0.25)
    target = d.frames[d.scrub_index]
    t_then, pos_then = target.elapsed, target.positions.copy()
    d._resume_here()
    assert d.sim.elapsed == pytest.approx(t_then)
    assert np.allclose(d.sim.positions, pos_then)
    assert not d.scrubbing
    # forces and neighbour list were rebuilt, so it keeps integrating
    for _ in range(5):
        d.update()
    assert np.isfinite(d.sim.e_pot)


def test_history_is_bounded():
    d = _dash(history=5)
    for _ in range(20):
        d.update()
    assert len(d.frames) == 5


# --- defects -------------------------------------------------------------- #
def test_add_vacancy_shrinks_the_system_consistently():
    d = _dash()
    d.update()
    n = d.sim.n_particles
    d._add_vacancy()
    assert d.sim.n_particles == n - 1
    for arr in (d.sim.velocities, d.sim.unwrapped, d.sim.masses, d.sim.types, d.sim.forces):
        assert arr.shape[0] == n - 1
    d.update()
    assert np.isfinite(d.sim.e_pot)


def test_add_interstitial_grows_it_and_stays_finite():
    d = _dash()
    d.update()
    n = d.sim.n_particles
    d._add_interstitial()
    assert d.sim.n_particles == n + 1
    assert np.max(np.abs(d.sim.forces)) < 1e3      # the capped relaxation did its job
    for _ in range(5):
        d.update()
    assert np.isfinite(d.sim.total_energy)


def test_hollow_sites_are_further_from_atoms_than_random_points():
    cfg = tl.triangular_lattice(9, 5, A)
    sites = hollow_sites(cfg.positions, cfg.box, A)
    d = cfg.box.minimum_image(sites[:, None, :] - cfg.positions[None, :, :])
    dmin = np.min(np.hypot(d[:, :, 0], d[:, :, 1]), axis=1)
    assert np.allclose(dmin, A / np.sqrt(3), atol=1e-9)


# --- panels and overlays -------------------------------------------------- #
@pytest.mark.parametrize("panel", PANELS)
def test_every_analysis_panel_renders(panel):
    d = _dash(panel=panel)
    for _ in range(8):
        d.update()
    d._draw_panel()
    assert d.ax_panel.get_xlabel() != "" or panel == "S(k)"


@pytest.mark.parametrize("mode", COLOR_MODES)
def test_every_colour_mode_gives_finite_values(mode):
    d = _dash(color_mode=mode)
    for _ in range(6):
        d.update()
    colors, _, _ = d._colors(d.sim.positions)
    assert colors.shape[0] == d.sim.n_particles
    assert np.all(np.isfinite(colors))


def test_bond_overlay_skips_bonds_that_cross_the_boundary():
    d = _dash()
    d.update()
    segs = d._bond_segments(d.sim.positions)
    lengths = np.linalg.norm(segs[:, 0] - segs[:, 1], axis=1)
    assert np.all(lengths < 1.35 * A + 1e-9)


def test_reset_restores_the_initial_state():
    d = _dash()
    start = d.sim.positions.copy()
    for _ in range(12):
        d.update()
    d._reset()
    assert np.allclose(d.sim.positions, start)
    assert d.sim.elapsed == 0.0
    assert len(d.frames) == 0


def test_quench_lowers_the_potential_energy():
    d = _dash()
    for _ in range(10):
        d.update()
    before = d.sim.e_pot
    d._quench()
    assert d.sim.e_pot <= before + 1e-9


def test_recording_writes_a_file(tmp_path):
    d = _dash()
    out = tmp_path / "d.gif"
    d.save(str(out), frames=3, fps=4, dpi=40)
    assert out.exists() and out.stat().st_size > 0


def test_record_script_can_drive_the_controls(tmp_path):
    d = _dash()
    seen = []
    d.save(str(tmp_path / "s.gif"), frames=4, fps=4, dpi=40,
           script=lambda dash, k: (seen.append(k), dash.s_temp.set_val(300.0 + 100 * k)))
    # FuncAnimation renders frame 0 twice (the initial draw, then the sequence)
    assert sorted(set(seen)) == [0, 1, 2, 3]
    assert d.sim.thermostat.temperature == pytest.approx(600.0)
