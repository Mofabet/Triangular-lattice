"""The dashboard is driven headlessly: callbacks are invoked directly."""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

import matplotlib.pyplot as plt

import trilattice as tl
from trilattice import observables as obs
from trilattice.dashboard import (
    COLOR_MODES,
    PANELS,
    RATES,
    Dashboard,
    DashboardConfig,
    hollow_sites,
)

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
    d.physics.set_active(d.PHYSICS.index("thermostat"))
    assert not d.thermostat_on
    t = d.target_temperature          # still readable while off
    assert t > 0
    d.physics.set_active(d.PHYSICS.index("thermostat"))
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


# --- lattice switching ----------------------------------------------------- #
def _dash_with_settings(lattice="triangular"):
    from trilattice.config import Settings
    from trilattice.dashboard import DashboardConfig

    s = Settings()
    s.system.lattice = lattice
    s.system.a = 3.3567
    s.system.nx, s.system.ny = 8, 5
    sim = tl.Simulation(s.build_configuration(), s.build_potential(),
                        timestep=s.run.timestep,
                        thermostat=s.build_thermostat(400.0), seed=6)
    sim.set_temperature(400.0)
    return Dashboard(sim, s.system.a, DashboardConfig(rate_index=RATES.index(5.0),
                                                      heavy_every=2),
                     lattice=lattice, settings=s)


@pytest.mark.parametrize("target", ["square", "honeycomb", "kagome"])
def test_switching_lattice_rebuilds_a_consistent_system(target):
    from trilattice.lattices import LATTICE_NAMES, lattice_spec

    d = _dash_with_settings()
    for _ in range(4):
        d.update()
    d.r_lattice.set_active(LATTICE_NAMES.index(target))
    spec = lattice_spec(target)
    assert d.spec.name == target
    n = d.sim.n_particles
    for arr in (d.sim.velocities, d.sim.unwrapped, d.sim.masses, d.sim.types, d.sim.forces):
        assert arr.shape[0] == n
    # the freshly built structure is ideal: right coordination, psi = 1
    assert d.sim.neighbors.r_list <= d.sim.box.max_cutoff
    value = obs.global_psi_n(d.sim.positions, d.sim.box, d.shell_cutoff, spec.psi_order)
    assert value == pytest.approx(1.0, abs=1e-6)
    for _ in range(4):
        d.update()
    assert np.isfinite(d.sim.e_pot)


def test_switching_keeps_the_particle_count_in_range():
    from trilattice.lattices import LATTICE_NAMES

    d = _dash_with_settings()
    n0 = d.sim.n_particles
    for name in ("square", "honeycomb", "kagome", "triangular"):
        d.r_lattice.set_active(LATTICE_NAMES.index(name))
        assert 0.7 * n0 <= d.sim.n_particles <= 1.4 * n0


def test_switching_to_a_metastable_lattice_lowers_the_target_temperature():
    """Otherwise the structure is gone before the first frame is drawn."""
    from trilattice.lattices import LATTICE_NAMES, lattice_spec

    d = _dash_with_settings()
    d.s_temp.set_val(1500.0)
    d.r_lattice.set_active(LATTICE_NAMES.index("kagome"))
    assert d.target_temperature <= lattice_spec("kagome").metastable_below
    assert d.sim.thermostat.temperature == pytest.approx(d.target_temperature)


def test_switching_updates_the_axes_to_the_new_box():
    from trilattice.lattices import LATTICE_NAMES

    d = _dash_with_settings()
    d.r_lattice.set_active(LATTICE_NAMES.index("honeycomb"))
    box = d.sim.box
    assert d.ax_cfg.get_xlim() == pytest.approx((-0.04 * box.lx, 1.04 * box.lx))
    assert d.ax_cfg.get_ylim() == pytest.approx((-0.04 * box.ly, 1.04 * box.ly))
    assert d.box_patch.get_width() == pytest.approx(box.lx)


def test_switching_without_settings_is_refused_cleanly():
    d = _dash()          # built without a Settings object
    from trilattice.lattices import LATTICE_NAMES

    before = d.spec.name
    d.r_lattice.set_active(LATTICE_NAMES.index("square"))
    assert d.spec.name == before
    assert d.sim.n_particles > 0


# --- the panel geometry regression ----------------------------------------- #
def test_panel_keeps_its_rectangle_through_every_switch():
    """S(k) sets aspect=equal; Axes.clear() does not undo it, and without an
    explicit reset every later panel gets squeezed to a sliver."""
    d = _dash()
    rect = d._panel_rect
    for k in (0, 1, 0, 2, 1, 3, 4, 2, 5, 1, 0):
        d.r_panel.set_active(k)
        for _ in range(3):
            d.update()
        d.fig.canvas.draw()
        p = d.ax_panel.get_position()
        assert (p.x0, p.y0, p.width, p.height) == pytest.approx(tuple(rect), abs=1e-9)


# --- layout ---------------------------------------------------------------- #
def test_nothing_overlaps_in_the_default_window():
    d = _dash_with_settings()
    for _ in range(6):
        d.update()
    assert d.layout_conflicts() == []


def test_nothing_overlaps_at_any_rate_setting():
    """The rate slider draws its value outside its own axes; the longest label
    used to run straight into the lattice selector."""
    d = _dash_with_settings()
    d.update()
    for k in range(len(RATES)):
        d.s_rate.set_val(k)
        d.update()
        assert d.layout_conflicts() == [], f"rate {RATES[k]}"


@pytest.mark.parametrize("panel", PANELS)
def test_nothing_overlaps_for_any_panel(panel):
    d = _dash_with_settings()
    d.r_panel.set_active(PANELS.index(panel))
    for _ in range(4):
        d.update()
    assert d.layout_conflicts() == []


@pytest.mark.parametrize("mode", COLOR_MODES)
def test_nothing_overlaps_for_any_colour_mode(mode):
    d = _dash_with_settings()
    d.r_color.set_active(COLOR_MODES.index(mode))
    for _ in range(4):
        d.update()
    assert d.layout_conflicts() == []


def test_nothing_overlaps_on_any_lattice():
    from trilattice.lattices import LATTICE_NAMES

    d = _dash_with_settings()
    for name in LATTICE_NAMES:
        d.r_lattice.set_active(LATTICE_NAMES.index(name))
        for _ in range(3):
            d.update()
        assert d.layout_conflicts() == [], name


def test_nothing_overlaps_at_the_extremes_of_the_sliders():
    d = _dash_with_settings()
    d.update()
    for value in (10.0, 4000.0):
        d.s_temp.set_val(value)
        d.update()
        assert d.layout_conflicts() == []
    d.s_scrub.set_val(0.3)
    d.update()
    assert d.layout_conflicts() == []


def test_the_layout_checker_actually_detects_an_overlap():
    """Guards the guard: a checker that never fires is worthless."""
    d = _dash_with_settings()
    d.update()
    assert d.layout_conflicts() == []
    d.r_lattice.ax.set_position(d.s_rate.ax.get_position())
    assert d.layout_conflicts() != []


def test_every_declared_rectangle_is_inside_the_figure():
    from trilattice.dashboard import LAYOUT

    for name, rect in LAYOUT.items():
        if not isinstance(rect, tuple) or len(rect) != 4:
            continue
        x, y, w, h = rect
        assert 0.0 <= x and 0.0 <= y and x + w <= 1.0 and y + h <= 1.0, name


# --- the physics switches --------------------------------------------------- #
def test_directional_switch_adds_and_removes_the_angular_term():
    """It lives in its own group: it changes the model, not the drawing."""
    d = _dash_with_settings()
    assert "directional" in d.PHYSICS and "directional" not in d.OVERLAYS
    assert d.sim.three_body is None
    d.physics.set_active(d.PHYSICS.index("directional"))
    assert d.sim.three_body is not None
    d.update()
    assert np.isfinite(d.sim.e_pot)
    d.physics.set_active(d.PHYSICS.index("directional"))
    assert d.sim.three_body is None


def test_the_angular_term_survives_a_lattice_switch_and_matches_the_new_order():
    from trilattice.lattices import LATTICE_NAMES, lattice_spec

    d = _dash_with_settings()
    d.physics.set_active(d.PHYSICS.index("directional"))
    for name in ("square", "honeycomb", "kagome"):
        d.r_lattice.set_active(LATTICE_NAMES.index(name))
        assert d.sim.three_body is not None
        assert d.sim.three_body.order == lattice_spec(name).angular_order


def test_directional_bonding_skips_the_cold_start_clamp():
    """With the angular term on, a metastable lattice no longer needs it."""
    from trilattice.lattices import LATTICE_NAMES

    d = _dash_with_settings()
    d.s_temp.set_val(900.0)
    d.physics.set_active(d.PHYSICS.index("directional"))
    d.r_lattice.set_active(LATTICE_NAMES.index("honeycomb"))
    assert d.target_temperature == pytest.approx(900.0)


def test_sk_window_is_square_and_reaches_the_first_bragg_ring():
    """The window used to be 2*pi*n_max/L, so it depended on the box rather than
    on the structure and cut the triangular and square peaks out of the frame."""
    from trilattice.lattices import LATTICE_NAMES, lattice_spec

    d = _dash_with_settings()
    d.r_panel.set_active(PANELS.index("S(k)"))
    for name in LATTICE_NAMES:
        d.r_lattice.set_active(LATTICE_NAMES.index(name))
        for _ in range(4):
            d.update()
        kx, ky, _ = d._heavy_sk
        spec = lattice_spec(name)
        assert kx.max() >= spec.first_bragg(d.a)
        assert ky.max() >= spec.first_bragg(d.a)
        assert kx.max() == pytest.approx(ky.max(), rel=0.15)   # square in k
