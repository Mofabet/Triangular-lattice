"""An interactive dashboard around a running simulation.

Where :class:`~trilattice.animate.LiveViewer` is a window that shows the
dynamics, this is an instrument you can drive.  It adds three things.

**Control.**  Sliders for the target temperature and for the simulation rate --
the rate goes below one step per frame, so the motion can be watched in genuine
slow motion, one integration step at a time if wanted.  Buttons quench, reheat,
reset, toggle the thermostat off into NVE, and inject vacancies or interstitials
into the running crystal.

**Time.**  Every frame is pushed into a ring buffer of full dynamical states.
The scrub slider plays back through it, and ``Resume here`` restores that state
and continues from it, discarding the future.  Rewinding a molecular dynamics
trajectory and taking a different branch is otherwise surprisingly fiddly: the
neighbour list and the force array have to be rebuilt in step with the
positions, which is what :meth:`~trilattice.simulation.Simulation.set_state`
handles.

**Measurement.**  The analysis panel cycles through live versions of the figures
this package produces offline -- g(r), S(k), the mean squared displacement, the
speed distribution against the exact 2-D Maxwell-Boltzmann curve, the Voronoi
coordination histogram, and the (T, |psi_6|) trajectory that traces out the
melting curve as the temperature slider is moved.

Expensive quantities are throttled: the Delaunay triangulation and the structure
factor are recomputed every few frames, not every frame, which keeps the whole
thing interactive at a few thousand integration steps per second.
"""

from __future__ import annotations

import time
from collections import deque
from dataclasses import dataclass, field

import matplotlib as mpl
import numpy as np

from . import observables as obs
from .lattice import SQRT3, Box
from .lattices import LATTICE_NAMES, emptiest_points, lattice_spec
from .plotting import COLORS, use_style
from .simulation import Simulation
from .thermostats import NoThermostat
from .units import KB, MVV2E

COLOR_MODES = ("psi6", "coord", "speed", "energy", "displacement")
PANELS = ("g(r)", "S(k)", "MSD", "speeds", "coordination", "T-psi6")

#: Simulation rate in integration steps per rendered frame.  Values below one
#: advance the system only every few frames, which is what slow motion means
#: here -- the timestep itself must never be changed for that purpose.
RATES = (0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0, 320.0)

#: Every rectangle in the window, in figure coordinates, in one place.
#:
#: Widgets are placed with ``add_axes`` rather than a grid, so nothing stops two
#: of them from being put on top of each other -- and a slider's value text is
#: drawn *outside* its own axes, to the right, which is how the lattice selector
#: ended up underneath "1 step / 8 frames".  Keeping the geometry in a single
#: table makes the arrangement reviewable, and :meth:`Dashboard.layout_conflicts`
#: checks it against the *rendered* extents of every control and label, so a
#: future move that overlaps something fails a test instead of shipping.
LAYOUT = {
    # main axes
    "ax_cfg":      (0.035, 0.300, 0.400, 0.620),
    "ax_panel":    (0.530, 0.665, 0.290, 0.255),
    "ax_t":        (0.530, 0.435, 0.290, 0.175),
    "ax_e":        (0.530, 0.205, 0.290, 0.175),
    # sliders (left of the bottom band); their labels sit left, values right
    "s_temp":      (0.065, 0.175, 0.270, 0.022),
    "s_rate":      (0.065, 0.130, 0.270, 0.022),
    "s_scrub":     (0.065, 0.085, 0.270, 0.022),
    # buttons, two rows of five
    "buttons":     (0.065, 0.038, 0.058, 0.024),
    "button_dx":   0.063,
    "button_dy":   0.030,
    # lattice selector: the empty band under the right-hand plots
    "r_lattice":   (0.545, 0.010, 0.115, 0.130),
    "t_lattice":   (0.545, 0.150),
    # right margin
    "r_color":     (0.850, 0.685, 0.140, 0.235),
    "t_color":     (0.850, 0.925),
    "r_panel":     (0.850, 0.395, 0.140, 0.255),
    "t_panel":     (0.850, 0.655),
    "checks":      (0.850, 0.215, 0.140, 0.115),
    "t_checks":    (0.850, 0.340),
    "physics":     (0.850, 0.055, 0.140, 0.090),
    "t_physics":   (0.850, 0.155),
}

_CMAP = {
    "psi6": ("viridis", r"$|\psi_6|$", 0.0, 1.0),
    "coord": ("coolwarm", "coordination", 4.0, 8.0),
    "speed": ("inferno", "speed (A/ps)", None, None),
    "energy": ("magma", "$E_\\mathrm{pot}$ (eV/atom)", None, None),
    "displacement": ("cividis", "displacement (A)", 0.0, None),
}


def hollow_sites(positions: np.ndarray, box: Box, a: float) -> np.ndarray:
    """The two three-fold hollows per atom of a triangular lattice.

    Same construction as :func:`trilattice.lattice.interstitial_sites`, taken on
    raw arrays so that it also works on a configuration that is no longer a
    perfect lattice.
    """
    offsets = np.array([[0.5 * a, a * SQRT3 / 6.0], [0.0, a * SQRT3 / 3.0]])
    return box.wrap((positions[:, None, :] + offsets[None, :, :]).reshape(-1, 2))


@dataclass
class Frame:
    """One stored state: enough to restart the dynamics exactly."""

    positions: np.ndarray
    velocities: np.ndarray
    unwrapped: np.ndarray
    elapsed: float
    temperature: float
    e_pot: float
    e_tot: float
    psi6: float
    defects: float
    target: float


@dataclass
class DashboardConfig:
    rate_index: int = 7               #: index into RATES; 7 -> 10 steps/frame
    history: int = 500                #: frames kept for scrubbing
    color_mode: str = "psi6"
    panel: str = "g(r)"
    psi6_cutoff_factor: float = 1.35
    point_size: float = 26.0
    interval_ms: int = 20
    heavy_every: int = 4              #: frames between Delaunay / S(k) updates
    trail_length: int = 40
    n_trails: int = 14
    t_min: float = 10.0
    t_max: float = 4000.0
    figsize: tuple[float, float] = (15.0, 8.6)


class Dashboard:
    """Interactive control panel for a :class:`Simulation`."""

    def __init__(self, sim: Simulation, a: float, config: DashboardConfig | None = None,
                 lattice: str = "triangular", settings=None):
        import matplotlib.pyplot as plt

        self.sim = sim
        self.a = a
        self.cfg = config or DashboardConfig()
        #: symmetry metadata -- decides the order parameter, the first-shell
        #: cutoff and what counts as a defect.  None of those transfer between
        #: lattices, so they are looked up rather than hard-coded.
        self.spec = lattice_spec(lattice)
        #: a Settings object enables switching lattice while the run is live
        self.settings = settings
        #: lambda used when the "directional" switch is turned on
        self.three_body_strength = (
            settings.potential.three_body if settings is not None else 0.0
        ) or 0.5
        self.rate_index = self.cfg.rate_index
        self.paused = False
        self.scrubbing = False
        self.scrub_index = -1
        self._step_accumulator = 0.0
        self._frame = 0
        self._sps = 0.0
        self._message = ""
        self._message_until = 0.0
        self._live_thermostat = sim.thermostat
        self._initial = Frame(
            sim.positions.copy(), sim.velocities.copy(), sim.unwrapped.copy(),
            0.0, sim.temperature, sim.e_pot, sim.total_energy, 1.0, 0.0,
            self.target_temperature,
        )

        self.frames: deque[Frame] = deque(maxlen=self.cfg.history)
        self.msd_origin = sim.unwrapped.copy()
        self.msd_origin_time = sim.elapsed
        self.rdf = obs.RDFAccumulator(sim.box, r_max=min(18.0, sim.box.max_cutoff), n_bins=220)
        self.trail_ids = np.linspace(0, sim.n_particles - 1,
                                     min(self.cfg.n_trails, sim.n_particles)).astype(int)
        self.trails: deque = deque(maxlen=self.cfg.trail_length)

        self._heavy_coord = np.full(sim.n_particles, 6.0)
        self._heavy_defects = 0.0
        self._heavy_sk = None

        use_style()
        self.fig = plt.figure(figsize=self.cfg.figsize)
        self._build_axes()
        self._build_widgets()
        self._draw_static()

    # ================================================================== #
    # properties
    # ================================================================== #
    @property
    def shell_cutoff(self) -> float:
        """First-neighbour-shell cutoff for this lattice, in angstrom."""
        return self.spec.cutoff(self.a)

    def _coordination(self) -> np.ndarray:
        """Delaunay on the triangular lattice, a distance cutoff otherwise.

        A square lattice has no unique Delaunay triangulation -- each plaquette
        can be cut along either diagonal -- so counting triangulation neighbours
        there produces defects that are not physical.
        """
        if self.spec.use_delaunay():
            return obs.coordination_by_delaunay(self.sim.positions, self.sim.box)
        return obs.coordination_by_cutoff(self.sim.positions, self.sim.box,
                                          self.shell_cutoff)

    @property
    def rate(self) -> float:
        return RATES[self.rate_index]

    @property
    def target_temperature(self) -> float:
        th = self.sim.thermostat
        if hasattr(th, "temperature"):
            return float(th.temperature)
        return float(getattr(self._live_thermostat, "temperature", self.sim.temperature))

    @target_temperature.setter
    def target_temperature(self, value: float) -> None:
        value = float(np.clip(value, self.cfg.t_min, self.cfg.t_max))
        for th in (self.sim.thermostat, self._live_thermostat):
            if hasattr(th, "temperature"):
                th.temperature = value

    @property
    def thermostat_on(self) -> bool:
        return not isinstance(self.sim.thermostat, NoThermostat)

    # ================================================================== #
    # layout
    # ================================================================== #
    def _build_axes(self) -> None:
        f = self.fig
        self.ax_cfg = f.add_axes(LAYOUT["ax_cfg"])
        self._panel_rect = list(LAYOUT["ax_panel"])
        self.ax_panel = f.add_axes(self._panel_rect)
        self.ax_t = f.add_axes(LAYOUT["ax_t"])
        self.ax_e = f.add_axes(LAYOUT["ax_e"])

    def _build_widgets(self) -> None:
        from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider

        f = self.fig
        self.s_temp = Slider(
            f.add_axes(LAYOUT["s_temp"]), "target T (K)",
            self.cfg.t_min, self.cfg.t_max, valinit=self.target_temperature, valfmt="%.0f",
            color=COLORS[1],
        )
        self.s_rate = Slider(
            f.add_axes(LAYOUT["s_rate"]), "rate", 0, len(RATES) - 1,
            valinit=self.rate_index, valstep=1, color=COLORS[0],
        )
        self.s_scrub = Slider(
            f.add_axes(LAYOUT["s_scrub"]), "time", 0.0, 1.0,
            valinit=1.0, color=COLORS[2],
        )
        self.s_rate.valtext.set_text(self._rate_label())
        self.s_scrub.valtext.set_text("live")
        self.s_temp.on_changed(self._on_temp)
        self.s_rate.on_changed(self._on_rate)
        self.s_scrub.on_changed(self._on_scrub)

        specs = [
            ("run", self._toggle_pause), ("live", self._go_live),
            ("resume here", self._resume_here), ("reset", self._reset),
            ("quench", self._quench), ("reheat", self._reheat),
            ("+vacancy", self._add_vacancy), ("+interstitial", self._add_interstitial),
            ("new origin", self._new_origin), ("clear stats", self._clear_stats),
        ]
        self.buttons = {}
        for k, (label, cb) in enumerate(specs):
            row, col = divmod(k, 5)
            bx, by, bw, bh = LAYOUT["buttons"]
            ax = f.add_axes([bx + col * LAYOUT["button_dx"],
                             by - row * LAYOUT["button_dy"], bw, bh])
            b = Button(ax, label, color="0.94", hovercolor="0.86")
            b.label.set_fontsize(7)
            b.on_clicked(cb)
            self.buttons[label] = b

        self.r_color = RadioButtons(
            f.add_axes(LAYOUT["r_color"]),
            COLOR_MODES, active=COLOR_MODES.index(self.cfg.color_mode),
            label_props={"fontsize": [9.5]},
        )
        self.r_panel = RadioButtons(
            f.add_axes(LAYOUT["r_panel"]),
            PANELS, active=PANELS.index(self.cfg.panel),
            label_props={"fontsize": [9.5]},
        )
        for r in (self.r_color, self.r_panel):
            r.ax.set_facecolor("none")
            r.ax.set_frame_on(False)
        self.r_color.on_clicked(self._on_color)
        self.r_panel.on_clicked(self._on_panel)

        self.r_lattice = RadioButtons(
            f.add_axes(LAYOUT["r_lattice"]), LATTICE_NAMES,
            active=LATTICE_NAMES.index(self.spec.name),
            label_props={"fontsize": [9.5]},
        )
        self.r_lattice.ax.set_facecolor("none")
        self.r_lattice.ax.set_frame_on(False)
        self.r_lattice.on_clicked(self._on_lattice)
        self.fig.text(*LAYOUT["t_lattice"], "lattice", fontsize=8, weight="bold")

        self.checks = CheckButtons(
            f.add_axes(LAYOUT["checks"]),
            list(self.OVERLAYS),
            [False, False, True],
            label_props={"fontsize": [9.5]},
        )
        self.checks.ax.set_facecolor("none")
        self.checks.ax.set_frame_on(False)
        self.checks.on_clicked(self._on_check)

        # Physics switches, kept apart from the drawing overlays: these change
        # what is being simulated, not what is being drawn.  Mixing the two was
        # the reason "directional" read as a display option.
        self.physics = CheckButtons(
            f.add_axes(LAYOUT["physics"]), list(self.PHYSICS),
            [True, self.sim.three_body is not None],
            label_props={"fontsize": [9.5]},
        )
        self.physics.ax.set_facecolor("none")
        self.physics.ax.set_frame_on(False)
        self.physics.on_clicked(self._on_physics)
        self.fig.text(*LAYOUT["t_physics"], "physics", fontsize=8, weight="bold")

        self.fig.text(*LAYOUT["t_color"], "colour by", fontsize=8, weight="bold")
        self.fig.text(*LAYOUT["t_panel"], "analysis panel", fontsize=8, weight="bold")
        self.fig.text(*LAYOUT["t_checks"], "overlays", fontsize=8, weight="bold")

    def _draw_static(self) -> None:
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        box = self.sim.box
        cmap, label, lo, hi = _CMAP[self.cfg.color_mode]
        self.bonds = LineCollection([], colors="0.72", linewidths=0.6, zorder=1)
        self.ax_cfg.add_collection(self.bonds)
        self.trail_lines = LineCollection([], colors=COLORS[1], linewidths=0.9,
                                          alpha=0.8, zorder=2)
        self.ax_cfg.add_collection(self.trail_lines)
        self.scatter = self.ax_cfg.scatter(
            self.sim.positions[:, 0], self.sim.positions[:, 1],
            c=np.zeros(self.sim.n_particles), cmap=cmap, vmin=lo, vmax=hi,
            s=self.cfg.point_size, linewidths=0.3, edgecolors="0.25", zorder=3,
        )
        self.box_patch = plt.Rectangle((0, 0), box.lx, box.ly, fill=False,
                                       ec="0.45", lw=0.8, ls="--", zorder=0)
        self.ax_cfg.add_patch(self.box_patch)
        self.ax_cfg.set_xlim(-0.04 * box.lx, 1.04 * box.lx)
        self.ax_cfg.set_ylim(-0.04 * box.ly, 1.04 * box.ly)
        self.ax_cfg.set_aspect("equal")
        self.ax_cfg.grid(False)
        self.ax_cfg.set_xlabel("x (A)")
        self.ax_cfg.set_ylabel("y (A)")
        self.cbar = self.fig.colorbar(self.scatter, ax=self.ax_cfg, fraction=0.040, pad=0.015)
        self.cbar.set_label(label, fontsize=8, rotation=270, labelpad=13)

        self.readout = self.ax_cfg.text(
            0.012, 0.985, "", transform=self.ax_cfg.transAxes, va="top", ha="left",
            fontsize=8, family="monospace", zorder=5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.75", alpha=0.9),
        )
        self.banner = self.fig.text(0.435, 0.945, "", fontsize=9, ha="right",
                                    color=COLORS[1], weight="bold")

        (self.line_t,) = self.ax_t.plot([], [], color=COLORS[1], lw=1.2, label="$T$")
        (self.line_ttarget,) = self.ax_t.plot([], [], color="0.55", lw=0.9, ls=":",
                                              label="target")
        self.ax_t.set_ylabel("$T$ (K)")
        self.ax_t.tick_params(labelbottom=False)
        self.ax_energy = self.ax_t.twinx()
        (self.line_etot,) = self.ax_energy.plot([], [], color=COLORS[2], lw=1.0)
        self.ax_energy.set_ylabel(r"$E_\mathrm{tot}$ (eV/atom)", color=COLORS[2], fontsize=8)
        self.ax_energy.tick_params(axis="y", labelcolor=COLORS[2], labelsize=7)
        self.ax_energy.grid(False)
        self.ax_t.legend(loc="lower right", fontsize=7, ncol=2, framealpha=0.85)

        (self.line_psi6,) = self.ax_e.plot([], [], color=COLORS[3], lw=1.2,
                                           label=r"$|\langle\psi_6\rangle|$")
        self.ax_e.set_ylabel(r"$|\langle\psi_n\rangle|$")
        self.ax_e.set_ylim(0, 1.02)
        self.ax_e.set_xlabel("time (ps)")
        self.ax_def = self.ax_e.twinx()
        (self.line_def,) = self.ax_def.plot([], [], color=COLORS[0], lw=1.0,
                                            label="defects")
        self.ax_def.set_ylabel("defects (%)", color=COLORS[0])
        self.ax_def.tick_params(axis="y", labelcolor=COLORS[0], labelsize=7)
        self.ax_def.grid(False)

        self._refresh_title()

    # ================================================================== #
    # callbacks
    # ================================================================== #
    def _notify(self, text: str, seconds: float = 2.5) -> None:
        self._message = text
        self._message_until = time.perf_counter() + seconds

    def _on_temp(self, value) -> None:
        self.target_temperature = float(value)

    def _on_rate(self, value) -> None:
        self.rate_index = int(value)
        self.s_rate.valtext.set_text(self._rate_label())

    def _rate_label(self) -> str:
        r = self.rate
        if r == 0:
            return "paused"
        if r < 1:
            return f"1/{int(round(1 / r))} step/frame"
        return f"{r:.0f} steps/frame"

    def _on_scrub(self, value) -> None:
        if not self.frames:
            return
        if value >= 0.999:
            self.scrubbing = False
            self.scrub_index = -1
            self.s_scrub.valtext.set_text("live")
            return
        self.scrubbing = True
        self.scrub_index = int(round(value * (len(self.frames) - 1)))
        fr = self.frames[self.scrub_index]
        self.s_scrub.valtext.set_text(f"{fr.elapsed:.2f} ps")

    def _go_live(self, _=None) -> None:
        self.scrubbing = False
        self.scrub_index = -1
        self.s_scrub.set_val(1.0)
        self.paused = False

    def _toggle_pause(self, _=None) -> None:
        self.paused = not self.paused
        self.buttons["run"].label.set_text("run" if self.paused else "pause")

    def _resume_here(self, _=None) -> None:
        if not self.scrubbing or not self.frames:
            self._notify("scrub first, then resume")
            return
        fr = self.frames[self.scrub_index]
        self.sim.set_state(fr.positions, fr.velocities, fr.elapsed, fr.unwrapped)
        self.target_temperature = fr.target
        self.s_temp.set_val(fr.target)
        for _ in range(len(self.frames) - self.scrub_index - 1):
            self.frames.pop()
        self._notify(f"resumed at t = {fr.elapsed:.2f} ps; future discarded")
        self._go_live()

    def _reset(self, _=None) -> None:
        fr = self._initial
        self.sim.set_state(fr.positions, fr.velocities, 0.0, fr.unwrapped)
        self.frames.clear()
        self.trails.clear()
        self.rdf = obs.RDFAccumulator(self.sim.box,
                                      r_max=min(18.0, self.sim.box.max_cutoff), n_bins=220)
        self._new_origin()
        self._notify("reset to the initial configuration")
        self._go_live()

    def _quench(self, _=None) -> None:
        info = self.sim.minimize(max_steps=400, f_tol=1e-7, max_move=0.05)
        self._notify(f"quenched: |F|max = {info['f_max']:.2e} eV/A")

    def _reheat(self, _=None) -> None:
        self.sim.set_temperature(self.target_temperature)
        self._notify(f"velocities redrawn at {self.target_temperature:.0f} K")

    def _add_vacancy(self, _=None) -> None:
        if self.sim.n_particles < 8:
            return
        idx = int(self.sim.rng.integers(self.sim.n_particles))
        self.sim.remove_particles([idx])
        self._resize_visuals()
        self._notify(f"removed atom {idx}: N = {self.sim.n_particles}")

    def _add_interstitial(self, _=None) -> None:
        """Insert at the emptiest point of the cell, then relax with capped FIRE.

        The insertion site is found by a grid search rather than from the ideal
        lattice geometry, so this works on a square lattice, a honeycomb or a
        half-melted configuration alike.  In a close-packed structure there is
        nowhere that does not overlap, so the relaxation is not optional; it
        zeroes the velocities, hence the re-thermalisation afterwards.
        """
        t_before = self.sim.temperature
        best = emptiest_points(self.sim.positions, self.sim.box, n=1)[0]
        self.sim.insert_particles(best, temperature=t_before)
        self.sim.minimize(max_steps=300, f_tol=1e-6, max_move=0.04)
        self.sim.set_temperature(max(t_before, 1.0))
        self._resize_visuals()
        self._notify(f"inserted interstitial: N = {self.sim.n_particles}, re-thermalised")

    def _new_origin(self, _=None) -> None:
        self.msd_origin = self.sim.unwrapped.copy()
        self.msd_origin_time = self.sim.elapsed
        self._initial_positions_for_displacement = self.sim.positions.copy()
        self._notify("MSD origin reset to now")

    def _clear_stats(self, _=None) -> None:
        self.rdf = obs.RDFAccumulator(self.sim.box,
                                      r_max=min(18.0, self.sim.box.max_cutoff), n_bins=220)
        self._notify("g(r) accumulator cleared")

    def _make_three_body(self):
        """The angular term matching the current lattice and spacing."""
        from .threebody import ThreeBodyAngular

        return ThreeBodyAngular.for_lattice(self.spec, self.a,
                                            strength=self.three_body_strength)

    def _on_lattice(self, label) -> None:
        if label == self.spec.name:
            return
        if self.settings is None:
            self._notify("lattice switching needs a Settings object")
            self.r_lattice.set_active(LATTICE_NAMES.index(self.spec.name))
            return
        self._switch_lattice(label)

    def _switch_lattice(self, name: str) -> None:
        """Rebuild the system on another lattice, live.

        The cell counts are chosen to keep the particle count roughly fixed
        (see :meth:`~trilattice.lattices.LatticeSpec.cells_for`) -- at fixed
        ``nx, ny`` a kagome cell would hold six times as many atoms as a square
        one and the dashboard would stall.  Everything sized by N or by the box
        is rebuilt: the neighbour list, the g(r) accumulator, the axes limits and
        the stored history, which cannot be carried across a change of geometry.
        """
        from .simulation import Simulation

        spec = lattice_spec(name)
        target_n = self.sim.n_particles
        target_t = self.target_temperature

        s = self.settings
        s.system.lattice = name
        s.system.nx, s.system.ny = spec.cells_for(target_n)
        config = s.build_configuration()
        potential = s.build_potential()

        if potential.cutoff + s.run.skin > config.box.max_cutoff:
            self._notify(f"{name}: box too small for the cutoff -- not switched")
            self.r_lattice.set_active(LATTICE_NAMES.index(self.spec.name))
            return

        directional = self.sim.three_body is not None
        new = Simulation(config, potential, timestep=s.run.timestep,
                         thermostat=self.sim.thermostat, skin=s.run.skin,
                         seed=s.run.seed)
        if directional:
            from .threebody import ThreeBodyAngular

            new.three_body = ThreeBodyAngular.for_lattice(
                spec, s.system.a, strength=self.three_body_strength)
        # A metastable structure vanishes within a few hundred femtoseconds at any
        # useful temperature, so start it somewhere it actually survives and let
        # the user raise the slider to trigger the collapse deliberately.
        if not spec.lj_ground_state and not directional and target_t > spec.metastable_below:
            target_t = max(10.0, spec.metastable_below)
            self.target_temperature = target_t
            self.s_temp.set_val(target_t)
        new.set_temperature(target_t)

        self.sim = new
        self.spec = spec
        self._live_thermostat = new.thermostat
        self._initial = Frame(
            new.positions.copy(), new.velocities.copy(), new.unwrapped.copy(),
            0.0, new.temperature, new.e_pot, new.total_energy, 1.0, 0.0, target_t,
        )
        self._rebuild_for_new_box()
        if directional:
            tail = "  -- held open by the angular term"
        elif spec.lj_ground_state:
            tail = ""
        elif spec.metastable_below <= 0.0:
            tail = "  -- unstable at any T, it will shear at once"
        else:
            tail = f"  -- metastable; started at {target_t:.0f} K, raise T to collapse it"
        self._notify(
            f"{name}: N = {new.n_particles}, z = {spec.coordination}, "
            f"psi_{spec.psi_order}" + tail, seconds=6.0,
        )

    def _rebuild_for_new_box(self) -> None:
        """Re-point every visual and accumulator at the current box and N."""
        box = self.sim.box
        self.frames.clear()
        self.trails.clear()
        self.rdf = obs.RDFAccumulator(box, r_max=min(18.0, box.max_cutoff), n_bins=220)
        self.msd_origin = self.sim.unwrapped.copy()
        self.msd_origin_time = self.sim.elapsed
        self._initial_positions_for_displacement = self.sim.positions.copy()
        self.trail_ids = np.linspace(
            0, self.sim.n_particles - 1,
            min(self.cfg.n_trails, self.sim.n_particles)).astype(int)
        self._heavy_coord = np.full(self.sim.n_particles, float(self.spec.coordination))
        self._heavy_defects = 0.0
        self.box_patch.set_bounds(0, 0, box.lx, box.ly)
        self.ax_cfg.set_xlim(-0.04 * box.lx, 1.04 * box.lx)
        self.ax_cfg.set_ylim(-0.04 * box.ly, 1.04 * box.ly)
        self.bonds.set_segments([])
        self.trail_lines.set_segments([])
        self._refresh_title()

    def _refresh_title(self) -> None:
        self.fig.suptitle(
            f"trilattice   {self.spec.name}   N = {self.sim.n_particles}   "
            f"{self.sim.box.lx:.1f} x {self.sim.box.ly:.1f} A   "
            f"dt = {self.sim.dt * 1000:.1f} fs", fontsize=10, x=0.235,
        )

    def _on_color(self, label) -> None:
        self.cfg.color_mode = label
        cmap, cl, lo, hi = _CMAP[label]
        self.scatter.set_cmap(cmap)
        self.cbar.set_label(cl, fontsize=8, rotation=270, labelpad=13)
        if lo is not None and hi is not None:
            self.scatter.set_clim(lo, hi)

    def _on_panel(self, label) -> None:
        self.cfg.panel = label
        self.ax_panel.clear()

    def _on_check(self, label) -> None:
        if label == "box":
            self.box_patch.set_visible(not self.box_patch.get_visible())
        return

    def _on_physics(self, label) -> None:
        if label == "directional":
            if self.sim.three_body is None:
                self.sim.three_body = self._make_three_body()
                self._notify("directional bonding on -- open lattices now hold", 4.0)
            else:
                self.sim.three_body = None
                self._notify("directional bonding off -- pair potential only", 4.0)
        elif label == "thermostat":
            if self.thermostat_on:
                self._live_thermostat = self.sim.thermostat
                self.sim.thermostat = NoThermostat()
                self._notify("thermostat off -- NVE")
            else:
                self.sim.thermostat = self._live_thermostat
                self._notify(f"thermostat on -- {self.sim.thermostat.name}")

    OVERLAYS = ("bonds", "trails", "box")
    PHYSICS = ("thermostat", "directional")

    def _checked(self, name: str) -> bool:
        if name in self.OVERLAYS:
            return bool(self.checks.get_status()[self.OVERLAYS.index(name)])
        return bool(self.physics.get_status()[self.PHYSICS.index(name)])

    def _resize_visuals(self) -> None:
        """Re-point everything that is sized by N after an insertion/removal."""
        n = self.sim.n_particles
        self.msd_origin = self.sim.unwrapped.copy()
        self.msd_origin_time = self.sim.elapsed
        self.trail_ids = np.linspace(0, n - 1, min(self.cfg.n_trails, n)).astype(int)
        self.trails.clear()
        self.frames.clear()
        self._heavy_coord = np.full(n, 6.0)
        self._refresh_title()

    # ================================================================== #
    # stepping
    # ================================================================== #
    def _advance(self) -> None:
        self._step_accumulator += self.rate
        n_steps = int(self._step_accumulator)
        self._step_accumulator -= n_steps
        if n_steps <= 0:
            return
        t0 = time.perf_counter()
        for _ in range(n_steps):
            self.sim.step()
        dt = time.perf_counter() - t0
        if dt > 0:
            self._sps = 0.75 * self._sps + 0.25 * n_steps / dt

    def _record(self, psi6_value: float) -> None:
        self.frames.append(Frame(
            self.sim.positions.copy(), self.sim.velocities.copy(),
            self.sim.unwrapped.copy(), self.sim.elapsed, self.sim.temperature,
            self.sim.e_pot / self.sim.n_particles,
            self.sim.total_energy / self.sim.n_particles,
            psi6_value, self._heavy_defects, self.target_temperature,
        ))

    # ================================================================== #
    # rendering
    # ================================================================== #
    def _display_positions(self) -> np.ndarray:
        if self.scrubbing and self.frames:
            return self.frames[self.scrub_index].positions
        return self.sim.positions

    def _colors(self, pos: np.ndarray):
        mode = self.cfg.color_mode
        box = self.sim.box
        if mode == "psi6":
            return (np.abs(obs.psi_n(pos, box, self.shell_cutoff, self.spec.psi_order)),
                    0.0, 1.0)
        if mode == "coord":
            if pos.shape[0] == self._heavy_coord.shape[0]:
                return self._heavy_coord, 4.0, 8.0
            return np.full(pos.shape[0], 6.0), 4.0, 8.0
        if mode == "speed":
            v = (self.frames[self.scrub_index].velocities if self.scrubbing and self.frames
                 else self.sim.velocities)
            s = np.linalg.norm(v, axis=1)
            return s, 0.0, max(float(np.percentile(s, 99)), 1e-6)
        if mode == "energy":
            e = self.sim.per_atom_energy
            if e.shape[0] != pos.shape[0]:
                e = np.zeros(pos.shape[0])
            return e, float(np.percentile(e, 1)), float(np.percentile(e, 99))
        ref = getattr(self, "_initial_positions_for_displacement", self._initial.positions)
        if ref.shape[0] != pos.shape[0]:
            ref = pos
        d = np.linalg.norm(self.sim.box.minimum_image(pos - ref), axis=1)
        return d, 0.0, max(float(np.percentile(d, 98)), 1e-6)

    def _bond_segments(self, pos: np.ndarray):
        i, j, d = obs.neighbour_pairs_within(self.sim.box, pos, self.shell_cutoff)
        raw = pos[i] - pos[j]
        keep = np.all(np.abs(raw - d) < 1e-9, axis=1)   # drop bonds crossing the boundary
        return np.stack([pos[i[keep]], pos[j[keep]]], axis=1)

    def _trail_segments(self):
        if len(self.trails) < 2:
            return []
        arr = np.asarray(self.trails)                       # (T, n_trail, 2)
        segs = []
        lengths = self.sim.box.lengths
        for k in range(arr.shape[1]):
            path = arr[:, k, :]
            jump = np.any(np.abs(np.diff(path, axis=0)) > 0.5 * lengths, axis=1)
            start = 0
            for b in np.flatnonzero(jump):
                if b + 1 - start > 1:
                    segs.append(path[start:b + 1])
                start = b + 1
            if len(path) - start > 1:
                segs.append(path[start:])
        return segs

    # ------------------------------------------------------------------ #
    def _draw_panel(self) -> None:
        ax = self.ax_panel
        name = self.cfg.panel
        ax.clear()
        # Axes.clear() does NOT reset the aspect ratio or the axis scales.  A
        # single visit to the S(k) panel would otherwise leave aspect=1 behind
        # for good, and matplotlib then squeezes every later panel's box to
        # satisfy it -- the histogram ends up a few thousandths of a figure
        # high.  Reassert all three, and the rectangle, before drawing.
        ax.set_aspect("auto", adjustable="box")
        ax.set_xscale("linear")
        ax.set_yscale("linear")
        ax.set_position(self._panel_rect)
        ax.set_autoscale_on(True)
        ax.grid(alpha=0.25, lw=0.5)
        sim = self.sim

        if name == "g(r)":
            r, g = self.rdf.result()
            ax.plot(r, g, color=COLORS[0], lw=1.1)
            ax.axhline(1.0, color="0.6", lw=0.7, ls=":")
            for shell in (1.0, self.spec.second_shell, 2.0):
                ax.axvline(shell * self.a, color="0.8", lw=0.6, ls=":")
            ax.set_xlabel("$r$ (A)")
            ax.set_ylabel("$g(r)$")
            ax.set_xlim(0, self.rdf.r_max)
            ax.set_title(f"g(r), {self.rdf.n_frames} frames accumulated", fontsize=8.5)

        elif name == "S(k)":
            if self._heavy_sk is not None:
                kx, ky, s = self._heavy_sk
                ax.pcolormesh(kx, ky, np.clip(s.T, 0.05, None),
                              norm=mpl.colors.LogNorm(vmin=0.05, vmax=max(2.0, sim.n_particles)),
                              cmap="magma", shading="auto")
                # adjustable="datalim" keeps the axes rectangle fixed and pads
                # the data range instead, so the panel never changes size
                ax.set_aspect("equal", adjustable="datalim")
            ax.set_xlabel(r"$k_x$ (A$^{-1}$)")
            ax.set_ylabel(r"$k_y$ (A$^{-1}$)")
            ax.set_title(r"$S(\mathbf{k})$", fontsize=8.5)
            ax.grid(False)

        elif name == "MSD":
            if len(self.frames) > 3:
                t = np.array([f.elapsed for f in self.frames]) - self.msd_origin_time
                m = t > 0
                if m.sum() > 2:
                    ref = self.msd_origin
                    msd = np.array([
                        np.mean(np.sum((f.unwrapped - ref) ** 2, axis=1))
                        if f.unwrapped.shape == ref.shape else np.nan
                        for f in self.frames
                    ])
                    ax.loglog(t[m], np.maximum(msd[m], 1e-6), color=COLORS[0], lw=1.2)
                    tt = t[m]
                    ax.loglog(tt, msd[m][-1] * tt / tt[-1], "k:", lw=0.8)
                    ax.axhline(self.a ** 2, color="0.6", lw=0.7, ls="--")
                    d_coef = obs.diffusion_coefficient(t[m], msd[m])
                    ax.set_title(f"MSD,  D = {d_coef:.4f} A$^2$/ps", fontsize=8.5)
            ax.set_xlabel("$t$ since origin (ps)")
            ax.set_ylabel(r"MSD (A$^2$)")

        elif name == "speeds":
            v = np.linalg.norm(sim.velocities, axis=1)
            ax.hist(v, bins=40, density=True, color=COLORS[0], alpha=0.75)
            t_now = max(sim.temperature, 1.0)
            m = float(sim.masses[0]) * MVV2E
            x = np.linspace(0, max(v.max(), 1e-6), 300)
            ax.plot(x, m * x / (KB * t_now) * np.exp(-m * x ** 2 / (2 * KB * t_now)),
                    color=COLORS[1], lw=1.3)
            ax.set_xlabel("speed (A/ps)")
            ax.set_ylabel("density")
            ax.set_title(f"2-D Maxwell-Boltzmann at {t_now:.0f} K", fontsize=8.5)

        elif name == "coordination":
            c = self._heavy_coord
            values, counts = np.unique(c.astype(int), return_counts=True)
            z0 = self.spec.coordination
            cols = ["0.7" if v == z0 else COLORS[1] if v < z0 else COLORS[0]
                    for v in values]
            ax.bar(values, 100 * counts / max(1, c.size), color=cols, width=0.7)
            ax.set_xlabel("coordination"
                          + (" (Delaunay)" if self.spec.use_delaunay() else " (cutoff)"))
            ax.set_ylabel("% of atoms")
            ax.set_xlim(max(0.5, z0 - 4.5), z0 + 4.5)
            ax.set_title(f"z != {z0}: {100 * self._heavy_defects:.1f} %", fontsize=8.5)

        else:  # T-psi6
            if self.frames:
                t = np.array([f.temperature for f in self.frames])
                p = np.array([f.psi6 for f in self.frames])
                ax.scatter(t, p, c=np.arange(len(t)), cmap="viridis", s=6, linewidths=0)
                ax.plot(t[-1], p[-1], "o", color=COLORS[1], ms=6)
            ax.set_xlabel("$T$ (K)")
            ax.set_ylabel(rf"$|\langle\psi_{{{self.spec.psi_order}}}\rangle|$")
            ax.set_ylim(0, 1.02)
            ax.set_title("order-disorder curve traced live", fontsize=8.5)

    # ------------------------------------------------------------------ #
    def update(self, _frame=None):
        run = not (self.paused or self.scrubbing)
        if run:
            self._advance()
        self._frame += 1

        sim = self.sim
        heavy = self._frame % self.cfg.heavy_every == 0
        if heavy and run:
            self._heavy_coord = self._coordination().astype(float)
            self._heavy_defects = float(
                np.mean(self._heavy_coord != self.spec.coordination))
            if self.cfg.panel == "S(k)":
                # the window is set by the lattice, not by the box: the same
                # n_max reaches a different |k| for every cell size, which is how
                # the triangular and square Bragg peaks ended up outside the frame
                self._heavy_sk = obs.structure_factor(
                    sim.positions, sim.box, k_max=self.spec.k_window(self.a))
            self.rdf.accumulate(sim.positions)

        p6 = obs.global_psi_n(sim.positions, sim.box, self.shell_cutoff,
                              self.spec.psi_order)
        if run:
            self._record(p6)
            self.trails.append(sim.positions[self.trail_ids].copy())
            if len(self.frames) > 1:
                self.s_scrub.eventson = False
                self.s_scrub.set_val(1.0)
                self.s_scrub.eventson = True
                self.s_scrub.valtext.set_text("live")

        pos = self._display_positions()
        colors, lo, hi = self._colors(pos)
        if colors.shape[0] != pos.shape[0]:
            colors = np.zeros(pos.shape[0])
        self.scatter.set_offsets(pos)
        self.scatter.set_array(colors)
        if lo is not None and hi is not None and hi > lo:
            self.scatter.set_clim(lo, hi)
        self.scatter.set_sizes(np.full(pos.shape[0], self.cfg.point_size))

        self.bonds.set_visible(self._checked("bonds"))
        if self._checked("bonds"):
            self.bonds.set_segments(list(self._bond_segments(pos)))
        self.trail_lines.set_visible(self._checked("trails"))
        if self._checked("trails"):
            self.trail_lines.set_segments(self._trail_segments())

        if self.frames:
            t = np.array([f.elapsed for f in self.frames])
            self.line_t.set_data(t, [f.temperature for f in self.frames])
            self.line_ttarget.set_data(t, [f.target for f in self.frames])
            self.line_psi6.set_data(t, [f.psi6 for f in self.frames])
            self.line_def.set_data(t, [100 * f.defects for f in self.frames])
            self.line_etot.set_data(t, [f.e_tot for f in self.frames])
            for ax in (self.ax_t, self.ax_e):
                ax.relim()
                ax.autoscale_view(scalex=True, scaley=(ax is self.ax_t))
                if t[-1] > t[0]:
                    ax.set_xlim(t[0], t[-1])
            for ax in (self.ax_def, self.ax_energy):
                ax.relim()
                ax.autoscale_view()
                if t[-1] > t[0]:
                    ax.set_xlim(t[0], t[-1])
            top = max(1.0, 100 * max(f.defects for f in self.frames) * 1.2)
            self.ax_def.set_ylim(0.0, top)

        if self._frame % 2 == 0 or not run:
            self._draw_panel()

        self._update_readout(p6)
        return (self.scatter, self.bonds, self.trail_lines, self.line_t,
                self.line_psi6, self.line_def, self.readout)

    def _update_readout(self, p6: float) -> None:
        sim = self.sim
        n = sim.n_particles
        if self.scrubbing and self.frames:
            fr = self.frames[self.scrub_index]
            head = (f"REWIND  frame {self.scrub_index + 1}/{len(self.frames)}\n"
                    f"t      {fr.elapsed:8.2f} ps\n"
                    f"T      {fr.temperature:8.1f} K\n"
                    f"E_pot  {fr.e_pot:8.4f} eV/at\n"
                    f"|psi6| {fr.psi6:8.3f}\n"
                    f"defects{100 * fr.defects:7.1f} %")
        else:
            head = (f"t      {sim.elapsed:8.2f} ps\n"
                    f"T      {sim.temperature:8.1f} K -> {self.target_temperature:.0f}\n"
                    f"E_pot  {sim.e_pot / n:8.4f} eV/at\n"
                    f"E_tot  {sim.total_energy / n:8.4f} eV/at\n"
                    f"|psi{self.spec.psi_order}| {p6:8.3f}\n"
                    f"defects{100 * self._heavy_defects:7.1f} % (z!={self.spec.coordination})\n"
                    f"P      {sim.pressure:8.4f} eV/A^2\n"
                    f"N      {n:8d}\n"
                    f"speed  {self._sps:8.0f} steps/s")
        state = []
        if self.paused:
            state.append("PAUSED")
        if not self.thermostat_on:
            state.append("NVE")
        if self.rate == 0:
            state.append("rate 0")
        if state:
            head += "\n" + "  ".join(state)
        self.readout.set_text(head)
        self.banner.set_text(
            self._message if time.perf_counter() < self._message_until else ""
        )

    # ================================================================== #
    # Layout verification
    # ================================================================== #
    def _layout_items(self) -> list[tuple[str, str, object]]:
        """``(group, label, artist)`` for everything the user has to be able to read.

        Items in the same group are allowed to overlap -- a radio button's text
        lives inside its own axes -- but nothing may cross a group boundary.
        """
        items: list[tuple[str, str, object]] = [
            ("config", "configuration axes", self.ax_cfg),
            ("colorbar", "colour bar", self.cbar.ax),
            ("panel", "analysis panel", self.ax_panel),
            ("traceT", "temperature trace", self.ax_t),
            ("tracePsi", "order/defect trace", self.ax_e),
        ]
        for group, name, sl in (("sTemp", "target T", self.s_temp),
                                ("sRate", "rate", self.s_rate),
                                ("sScrub", "time", self.s_scrub)):
            items += [(group, f"slider {name}", sl.ax),
                      (group, f"slider {name} label", sl.label),
                      (group, f"slider {name} value", sl.valtext)]
        for label, b in self.buttons.items():
            items.append((f"btn:{label}", f"button '{label}'", b.ax))
        for group, name, w in (("rColor", "colour selector", self.r_color),
                               ("rPanel", "panel selector", self.r_panel),
                               ("rLattice", "lattice selector", self.r_lattice),
                               ("checks", "overlay switches", self.checks),
                               ("physics", "physics switches", self.physics)):
            items.append((group, name, w.ax))
            for lab in w.labels:
                items.append((group, f"{name}: {lab.get_text()}", lab))
        for t in self.fig.texts:
            text = t.get_text()
            if not text or t is getattr(self, "banner", None):
                continue
            items.append((f"text:{text}", f"label '{text}'", t))
        return items

    def layout_conflicts(self, tol: float = 0.002) -> list[tuple[str, str, float, float]]:
        """Overlapping pairs of controls, as ``(a, b, width, height)`` in figure units.

        Measures what is actually *rendered* -- including axis labels and the
        value text a slider draws outside its own axes -- rather than the
        declared rectangles, because the rendered extent is what a user sees
        collide.  An empty list means the window is legible.
        """
        from matplotlib.transforms import Bbox

        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        inv = self.fig.transFigure.inverted()

        boxes = []
        for group, name, art in self._layout_items():
            try:
                bb = (art.get_tightbbox(renderer) if hasattr(art, "get_tightbbox")
                      else art.get_window_extent(renderer))
            except Exception:                      # pragma: no cover - backend quirks
                continue
            if bb is None or bb.width <= 0 or bb.height <= 0:
                continue
            boxes.append((group, name, Bbox(inv.transform(bb))))

        conflicts = []
        for i in range(len(boxes)):
            gi, ni, bi = boxes[i]
            for j in range(i + 1, len(boxes)):
                gj, nj, bj = boxes[j]
                if gi == gj:
                    continue
                inter = Bbox.intersection(bi, bj)
                if inter is None:
                    continue
                if inter.width > tol and inter.height > tol:
                    conflicts.append((ni, nj, inter.width, inter.height))
        return conflicts

    def layout_report(self) -> str:                # pragma: no cover - diagnostic
        c = self.layout_conflicts()
        if not c:
            return "layout OK: no overlapping controls"
        lines = [f"{len(c)} overlapping pair(s):"]
        lines += [f"  {a}  <->  {b}   ({w:.3f} x {h:.3f})" for a, b, w, h in c]
        return "\n".join(lines)

    # ================================================================== #
    def show(self) -> None:
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation

        self.s_rate.valtext.set_text(self._rate_label())
        self.buttons["run"].label.set_text("pause")
        self.anim = FuncAnimation(self.fig, self.update, interval=self.cfg.interval_ms,
                                  blit=False, cache_frame_data=False)
        plt.show()

    def save(self, path: str, frames: int = 200, fps: int = 25, dpi: int = 100,
             script=None) -> None:
        """Record without a window.  ``script(dashboard, frame)`` can drive the controls."""
        from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter

        self.s_rate.valtext.set_text(self._rate_label())

        def step(k):
            if script is not None:
                script(self, k)
            return self.update(k)

        anim = FuncAnimation(self.fig, step, frames=frames, blit=False,
                             cache_frame_data=False)
        writer = (PillowWriter(fps=fps) if str(path).lower().endswith(".gif")
                  else FFMpegWriter(fps=fps, bitrate=2800))
        anim.save(path, writer=writer, dpi=dpi)


__all__ = ["Dashboard", "DashboardConfig", "Frame", "COLOR_MODES", "PANELS", "RATES"]
