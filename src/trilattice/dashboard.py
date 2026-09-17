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

    def __init__(self, sim: Simulation, a: float, config: DashboardConfig | None = None):
        import matplotlib.pyplot as plt

        self.sim = sim
        self.a = a
        self.cfg = config or DashboardConfig()
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
        self.ax_cfg = f.add_axes([0.035, 0.30, 0.40, 0.62])
        self.ax_panel = f.add_axes([0.53, 0.665, 0.29, 0.255])
        self.ax_t = f.add_axes([0.53, 0.435, 0.29, 0.175])
        self.ax_e = f.add_axes([0.53, 0.205, 0.29, 0.175])

    def _build_widgets(self) -> None:
        from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider

        f = self.fig
        self.s_temp = Slider(
            f.add_axes([0.075, 0.175, 0.30, 0.022]), "target T (K)",
            self.cfg.t_min, self.cfg.t_max, valinit=self.target_temperature, valfmt="%.0f",
            color=COLORS[1],
        )
        self.s_rate = Slider(
            f.add_axes([0.075, 0.130, 0.30, 0.022]), "rate", 0, len(RATES) - 1,
            valinit=self.rate_index, valstep=1, color=COLORS[0],
        )
        self.s_scrub = Slider(
            f.add_axes([0.075, 0.085, 0.30, 0.022]), "time", 0.0, 1.0,
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
            ax = f.add_axes([0.075 + col * 0.063, 0.038 - row * 0.030, 0.058, 0.024])
            b = Button(ax, label, color="0.94", hovercolor="0.86")
            b.label.set_fontsize(7)
            b.on_clicked(cb)
            self.buttons[label] = b

        self.r_color = RadioButtons(
            f.add_axes([0.850, 0.685, 0.140, 0.235]),
            COLOR_MODES, active=COLOR_MODES.index(self.cfg.color_mode),
            label_props={"fontsize": [9.5]},
        )
        self.r_panel = RadioButtons(
            f.add_axes([0.850, 0.395, 0.140, 0.255]),
            PANELS, active=PANELS.index(self.cfg.panel),
            label_props={"fontsize": [9.5]},
        )
        for r in (self.r_color, self.r_panel):
            r.ax.set_facecolor("none")
            r.ax.set_frame_on(False)
        self.r_color.on_clicked(self._on_color)
        self.r_panel.on_clicked(self._on_panel)

        self.checks = CheckButtons(
            f.add_axes([0.850, 0.215, 0.140, 0.145]),
            ["bonds", "trails", "box", "thermostat"], [False, False, True, True],
            label_props={"fontsize": [9.5]},
        )
        self.checks.ax.set_facecolor("none")
        self.checks.ax.set_frame_on(False)
        self.checks.on_clicked(self._on_check)

        self.fig.text(0.850, 0.925, "colour by", fontsize=8, weight="bold")
        self.fig.text(0.850, 0.655, "analysis panel", fontsize=8, weight="bold")
        self.fig.text(0.850, 0.365, "overlays", fontsize=8, weight="bold")

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
        self.ax_e.set_ylabel(r"$|\langle\psi_6\rangle|$")
        self.ax_e.set_ylim(0, 1.02)
        self.ax_e.set_xlabel("time (ps)")
        self.ax_def = self.ax_e.twinx()
        (self.line_def,) = self.ax_def.plot([], [], color=COLORS[0], lw=1.0,
                                            label="defects")
        self.ax_def.set_ylabel("defects (%)", color=COLORS[0])
        self.ax_def.tick_params(axis="y", labelcolor=COLORS[0], labelsize=7)
        self.ax_def.grid(False)

        self.fig.suptitle(
            f"trilattice   N = {self.sim.n_particles}   "
            f"{self.sim.box.lx:.1f} x {self.sim.box.ly:.1f} A   "
            f"dt = {self.sim.dt * 1000:.1f} fs",
            fontsize=10, x=0.235,
        )

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
            return f"1 step / {int(round(1 / r))} frames"
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
        """Insert at the emptiest hollow, then relax with a capped FIRE step.

        Insertion always overlaps -- the deepest hole of a triangular lattice is
        a/sqrt(3) from its neighbours, inside sigma -- so the relaxation is not
        optional.  It zeroes the velocities, hence the re-thermalisation.
        """
        t_before = self.sim.temperature
        sites = hollow_sites(self.sim.positions, self.sim.box, self.a)
        d = self.sim.box.minimum_image(sites[:, None, :] - self.sim.positions[None, :, :])
        best = sites[np.argmax(np.min(np.hypot(d[:, :, 0], d[:, :, 1]), axis=1))]
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
        elif label == "thermostat":
            if self.thermostat_on:
                self._live_thermostat = self.sim.thermostat
                self.sim.thermostat = NoThermostat()
                self._notify("thermostat off -- NVE")
            else:
                self.sim.thermostat = self._live_thermostat
                self._notify(f"thermostat on -- {self.sim.thermostat.name}")

    def _checked(self, name: str) -> bool:
        return bool(self.checks.get_status()[["bonds", "trails", "box", "thermostat"].index(name)])

    def _resize_visuals(self) -> None:
        """Re-point everything that is sized by N after an insertion/removal."""
        n = self.sim.n_particles
        self.msd_origin = self.sim.unwrapped.copy()
        self.msd_origin_time = self.sim.elapsed
        self.trail_ids = np.linspace(0, n - 1, min(self.cfg.n_trails, n)).astype(int)
        self.trails.clear()
        self.frames.clear()
        self._heavy_coord = np.full(n, 6.0)
        self.fig.suptitle(
            f"trilattice   N = {n}   {self.sim.box.lx:.1f} x {self.sim.box.ly:.1f} A   "
            f"dt = {self.sim.dt * 1000:.1f} fs", fontsize=10, x=0.235,
        )

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
            return np.abs(obs.psi6(pos, box, self.cfg.psi6_cutoff_factor * self.a)), 0.0, 1.0
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
        i, j, d = obs.neighbour_pairs_within(self.sim.box, pos, 1.35 * self.a)
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
        ax.grid(alpha=0.25, lw=0.5)
        sim = self.sim

        if name == "g(r)":
            r, g = self.rdf.result()
            ax.plot(r, g, color=COLORS[0], lw=1.1)
            ax.axhline(1.0, color="0.6", lw=0.7, ls=":")
            for shell in (1.0, np.sqrt(3), 2.0):
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
                ax.set_aspect("equal")
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
            cols = [COLORS[1] if v == 5 else COLORS[0] if v == 7
                    else "0.7" if v == 6 else COLORS[3] for v in values]
            ax.bar(values, 100 * counts / max(1, c.size), color=cols, width=0.7)
            ax.set_xlabel("Voronoi coordination")
            ax.set_ylabel("% of atoms")
            ax.set_xlim(2.5, 9.5)
            ax.set_title(f"non-six-coordinated: {100 * self._heavy_defects:.1f} %", fontsize=8.5)

        else:  # T-psi6
            if self.frames:
                t = np.array([f.temperature for f in self.frames])
                p = np.array([f.psi6 for f in self.frames])
                ax.scatter(t, p, c=np.arange(len(t)), cmap="viridis", s=6, linewidths=0)
                ax.plot(t[-1], p[-1], "o", color=COLORS[1], ms=6)
            ax.set_xlabel("$T$ (K)")
            ax.set_ylabel(r"$|\langle\psi_6\rangle|$")
            ax.set_ylim(0, 1.02)
            ax.set_title("melting curve traced live", fontsize=8.5)

    # ------------------------------------------------------------------ #
    def update(self, _frame=None):
        run = not (self.paused or self.scrubbing)
        if run:
            self._advance()
        self._frame += 1

        sim = self.sim
        heavy = self._frame % self.cfg.heavy_every == 0
        if heavy and run:
            self._heavy_coord = obs.coordination_by_delaunay(
                sim.positions, sim.box).astype(float)
            self._heavy_defects = float(np.mean(self._heavy_coord != 6))
            if self.cfg.panel == "S(k)":
                self._heavy_sk = obs.structure_factor(sim.positions, sim.box, n_max=16)
            self.rdf.accumulate(sim.positions)

        p6 = obs.global_psi6(sim.positions, sim.box, self.cfg.psi6_cutoff_factor * self.a)
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
                    f"|psi6| {p6:8.3f}\n"
                    f"defects{100 * self._heavy_defects:7.1f} %\n"
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
