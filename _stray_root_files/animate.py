"""Live viewer: run the simulation and watch it, with controls.

This restores -- and extends -- the one thing the original script did that this
package had dropped, its ``plt.ion()`` scatter loop.  The differences that
matter:

* the physics runs at the correct timestep, so what is on screen is a real
  trajectory rather than an exploding one;
* the configuration panel sits beside live traces of temperature, energy and
  hexatic order, so the thing being watched can be *read*;
* atoms are coloured by a physical quantity (local hexatic order, Voronoi
  coordination, speed, per-atom potential energy) instead of at random;
* the target temperature can be driven from the keyboard while it runs, which
  makes melting and recrystallisation something you can do by hand.

Headless use is supported: ``--save out.gif`` records a fixed number of frames
through Pillow or ffmpeg and never opens a window.

Keys
----
``up`` / ``down``   target temperature +/- 50 K   (``shift`` for 250 K)
``c``               cycle the colouring
``space``           pause / resume
``t``               re-draw Maxwell-Boltzmann velocities at the target T
``m``               minimise (quench to the nearest local minimum)
``r``               reset to the starting configuration
``s``               save a PNG snapshot
``q``               quit
"""

from __future__ import annotations

import time
from collections import deque
from dataclasses import dataclass

import matplotlib as mpl
import numpy as np

from . import observables as obs
from .plotting import COLORS, use_style
from .simulation import Simulation

COLOR_MODES = ("psi6", "coord", "speed", "energy")

_CMAP = {
    "psi6": ("viridis", r"$|\psi_6|$", 0.0, 1.0),
    "coord": ("coolwarm", "Voronoi coordination", 4.0, 8.0),
    "speed": ("inferno", "speed (A/ps)", None, None),
    "energy": ("magma", "potential energy (eV/atom)", None, None),
}


@dataclass
class ViewerConfig:
    steps_per_frame: int = 20
    history: int = 400
    color_mode: str = "psi6"
    psi6_cutoff_factor: float = 1.35
    point_size: float = 26.0
    interval_ms: int = 20


class LiveViewer:
    """Drives a :class:`~trilattice.simulation.Simulation` and renders it."""

    def __init__(self, sim: Simulation, a: float, config: ViewerConfig | None = None):
        import matplotlib.pyplot as plt

        self.sim = sim
        self.a = a
        self.cfg = config or ViewerConfig()
        if self.cfg.color_mode not in COLOR_MODES:
            raise ValueError(f"color_mode must be one of {COLOR_MODES}")
        self.paused = False
        self._quit = False
        self._initial = (sim.positions.copy(), sim.velocities.copy())
        self._last_wall = time.perf_counter()
        self._sps = 0.0

        n = self.cfg.history
        self.t_hist: deque = deque(maxlen=n)
        self.temp_hist: deque = deque(maxlen=n)
        self.epot_hist: deque = deque(maxlen=n)
        self.etot_hist: deque = deque(maxlen=n)
        self.psi6_hist: deque = deque(maxlen=n)

        use_style()
        self.fig = plt.figure(figsize=(12.0, 6.2))
        gs = self.fig.add_gridspec(3, 2, width_ratios=[1.35, 1.0], hspace=0.42, wspace=0.28)
        self.ax_cfg = self.fig.add_subplot(gs[:, 0])
        self.ax_t = self.fig.add_subplot(gs[0, 1])
        self.ax_e = self.fig.add_subplot(gs[1, 1])
        self.ax_p = self.fig.add_subplot(gs[2, 1])

        self._setup_axes()
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

    # ------------------------------------------------------------------ #
    @property
    def target_temperature(self) -> float:
        return float(getattr(self.sim.thermostat, "temperature", self.sim.temperature))

    @target_temperature.setter
    def target_temperature(self, value: float) -> None:
        if hasattr(self.sim.thermostat, "temperature"):
            self.sim.thermostat.temperature = max(1.0, float(value))

    # ------------------------------------------------------------------ #
    def _setup_axes(self) -> None:
        import matplotlib.pyplot as plt

        box = self.sim.box
        cmap, label, lo, hi = _CMAP[self.cfg.color_mode]
        self.scatter = self.ax_cfg.scatter(
            self.sim.positions[:, 0],
            self.sim.positions[:, 1],
            c=np.zeros(self.sim.n_particles),
            cmap=cmap,
            vmin=lo,
            vmax=hi,
            s=self.cfg.point_size,
            linewidths=0.3,
            edgecolors="0.25",
        )
        self.ax_cfg.add_patch(
            plt.Rectangle((0, 0), box.lx, box.ly, fill=False, ec="0.45", lw=0.8, ls="--")
        )
        self.ax_cfg.set_xlim(-0.04 * box.lx, 1.04 * box.lx)
        self.ax_cfg.set_ylim(-0.04 * box.ly, 1.04 * box.ly)
        self.ax_cfg.set_aspect("equal")
        self.ax_cfg.grid(False)
        self.ax_cfg.set_xlabel("x (A)")
        self.ax_cfg.set_ylabel("y (A)")
        self.cbar = self.fig.colorbar(self.scatter, ax=self.ax_cfg, fraction=0.040, pad=0.015)
        # the label goes *above* the bar: as a side label it lands on the middle
        # panel's y-axis in this layout
        self.cbar.ax.set_title(label, fontsize=8, pad=7)

        self.readout = self.ax_cfg.text(
            0.012, 0.985, "", transform=self.ax_cfg.transAxes, va="top", ha="left",
            fontsize=8, family="monospace",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.75", alpha=0.88),
        )

        (self.line_t,) = self.ax_t.plot([], [], color=COLORS[1], lw=1.2)
        self.line_target = self.ax_t.axhline(self.target_temperature, color="0.5", ls=":", lw=1.0)
        self.ax_t.set_ylabel("$T$ (K)")

        (self.line_epot,) = self.ax_e.plot([], [], color=COLORS[0], lw=1.2, label=r"$E_\mathrm{pot}$")
        (self.line_etot,) = self.ax_e.plot([], [], color=COLORS[2], lw=1.0, label=r"$E_\mathrm{tot}$")
        self.ax_e.set_ylabel("eV/atom")
        self.ax_e.legend(loc="lower right", ncol=2, fontsize=7, framealpha=0.85)

        (self.line_p,) = self.ax_p.plot([], [], color=COLORS[3], lw=1.2)
        self.ax_p.set_ylabel(r"$|\langle\psi_6\rangle|$")
        self.ax_p.set_ylim(0, 1.02)
        self.ax_p.set_xlabel("time (ps)")

        self.fig.suptitle(
            f"N = {self.sim.n_particles},  "
            f"{self.sim.box.lx:.1f} x {self.sim.box.ly:.1f} A,  "
            f"dt = {self.sim.dt * 1000:.1f} fs,  thermostat: {self.sim.thermostat.name}"
            "     [up/down] T   [c] colour   [space] pause   [m] quench   [q] quit",
            fontsize=9,
        )

    # ------------------------------------------------------------------ #
    def _colors(self) -> tuple[np.ndarray, float | None, float | None]:
        mode = self.cfg.color_mode
        if mode == "psi6":
            v = np.abs(obs.psi6(self.sim.positions, self.sim.box,
                                self.cfg.psi6_cutoff_factor * self.a))
            return v, 0.0, 1.0
        if mode == "coord":
            v = obs.coordination_by_delaunay(self.sim.positions, self.sim.box).astype(float)
            return v, 4.0, 8.0
        if mode == "speed":
            v = np.linalg.norm(self.sim.velocities, axis=1)
            return v, 0.0, float(np.percentile(v, 99)) or 1.0
        v = self.sim.per_atom_energy
        return v, float(np.percentile(v, 1)), float(np.percentile(v, 99))

    def _rescale_colorbar(self, lo, hi) -> None:
        if lo is not None and hi is not None and hi > lo:
            self.scatter.set_clim(lo, hi)

    # ------------------------------------------------------------------ #
    def step(self) -> None:
        """Advance the simulation by one frame's worth of steps."""
        t0 = time.perf_counter()
        for _ in range(self.cfg.steps_per_frame):
            self.sim.step()
        dt_wall = time.perf_counter() - t0
        if dt_wall > 0:
            self._sps = 0.8 * self._sps + 0.2 * self.cfg.steps_per_frame / dt_wall

    def update(self, _frame=None):
        if not self.paused and not self._quit:
            self.step()

        sim = self.sim
        n = sim.n_particles
        p6 = obs.global_psi6(sim.positions, sim.box, self.cfg.psi6_cutoff_factor * self.a)

        self.t_hist.append(sim.elapsed)
        self.temp_hist.append(sim.temperature)
        self.epot_hist.append(sim.e_pot / n)
        self.etot_hist.append(sim.total_energy / n)
        self.psi6_hist.append(p6)

        colors, lo, hi = self._colors()
        self.scatter.set_offsets(sim.positions)
        self.scatter.set_array(colors)
        self._rescale_colorbar(lo, hi)

        t = np.asarray(self.t_hist)
        self.line_t.set_data(t, self.temp_hist)
        self.line_target.set_ydata([self.target_temperature] * 2)
        self.line_epot.set_data(t, self.epot_hist)
        self.line_etot.set_data(t, self.etot_hist)
        self.line_p.set_data(t, self.psi6_hist)
        for ax in (self.ax_t, self.ax_e, self.ax_p):
            ax.relim()
            ax.autoscale_view(scalex=True, scaley=(ax is not self.ax_p))
            if len(t) > 1 and t[-1] > t[0]:
                ax.set_xlim(t[0], t[-1])

        self.readout.set_text(
            f"t      {sim.elapsed:8.2f} ps\n"
            f"T      {sim.temperature:8.1f} K  -> {self.target_temperature:.0f}\n"
            f"E_pot  {sim.e_pot / n:8.4f} eV/at\n"
            f"E_tot  {sim.total_energy / n:8.4f} eV/at\n"
            f"|psi6| {p6:8.3f}\n"
            f"P      {sim.pressure:8.4f} eV/A^2\n"
            f"speed  {self._sps:8.0f} steps/s"
            + ("\nPAUSED" if self.paused else "")
        )
        return (self.scatter, self.line_t, self.line_epot, self.line_etot,
                self.line_p, self.readout)

    # ------------------------------------------------------------------ #
    def _on_key(self, event) -> None:
        import matplotlib.pyplot as plt

        key = event.key or ""
        if key in ("up", "down", "shift+up", "shift+down"):
            delta = 250.0 if "shift" in key else 50.0
            self.target_temperature += delta if "up" in key else -delta
        elif key == "c":
            i = COLOR_MODES.index(self.cfg.color_mode)
            self.cfg.color_mode = COLOR_MODES[(i + 1) % len(COLOR_MODES)]
            cmap, label, lo, hi = _CMAP[self.cfg.color_mode]
            self.scatter.set_cmap(cmap)
            self.cbar.ax.set_title(label, fontsize=8, pad=7)
            if lo is not None:
                self.scatter.set_clim(lo, hi)
        elif key == " ":
            self.paused = not self.paused
        elif key == "t":
            self.sim.set_temperature(self.target_temperature)
        elif key == "m":
            self.sim.minimize(max_steps=400, f_tol=1e-6)
        elif key == "r":
            self.sim.positions[:] = self._initial[0]
            self.sim.velocities[:] = self._initial[1]
            self.sim.unwrapped[:] = self._initial[0]
            self.sim.neighbors.update(self.sim.positions, force=True)
            self.sim._evaluate()
            for d in (self.t_hist, self.temp_hist, self.epot_hist,
                      self.etot_hist, self.psi6_hist):
                d.clear()
            self.sim.elapsed = 0.0
        elif key == "s":
            name = f"snapshot_{self.sim.step_count:08d}.png"
            self.fig.savefig(name, dpi=160)
            print(f"wrote {name}")
        elif key == "q":
            self._quit = True
            plt.close(self.fig)

    # ------------------------------------------------------------------ #
    def show(self) -> None:
        """Open the interactive window (blocks until closed)."""
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation

        self.anim = FuncAnimation(
            self.fig, self.update, interval=self.cfg.interval_ms,
            blit=False, cache_frame_data=False,
        )
        plt.show()

    def save(self, path: str, frames: int = 200, fps: int = 25, dpi: int = 110) -> None:
        """Record ``frames`` frames to a GIF or MP4 without opening a window."""
        from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter

        anim = FuncAnimation(self.fig, self.update, frames=frames,
                             blit=False, cache_frame_data=False)
        writer = PillowWriter(fps=fps) if str(path).lower().endswith(".gif") \
            else FFMpegWriter(fps=fps, bitrate=2400)
        anim.save(path, writer=writer, dpi=dpi)


def is_headless() -> bool:
    """True when matplotlib has no interactive backend available."""
    return mpl.get_backend().lower() in ("agg", "pdf", "ps", "svg", "template")


__all__ = ["LiveViewer", "ViewerConfig", "COLOR_MODES", "is_headless"]
