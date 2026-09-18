"""Record the README tour: every lattice, every panel, every colour mode.

Six scenes.  Each one switches lattice, analysis panel, colouring and overlays
together, holds the fresh structure for a moment, then raises the temperature to
show what happens to it.  Between them the five colour modes and all six
analysis panels are covered once each, and the directional-bonding switch is
shown holding the open lattices together and then let go, so that honeycomb
collapses on camera.

    python examples/record_tour.py [output.mp4]
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import trilattice as tl
from trilattice.config import Settings
from trilattice.dashboard import COLOR_MODES, PANELS, Dashboard, DashboardConfig
from trilattice.lattices import LATTICE_NAMES

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path(__file__).resolve().parents[1] / "figures" / "lattices.mp4"

#  lattice,     colour,          panel,          T_hot,  overlays
SCENES = [
    ("triangular", "psi6",         "T-psi6",       2400.0, ("bonds",)),
    ("square",     "coord",        "coordination",  900.0, ("bonds", "directional")),
    ("honeycomb",  "energy",       "g(r)",          900.0, ("bonds", "directional")),
    ("kagome",     "displacement", "S(k)",          900.0, ("directional",)),
    ("honeycomb",  "coord",        "coordination",  900.0, ("bonds",)),
    ("triangular", "speed",        "speeds",       2600.0, ("trails",)),
]
HOLD, HEAT = 9, 19
BLOCK = HOLD + HEAT


s = Settings()
s.system.a = 3.3567
s.run.tau = 0.15
sim = tl.Simulation(s.build_configuration(), s.build_potential(), timestep=0.002,
                    thermostat=s.build_thermostat(300.0), seed=1)
sim.set_temperature(300.0)
dash = Dashboard(sim, s.system.a,
                 DashboardConfig(rate_index=8, heavy_every=2),
                 lattice="triangular", settings=s)


def set_overlays(d, wanted):
    for name in ("bonds", "trails"):
        if d._checked(name) != (name in wanted):
            d.checks.set_active(d.OVERLAYS.index(name))
    if d._checked("directional") != ("directional" in wanted):
        d.physics.set_active(d.PHYSICS.index("directional"))


def script(d, frame):
    scene, phase = divmod(frame, BLOCK)
    if scene >= len(SCENES):
        return
    lattice, colour, panel, t_hot, overlays = SCENES[scene]
    if phase == 0:
        d.r_lattice.set_active(LATTICE_NAMES.index(lattice))
        d.r_color.set_active(COLOR_MODES.index(colour))
        d.r_panel.set_active(PANELS.index(panel))
        set_overlays(d, overlays)
        tag = "  |  directional ON" if "directional" in overlays else "  |  pair only"
        d._notify(f"{lattice}  |  colour: {colour}  |  panel: {panel}{tag}", seconds=3.5)
    elif phase == HOLD:
        d.s_temp.set_val(t_hot)


frames = len(SCENES) * BLOCK
dash.save(str(OUT), frames=frames, fps=12, dpi=110, script=script)
print(f"{frames} frames -> {OUT}")
