"""Command-line interface.

    trilattice init           write a commented input deck
    trilattice run            run a simulation from a deck (TOML or legacy start.txt)
    trilattice scan           sweep temperature and report the melting diagnostics
    trilattice animate        interactive dashboard (--minimal for the small viewer)
    trilattice bench          measure throughput and the cost of the old O(N^2) scheme
    trilattice info           show derived quantities for a parameter set
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from . import observables as obs
from .analysis import estimate
from .config import DEFAULT_TOML, Settings
from .forces import HAVE_NUMBA
from .lattice import replicate, required_replication, triangular_lattice
from .potentials import LennardJones
from .simulation import Simulation
from .trajectory import Trajectory
from .units import KB, ReducedUnits


def _build(settings: Settings, temperature: float | None = None,
           replicate_to: tuple[int, int] | None = None, auto: bool = False) -> Simulation:
    t = settings.run.temperature if temperature is None else temperature
    cfg = settings.build_configuration()
    pot = settings.build_potential()

    r_list = pot.cutoff + settings.run.skin
    if replicate_to is None and auto:
        replicate_to = required_replication(cfg.box, r_list)
    if replicate_to is not None and replicate_to != (1, 1):
        cfg = replicate(cfg, *replicate_to)
        print(f"replicated {replicate_to[0]}x{replicate_to[1]} -> N = {cfg.n_particles}, "
              f"box {cfg.box.lx:.2f} x {cfg.box.ly:.2f} A")
    if r_list > cfg.box.max_cutoff:
        nx, ny = required_replication(cfg.box, r_list)
        raise SystemExit(
            f"\nThis cell is too small for its own cutoff.\n"
            f"  box        {cfg.box.lx:.2f} x {cfg.box.ly:.2f} A  "
            f"(minimum image allows r < {cfg.box.max_cutoff:.2f} A)\n"
            f"  r_cut+skin {r_list:.2f} A\n"
            f"An atom would interact with a neighbour and with that neighbour's own\n"
            f"periodic image at the same time.  Either lower the cutoff, or enlarge\n"
            f"the cell:\n\n    --replicate {nx} {ny}      (or --auto-replicate)\n"
        )
    sim = Simulation(
        cfg,
        pot,
        timestep=settings.run.timestep,
        thermostat=settings.build_thermostat(t),
        skin=settings.run.skin,
        seed=settings.run.seed,
    )
    sim.set_temperature(t)
    return sim


# --------------------------------------------------------------------------- #
def cmd_init(args) -> int:
    path = Path(args.output)
    if path.exists() and not args.force:
        print(f"{path} exists; pass --force to overwrite", file=sys.stderr)
        return 1
    path.write_text(DEFAULT_TOML, encoding="utf-8")
    print(f"wrote {path}")
    return 0


def cmd_info(args) -> int:
    s = Settings.load(args.config) if args.config else Settings()
    cfg = s.build_configuration()
    pot = s.build_potential()
    u = ReducedUnits(s.potential.epsilon, s.potential.sigma, s.system.mass)
    print(f"N            = {cfg.n_particles}")
    print(f"box          = {cfg.box.lx:.3f} x {cfg.box.ly:.3f} A   (area {cfg.box.area:.1f} A^2)")
    print(f"density      = {cfg.density:.5f} A^-2   ->  rho* = {u.density(cfg.density):.4f}")
    print(f"potential    = {pot.summary()}")
    r_list = pot.cutoff + s.run.skin
    ok = "ok" if r_list <= cfg.box.max_cutoff else "TOO SMALL -- use --auto-replicate"
    print(f"max cutoff   = {cfg.box.max_cutoff:.2f} A (minimum image limit); "
          f"r_cut+skin = {r_list:.2f} A  [{ok}]")
    print(f"LJ time unit = {u.tau * 1000:.2f} fs;  dt = {s.run.timestep * 1000:.2f} fs "
          f"= tau/{u.tau / s.run.timestep:.0f}")
    print(f"T            = {s.run.temperature:.1f} K  ->  T* = {u.temperature(s.run.temperature):.4f}")
    print(f"E_lattice    = {pot.lattice_energy(s.system.a):.6f} eV/atom at a = {s.system.a} A")
    print(f"backend      = {'numba' if HAVE_NUMBA else 'numpy'}")
    return 0


def cmd_run(args) -> int:
    s = Settings.load(args.config) if args.config else Settings()
    if args.temperature is not None:
        s.run.temperature = args.temperature
    if args.steps is not None:
        s.run.steps = args.steps
    rep = tuple(args.replicate) if args.replicate else None
    sim = _build(s, replicate_to=rep, auto=args.auto_replicate)
    print(sim.summary())

    print(f"\nequilibrating {s.run.equilibration} steps ...")
    sim.run(s.run.equilibration, log_every=0, progress=True)
    sim.clear_log()
    sim.reset_origin()

    traj = Trajectory(sim.box)
    print(f"production {s.run.steps} steps ...")
    log = sim.run(
        s.run.steps,
        log_every=s.run.log_every,
        callback=traj.append,
        callback_every=s.run.sample_every,
        progress=True,
    )

    n = sim.n_particles
    a = s.system.a
    e = estimate(np.array(log.e_pot) / n)
    t = estimate(np.array(log.temperature))
    p = estimate(np.array(log.pressure))
    msd = obs.mean_squared_displacement(traj.u)
    d_coef = obs.diffusion_coefficient(traj.t, msd)

    print("\n--- results " + "-" * 50)
    print(f"temperature      {t.mean:10.2f} +/- {t.error:.2f} K   (tau_int = {t.tau_int:.1f})")
    print(f"potential energy {e.mean:10.6f} +/- {e.error:.6f} eV/atom")
    print(f"pressure         {p.mean:10.6f} +/- {p.error:.6f} eV/A^2")
    print(f"|<psi_6>|        {obs.global_psi6(sim.positions, sim.box, 1.35 * a):10.4f}")
    print(f"defect fraction  {obs.defect_fraction(sim.positions, sim.box):10.4f}")
    print(f"Lindemann        {obs.lindemann_2d(traj.u, sim.box, a):10.4f}")
    print(f"diffusion        {d_coef:10.5f} A^2/ps")
    perf = sim.performance()
    print(f"\nbackend {perf['backend']}, {perf['us_per_step']:.0f} us/step, "
          f"{perf['pairs']} pairs, {perf['rebuild_fraction'] * 100:.1f}% of steps rebuilt the list")

    if args.xyz:
        extra = {
            "psi6": [np.abs(obs.psi6(f, sim.box, 1.35 * a)) for f in traj.positions],
            "coord": [obs.coordination_by_delaunay(f, sim.box).astype(float) for f in traj.positions],
        }
        traj.write_extxyz(args.xyz, extra=extra)
        print(f"wrote {args.xyz} ({len(traj)} frames)")
    if args.npz:
        traj.save_npz(args.npz)
        print(f"wrote {args.npz}")
    return 0


def cmd_scan(args) -> int:
    s = Settings.load(args.config) if args.config else Settings()
    temps = np.linspace(args.t_min, args.t_max, args.n_points)
    a = s.system.a
    rep = tuple(args.replicate) if args.replicate else None
    print(f"{'T/K':>8} {'E/atom':>12} {'P':>10} {'|psi6|':>8} {'defects':>9} {'D':>10}")
    sim = None
    for T in temps:
        if sim is None:
            sim = _build(s, float(T), replicate_to=rep, auto=args.auto_replicate)
        else:
            sim.thermostat = s.build_thermostat(float(T))
        sim.run(s.run.equilibration, log_every=0)
        sim.clear_log()
        sim.reset_origin()
        frames: list[np.ndarray] = []
        log = sim.run(
            s.run.steps,
            log_every=s.run.log_every,
            callback=lambda x: frames.append(x.unwrapped.copy()),
            callback_every=s.run.sample_every,
        )
        traj = np.asarray(frames)
        times = np.arange(traj.shape[0]) * s.run.sample_every * s.run.timestep
        d_coef = obs.diffusion_coefficient(times, obs.mean_squared_displacement(traj))
        print(
            f"{np.mean(log.temperature):8.1f} {np.mean(log.e_pot) / sim.n_particles:12.6f} "
            f"{np.mean(log.pressure):10.5f} {obs.global_psi6(sim.positions, sim.box, 1.35 * a):8.4f} "
            f"{obs.defect_fraction(sim.positions, sim.box):9.4f} {d_coef:10.5f}",
            flush=True,
        )
    return 0


def cmd_animate(args) -> int:
    from .animate import LiveViewer, ViewerConfig, is_headless
    from .dashboard import Dashboard, DashboardConfig, RATES

    s = Settings.load(args.config) if args.config else Settings()
    if args.temperature is not None:
        s.run.temperature = args.temperature
    if s.run.thermostat in ("none", "nve"):
        s.run.thermostat = "bussi"   # the temperature keys need something to drive
    rep = tuple(args.replicate) if args.replicate else None
    sim = _build(s, replicate_to=rep, auto=args.auto_replicate)
    print(sim.summary())

    if args.minimal:
        viewer = LiveViewer(
            sim, s.system.a,
            ViewerConfig(steps_per_frame=args.steps_per_frame,
                         color_mode=args.color, history=args.history),
        )
    else:
        rate_index = int(np.argmin(np.abs(np.array(RATES) - args.steps_per_frame)))
        viewer = Dashboard(
            sim, s.system.a,
            DashboardConfig(rate_index=rate_index, color_mode=args.color,
                            panel=args.panel, history=args.history),
        )
    if args.save:
        print(f"recording {args.frames} frames to {args.save} ...")
        viewer.save(args.save, frames=args.frames, fps=args.fps, dpi=args.dpi)
        print(f"wrote {args.save}")
        return 0
    if is_headless():
        raise SystemExit(
            f"\nmatplotlib is using the non-interactive '{__import__('matplotlib').get_backend()}' "
            "backend, so no window can be opened.\n"
            "On a machine with a display, install a GUI backend (for example\n"
            "`pip install PyQt5`).  Otherwise record a file instead:\n\n"
            "    trilattice animate --save melting.gif --frames 240\n"
        )
    if args.minimal:
        print("\nkeys: up/down = target T (+shift for 250 K), c = colouring, "
              "space = pause,\n      t = redraw velocities, m = quench, r = reset, "
              "s = snapshot, q = quit\n")
    else:
        print("\nsliders: target T, rate (down to one step per eight frames), "
              "time scrub\nbuttons: pause, live, resume here, reset, quench, reheat,"
              " +vacancy,\n         +interstitial, new MSD origin, clear g(r)\n"
              "panels:  g(r), S(k), MSD, speeds, coordination, T-psi6\n")
    viewer.show()
    return 0


def cmd_bench(args) -> int:
    print(f"backend: {'numba' if HAVE_NUMBA else 'numpy'}\n")
    print(f"{'N':>7} {'pairs':>9} {'us/step':>10} {'ns/step/atom':>13} {'naive pairs':>13} {'speedup':>9}")
    pot = LennardJones(0.27, 2.88, 8.0)
    for nx, ny in [(6, 4), (10, 6), (16, 9), (22, 13), (30, 18), (40, 24)]:
        cfg = triangular_lattice(nx, ny, 3.2)
        sim = Simulation(cfg, pot, timestep=0.002, seed=1)
        sim.set_temperature(300.0)
        sim.run(20, log_every=0)
        t0 = time.perf_counter()
        sim.run(args.steps, log_every=0)
        dt = time.perf_counter() - t0
        n = cfg.n_particles
        naive = 81 * n * n  # the original built 9 replicas and a full 9N x 9N matrix
        print(
            f"{n:7d} {sim.neighbors.n_pairs:9d} {1e6 * dt / args.steps:10.1f} "
            f"{1e9 * dt / args.steps / n:13.1f} {naive:13d} {naive / max(1, sim.neighbors.n_pairs):9.0f}x"
        )
    return 0


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="trilattice", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    q = sub.add_parser("init", help="write a commented input deck")
    q.add_argument("-o", "--output", default="config.toml")
    q.add_argument("--force", action="store_true")
    q.set_defaults(func=cmd_init)

    q = sub.add_parser("info", help="derived quantities for a parameter set")
    q.add_argument("config", nargs="?")
    q.set_defaults(func=cmd_info)

    q = sub.add_parser("run", help="run one simulation")
    q.add_argument("config", nargs="?")
    q.add_argument("-T", "--temperature", type=float)
    q.add_argument("-n", "--steps", type=int)
    q.add_argument("--xyz", help="write an extended-XYZ trajectory")
    q.add_argument("--npz", help="write a compressed npz trajectory")
    q.add_argument("--replicate", nargs=2, type=int, metavar=("NX", "NY"),
                   help="tile the cell before running")
    q.add_argument("--auto-replicate", action="store_true",
                   help="tile just enough to satisfy the minimum-image condition")
    q.set_defaults(func=cmd_run)

    q = sub.add_parser("scan", help="sweep temperature")
    q.add_argument("config", nargs="?")
    q.add_argument("--t-min", type=float, default=200.0)
    q.add_argument("--t-max", type=float, default=3000.0)
    q.add_argument("--n-points", type=int, default=15)
    q.add_argument("--replicate", nargs=2, type=int, metavar=("NX", "NY"))
    q.add_argument("--auto-replicate", action="store_true")
    q.set_defaults(func=cmd_scan)

    q = sub.add_parser("animate", help="live window with dynamics and traces")
    q.add_argument("config", nargs="?")
    q.add_argument("-T", "--temperature", type=float)
    q.add_argument("--color", default="psi6",
                   choices=("psi6", "coord", "speed", "energy", "displacement"),
                   help="what to colour the atoms by")
    q.add_argument("--panel", default="g(r)",
                   choices=("g(r)", "S(k)", "MSD", "speeds", "coordination", "T-psi6"),
                   help="which live analysis panel to open with")
    q.add_argument("--minimal", action="store_true",
                   help="the small three-panel viewer instead of the full dashboard")
    q.add_argument("--steps-per-frame", type=int, default=20)
    q.add_argument("--history", type=int, default=400, help="points kept in the traces")
    q.add_argument("--replicate", nargs=2, type=int, metavar=("NX", "NY"))
    q.add_argument("--auto-replicate", action="store_true")
    q.add_argument("--save", help="record to a .gif or .mp4 instead of opening a window")
    q.add_argument("--frames", type=int, default=200)
    q.add_argument("--fps", type=int, default=25)
    q.add_argument("--dpi", type=int, default=100)
    q.set_defaults(func=cmd_animate)

    q = sub.add_parser("bench", help="throughput vs system size")
    q.add_argument("--steps", type=int, default=200)
    q.set_defaults(func=cmd_bench)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
