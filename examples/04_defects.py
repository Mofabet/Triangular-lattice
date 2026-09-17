"""Figure 6 -- point defects, and the throughput of the neighbour machinery.

(a) A relaxed vacancy: the six neighbours of an empty site relax inward, and the
    Voronoi construction reports the ring of 5-coordinated atoms that any
    dislocation analysis keys on.
(b) A relaxed interstitial placed on a three-fold hollow.  This lattice has no
    empty space at all: the deepest hole is only a/sqrt(3) = 1.85 A from its
    three neighbours, well inside sigma = 2.88 A.  So an inserted atom always
    starts on the steep repulsive wall, the hollow site merely being the least
    bad choice -- it is the point that maximises the minimum distance.  Dropping
    it at a uniformly random position, as the original did, is typically an order
    of magnitude worse still and hands the integrator a force it cannot survive.
(c) Formation energies of both defects from a FIRE relaxation.
(d) Cost per step against system size: the linked-cell + Verlet scheme is O(N),
    while the original's nine-replica full matrix is O(81 N^2).
"""
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt

import trilattice as tl
from trilattice import observables as obs
from trilattice.plotting import use_style, COLORS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
use_style()

EPS, SIG, RC, A0 = 0.27, 2.88, 8.0, 3.2010
pot = tl.LennardJones(EPS, SIG, RC, mode="shifted-force")
rng = np.random.default_rng(3)

fig, ax = plt.subplots(1, 4, figsize=(13.5, 3.6))

# ---- reference perfect crystal ------------------------------------------- #
perfect = tl.triangular_lattice(12, 7, A0)
sim0 = tl.Simulation(perfect, pot, timestep=0.002, seed=1)
sim0.minimize(max_steps=4000, f_tol=1e-9)
e_perfect = sim0.e_pot / perfect.n_particles

results = {}
for k, (name, kw) in enumerate([("vacancy", dict(n_vacancies=1)),
                                ("interstitial", dict(n_interstitials=1))]):
    cfg = tl.add_defects(perfect, a=A0, rng=np.random.default_rng(7),
                         min_separation=0.5 * A0, **kw)
    sim = tl.Simulation(cfg, pot, timestep=0.002, seed=2)
    e_before = sim.e_pot
    info = sim.minimize(max_steps=6000, f_tol=1e-9)
    n = sim.n_particles
    # E_f = E(defect) - N_defect/N_perfect * E(perfect)
    e_form = sim.e_pot - n * e_perfect
    results[name] = dict(e_form=e_form, n=n, relax=e_before - sim.e_pot,
                         iters=info["iterations"])

    coord = obs.coordination_by_delaunay(sim.positions, sim.box)
    disp = np.linalg.norm(sim.box.minimum_image(sim.positions - cfg.positions), axis=1)
    a = ax[k]
    col = np.full(n, "0.8", dtype=object)
    col[coord == 5] = COLORS[1]
    col[coord == 7] = COLORS[0]
    col[(coord != 5) & (coord != 6) & (coord != 7)] = COLORS[3]
    a.scatter(sim.positions[:, 0], sim.positions[:, 1], c=list(col),
              s=40 + 900 * disp, linewidths=0.4, edgecolors="0.35")
    a.add_patch(plt.Rectangle((0, 0), sim.box.lx, sim.box.ly, fill=False,
                              ec="0.4", lw=0.8, ls="--"))
    a.set_aspect("equal"); a.grid(False)
    a.set_xlabel("x (A)"); a.set_ylabel("y (A)" if k == 0 else "")
    a.set_title(f"({'ab'[k]}) relaxed {name}\n"
                f"$E_f$ = {e_form:.3f} eV, max relax = {disp.max():.3f} A", fontsize=9)
    print(f"{name:13s} N={n:4d}  E_f = {e_form:+.4f} eV  relaxation = "
          f"{results[name]['relax']:.4f} eV in {info['iterations']} FIRE steps  "
          f"5/7 defects = {np.sum(coord != 6)}")

# ---- (c) what a random insertion costs ------------------------------------ #
a = ax[2]
trials = rng.uniform(0, 1, size=(20000, 2)) * np.array([perfect.box.lx, perfect.box.ly])
d = perfect.box.minimum_image(trials[:, None, :] - perfect.positions[None, :, :])
dmin = np.min(np.hypot(d[:, :, 0], d[:, :, 1]), axis=1)
a.hist(dmin, bins=60, density=True, color=COLORS[0], alpha=0.75,
       label="uniformly random insertion")
hollow = tl.interstitial_sites(perfect, A0)
dh = perfect.box.minimum_image(hollow[:, None, :] - perfect.positions[None, :, :])
dhmin = np.min(np.hypot(dh[:, :, 0], dh[:, :, 1]), axis=1)
a.axvline(float(np.mean(dhmin)), color=COLORS[2], lw=1.4,
          label=f"hollow site ({np.mean(dhmin):.2f} A)")
a.axvline(SIG, color=COLORS[1], ls="--", lw=1.0, label=r"$\sigma$")
u_random = pot.bare_energy(dmin)
u_hollow = float(pot.bare_energy(np.mean(dhmin)))
a.set_xlabel("distance to the nearest host atom (A)")
a.set_ylabel("probability density")
a.set_title("(c) there is no empty space\n"
            f"random median $u$ = {np.median(u_random):.0f} eV vs "
            f"{u_hollow:.0f} eV on the hollow", fontsize=9)
a.legend(fontsize=7)
print(f"\nevery point of the cell is within {dmin.max():.2f} A of an atom "
      f"(sigma = {SIG} A), so 100% of insertions start repulsive")
print(f"bare pair energy at insertion: hollow site {u_hollow:.0f} eV, "
      f"random median {np.median(u_random):.0f} eV, random 90th pct "
      f"{np.percentile(u_random, 90):.3g} eV")

# ---- (d) scaling ---------------------------------------------------------- #
a = ax[3]
sizes, us_step, pairs = [], [], []
for nx, ny in [(6, 4), (10, 6), (16, 9), (22, 13), (30, 18), (40, 24), (52, 30)]:
    cfg = tl.triangular_lattice(nx, ny, A0)
    sim = tl.Simulation(cfg, pot, timestep=0.002, seed=1)
    sim.set_temperature(300.0)
    sim.run(20, log_every=0)
    t0 = time.perf_counter()
    sim.run(120, log_every=0)
    dt = time.perf_counter() - t0
    sizes.append(cfg.n_particles)
    us_step.append(1e6 * dt / 120)
    pairs.append(sim.neighbors.n_pairs)
    print(f"N = {cfg.n_particles:5d}  {us_step[-1]:9.1f} us/step  "
          f"{pairs[-1]:7d} pairs  vs 81N^2 = {81 * cfg.n_particles**2:10d}")
sizes = np.array(sizes, dtype=float)
a.loglog(sizes, us_step, "o-", color=COLORS[0], ms=4, label="this code (cell + Verlet)")
a.loglog(sizes, us_step[0] * (sizes / sizes[0]), "k:", lw=0.9, label=r"$\propto N$")
a.loglog(sizes, us_step[0] * (sizes / sizes[0]) ** 2, "k--", lw=0.9, label=r"$\propto N^2$")
a.set_xlabel("$N$"); a.set_ylabel(r"$\mu$s per step")
a.set_title("(d) cost per step")
a.legend(loc="upper left")

fig.tight_layout()
fig.savefig(OUT / "fig6_defects_and_scaling.png")
slope = np.polyfit(np.log(sizes), np.log(us_step), 1)[0]
print(f"\nmeasured scaling exponent: N^{slope:.2f}")
print("wrote", OUT / "fig6_defects_and_scaling.png")
