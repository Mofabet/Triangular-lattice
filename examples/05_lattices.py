"""Figure 7 -- the four lattice types, with and without directional bonding.

Row 1: the ideal structures, bonds drawn, atoms coloured by coordination.
Row 2: after 5 ps at 300 K with the pair potential alone.
Row 3: after 5 ps at 300 K with the three-body angular term switched on.

Only the triangular lattice is a mechanically stable ground state of an isotropic
pair potential, so in row 2 the other three break into close-packed islands --
and the energy *falls* every time, because close packing simply wins when there
is nothing to pay for bond angles.  Row 3 adds that cost: a Stillinger-Weber
style penalty on bond angles away from the ones each structure wants.  It is
exactly zero for the ideal lattice, so it changes no equilibrium property; it
only makes departures expensive.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import trilattice as tl
from trilattice import observables as obs
from trilattice.lattices import LATTICES, build_lattice, cells_for
from trilattice.plotting import COLORS, use_style
from trilattice.threebody import ThreeBodyAngular

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
use_style()

A = 3.3567          # rho* = 0.85 for the triangular case
T_TEST = 300.0
LAMBDA = 0.5        # eV, the angular term strength
PS = 5.0
pot = tl.LennardJones(0.27, 2.88, 8.0, mode="shifted-force")

fig, axes = plt.subplots(3, 4, figsize=(13.5, 10.4))


def draw(ax, positions, box, coord, z0, title, cut):
    i, j, d = obs.neighbour_pairs_within(box, positions, cut)
    raw = positions[i] - positions[j]
    keep = np.all(np.abs(raw - d) < 1e-9, axis=1)
    ax.add_collection(LineCollection(
        np.stack([positions[i[keep]], positions[j[keep]]], axis=1),
        colors="0.75", linewidths=0.7, zorder=1))
    col = ["0.72" if c == z0 else COLORS[1] if c < z0 else COLORS[0] for c in coord]
    ax.scatter(positions[:, 0], positions[:, 1], c=col, s=16,
               linewidths=0.3, edgecolors="0.3", zorder=2)
    ax.add_patch(plt.Rectangle((0, 0), box.lx, box.ly, fill=False,
                               ec="0.45", lw=0.8, ls="--"))
    ax.set_xlim(-2, box.lx + 2)
    ax.set_ylim(-2, box.ly + 2)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=8.5)


print(f"{PS:.0f} ps at {T_TEST:.0f} K, lambda = {LAMBDA} eV\n")
for k, (name, spec) in enumerate(LATTICES.items()):
    nx, ny = cells_for(name, 300)
    cfg = build_lattice(name, nx, ny, A)
    cut = spec.cutoff(A)

    z_ideal = obs.coordination_by_cutoff(cfg.positions, cfg.box, cut)
    psi0 = obs.global_psi_n(cfg.positions, cfg.box, cut, spec.psi_order)
    sim = tl.Simulation(cfg, pot, timestep=0.002,
                        thermostat=tl.Bussi(temperature=T_TEST, tau=0.2), seed=3)
    e0 = sim.e_pot / sim.n_particles
    draw(axes[0, k], cfg.positions, cfg.box, z_ideal, spec.coordination,
         f"{name}   N = {cfg.n_particles}\n"
         rf"$z$ = {spec.coordination},  $|\psi_{{{spec.psi_order}}}|$ = {psi0:.3f},  "
         f"$E$ = {e0:.3f} eV/at", cut)

    for row, term in ((1, None), (2, ThreeBodyAngular.for_lattice(spec, A, LAMBDA))):
        run = tl.Simulation(build_lattice(name, nx, ny, A), pot, timestep=0.002,
                            three_body=term,
                            thermostat=tl.Bussi(temperature=T_TEST, tau=0.2), seed=3)
        run.set_temperature(T_TEST)
        run.run(int(PS / 0.002), log_every=0)
        z_end = obs.coordination_by_cutoff(run.positions, run.box, cut)
        psi1 = obs.global_psi_n(run.positions, run.box, cut, spec.psi_order)
        e1 = run.e_pot / run.n_particles
        verdict = ("intact" if abs(z_end.mean() - spec.coordination) < 0.3 and psi1 > 0.8
                   else "collapsed")
        tag = "pair potential only" if term is None else rf"$\lambda$ = {LAMBDA} eV"
        draw(axes[row, k], run.positions, run.box, z_end, spec.coordination,
             rf"{PS:.0f} ps at {T_TEST:.0f} K, {tag} -- {verdict}" "\n"
             rf"$\bar z$ = {z_end.mean():.2f},  $|\psi_{{{spec.psi_order}}}|$ = {psi1:.3f},  "
             f"$E$ = {e1:.3f} eV/at", cut)
        print(f"{name:11s} {'3-body' if term else 'pair  '} "
              f"z={z_end.mean():5.2f} psi={psi1:6.3f} E={e1:8.4f}  {verdict}")

fig.suptitle(
    "2-D lattices: a pair potential keeps only close packing (row 2); "
    "an angular term keeps them all (row 3)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT / "fig7_lattices.png")
print("\nwrote", OUT / "fig7_lattices.png")
