"""Figure 2 -- does the integrator do what it claims?

(a) NVE energy fluctuation vs timestep: velocity Verlet is second order, so the
    RMS energy error must fall as dt^2.  A straight line of slope 2 on log-log
    is the signature; anything else means the force is not the gradient of the
    energy.
(b) Long NVE run: the energy error of a symplectic integrator is *bounded*, it
    does not accumulate.  Total momentum stays at machine zero.
(c) Kinetic-energy distribution for each thermostat against the exact canonical
    result, chi^2 with N_f degrees of freedom.  Berendsen visibly fails: its
    distribution is far too narrow, which is why it gives the wrong heat
    capacity.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

import trilattice as tl
from trilattice.plotting import use_style, COLORS
from trilattice.units import KB

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
use_style()

EPS, SIG, RC = 0.27, 2.88, 8.0
A0 = 3.2010
cfg = tl.triangular_lattice(12, 7, A0)
N = cfg.n_particles
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5))

# ---- (a) order of the integrator ----------------------------------------- #
ax = axes[0]
dts = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
for k, mode in enumerate(("shifted", "shifted-force")):
    pot = tl.LennardJones(EPS, SIG, RC, mode=mode)
    rms = []
    for dtfs in dts:
        sim = tl.Simulation(cfg, pot, timestep=dtfs / 1000, seed=7)
        sim.set_temperature(400.0)
        sim.run(3000, log_every=0)
        sim.clear_log()
        log = sim.run(8000, log_every=10)
        rms.append(np.std(log.e_tot) / N)
    rms = np.array(rms)
    ax.loglog(dts, rms, "o-", color=COLORS[k], ms=4, label=mode)
    p = np.polyfit(np.log(dts), np.log(rms), 1)
    ax.text(0.05, 0.9 - 0.09 * k, f"{mode}: slope = {p[0]:.2f}",
            transform=ax.transAxes, fontsize=8, color=COLORS[k])
ax.loglog(dts, 3e-7 * (dts / 0.25) ** 2, "k:", lw=0.9, label=r"$\propto dt^2$")
ax.set_xlabel("timestep (fs)")
ax.set_ylabel(r"RMS $\delta E_\mathrm{tot}$ (eV/atom)")
ax.set_title(f"(a) NVE, velocity Verlet, N = {N}")
ax.legend(loc="lower right")

# ---- (b) bounded drift ---------------------------------------------------- #
ax = axes[1]
pot = tl.LennardJones(EPS, SIG, RC, mode="shifted-force")
sim = tl.Simulation(cfg, pot, timestep=0.002, seed=11)
sim.set_temperature(400.0)
sim.run(5000, log_every=0)
sim.clear_log()
log = sim.run(150000, log_every=100)
t = np.array(log.time)
e = np.array(log.e_tot)
ax.plot(t, 1e6 * (e - e[0]) / N, color=COLORS[0], lw=0.8)
slope = np.polyfit(t, e, 1)[0] / N * 1e6 * 1000
ax.set_xlabel("time (ps)")
ax.set_ylabel(r"$\delta E_\mathrm{tot}$ ($\mu$eV/atom)")
ax.set_title("(b) 300 ps NVE at 2 fs")
ax.text(0.04, 0.06,
        f"drift = {slope:+.3f} $\\mu$eV/atom/ns\n"
        f"max |p| = {max(log.momentum):.1e} u A/ps",
        transform=ax.transAxes, fontsize=8)

# ---- (c) canonical sampling ---------------------------------------------- #
ax = axes[2]
T0 = 600.0
specs = [
    ("berendsen", dict(tau=0.2)),
    ("bussi", dict(tau=0.2)),
    ("langevin", dict(friction=5.0)),
    ("nose-hoover", dict(tau=0.2)),
]
print(f"{'thermostat':>13} {'<T>/K':>9} {'var(T)/T^2':>12} {'exact 2/Nf':>11}")
for k, (name, kw) in enumerate(specs):
    th = tl.make_thermostat(name, T0, **kw)
    sim = tl.Simulation(cfg, pot, timestep=0.002, thermostat=th, seed=3 + k)
    sim.set_temperature(T0)
    sim.run(20000, log_every=0)
    sim.clear_log()
    log = sim.run(200000, log_every=20)
    temps = np.array(log.temperature)
    x = temps / T0
    ax.hist(x, bins=70, density=True, histtype="step", color=COLORS[k], label=name, lw=1.1)
    print(f"{name:>13} {temps.mean():9.2f} {temps.var()/T0**2:12.5f} {2/sim.dof:11.5f}")
nf = sim.dof
xx = np.linspace(0.7, 1.3, 400)
ax.plot(xx, chi2.pdf(xx * nf, nf) * nf, "k--", lw=1.1, label=r"exact $\chi^2_{N_f}$")
ax.set_xlim(0.78, 1.22)
ax.set_xlabel(r"$T_\mathrm{inst}/T_0$")
ax.set_ylabel("probability density")
ax.set_title(f"(c) canonical test, $T_0$ = {T0:.0f} K")
ax.legend(loc="upper right")

fig.tight_layout()
fig.savefig(OUT / "fig2_validation.png")
print("wrote", OUT / "fig2_validation.png")
