"""Figures 3-5: melting diagnostics, structure, and dynamics.

Reads ``data/melting_{heat,cool}.npz`` produced by ``02_melting.py``; no
simulation is repeated here.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import trilattice as tl
from trilattice.plotting import use_style, COLORS, plot_configuration, annotate_shells

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
use_style()

h = np.load(ROOT / "data" / "melting_heat.npz")
c = np.load(ROOT / "data" / "melting_cool.npz")
N = int(h["N"])
A0 = float(h["a0"])
box = tl.Box(*h["box"])
U = tl.ReducedUnits(0.27, 2.88, 24.305)
TSTAR = lambda t: U.temperature(t)

# =========================================================== figure 3 ===== #
fig, ax = plt.subplots(2, 3, figsize=(11.5, 6.4))
BR = [("heat", h, COLORS[1], "o-", "heating from the crystal"),
      ("cool", c, COLORS[0], "s--", "cooling from the melt")]


def top_axis(a):
    sec = a.secondary_xaxis("top", functions=(TSTAR, U.to_kelvin))
    sec.set_xlabel(r"$T^* = k_BT/\varepsilon$", fontsize=8)
    sec.tick_params(labelsize=7)


# (a) hexatic order parameter
for name, d, col, st, lab in BR:
    ax[0, 0].errorbar(d[f"{name}_T"], d[f"{name}_psi6"], yerr=d[f"{name}_psi6_err"],
                      fmt=st, color=col, ms=4, capsize=2, label=lab)
ax[0, 0].set_ylabel(r"$|\langle\psi_6\rangle|$")
ax[0, 0].set_title("(a) hexatic order")
ax[0, 0].legend(loc="lower left")

# (b) potential energy + the latent-heat-like step
for name, d, col, st, lab in BR:
    ax[0, 1].plot(d[f"{name}_T"], d[f"{name}_e_pot"], st, color=col, ms=4)
ax[0, 1].set_ylabel(r"$E_\mathrm{pot}$ (eV/atom)")
ax[0, 1].set_title("(b) potential energy")

# (c) heat capacity from energy fluctuations
for name, d, col, st, lab in BR:
    ax[0, 2].plot(d[f"{name}_T"], d[f"{name}_cv"], st, color=col, ms=4)
ax[0, 2].axhline(2.0, color="0.5", lw=0.8, ls=":")
ax[0, 2].text(250, 2.06, r"Dulong-Petit, 2-D: $C_v = 2k_B$", fontsize=7, color="0.4")
ax[0, 2].set_ylabel(r"$C_v$ ($k_B$/atom)")
ax[0, 2].set_title("(c) heat capacity")

# (d) topological defects
for name, d, col, st, lab in BR:
    ax[1, 0].plot(d[f"{name}_T"], 100 * d[f"{name}_defects"], st, color=col, ms=4)
ax[1, 0].set_ylabel("non-six-coordinated (%)")
ax[1, 0].set_title("(d) topological defects")

# (e) diffusion
for name, d, col, st, lab in BR:
    y = np.maximum(d[f"{name}_diffusion"], 1e-4)
    ax[1, 1].semilogy(d[f"{name}_T"], y, st, color=col, ms=4)
ax[1, 1].set_ylabel(r"$D$ (A$^2$/ps)")
ax[1, 1].set_title("(e) self-diffusion")
ax[1, 1].set_ylim(5e-5, 5)

# (f) pressure along the isochore
for name, d, col, st, lab in BR:
    ax[1, 2].plot(d[f"{name}_T"], d[f"{name}_pressure"], st, color=col, ms=4)
ax[1, 2].axhline(0.0, color="0.5", lw=0.8, ls=":")
ax[1, 2].set_ylabel(r"$P$ (eV/A$^2$)")
ax[1, 2].set_title("(f) pressure")

# melting window from the steepest drop of psi6 on the heating branch
dpsi = np.gradient(h["heat_psi6"], h["heat_T"])
tm = float(h["heat_T"][np.argmin(dpsi)])
for a in ax.ravel():
    a.axvspan(2200, 2450, color="0.85", zorder=0)
    a.set_xlabel("$T$ (K)")
    top_axis(a)
ax[0, 0].annotate(f"$T_m \\approx$ {tm:.0f} K\n($T^*$ = {TSTAR(tm):.2f})",
                  xy=(tm, 0.45), xytext=(1150, 0.32), fontsize=8,
                  arrowprops=dict(arrowstyle="->", lw=0.7, color="0.4"))
fig.suptitle(
    f"2-D Lennard-Jones triangular lattice, N = {N}, "
    r"$\rho^*$ = " + f"{U.density(N / box.area):.2f}, 56 ps per temperature",
    fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT / "fig3_melting.png")
print("wrote fig3_melting.png   T_m ~", round(tm), "K")

# =========================================================== figure 4 ===== #
snaps = [400, 1400, 1800, 2600]
fig, ax = plt.subplots(3, 4, figsize=(12.5, 8.6))
radii, counts = tl.perfect_lattice_neighbour_shells(A0, 17.0)
for k, T in enumerate(snaps):
    pos = h[f"snap_{T}_positions"]
    p6 = np.abs(h[f"snap_{T}_psi6"])
    coord = h[f"snap_{T}_coord"]

    plot_configuration(ax[0, k], pos, box, color=p6, cmap="viridis", vmin=0, vmax=1,
                       size=13, title=f"T = {T} K  ($T^*$ = {TSTAR(T):.2f})",
                       cbar_label=r"$|\psi_6|$" if k == 3 else None)

    # colour by Voronoi coordination: 5 red, 6 grey, 7 blue
    col = np.full(len(coord), "0.78", dtype=object)
    col[coord == 5] = COLORS[1]
    col[coord == 7] = COLORS[0]
    col[(coord != 5) & (coord != 6) & (coord != 7)] = COLORS[3]
    ax[1, k].scatter(pos[:, 0], pos[:, 1], c=list(col), s=13, linewidths=0.3, edgecolors="0.4")
    ax[1, k].add_patch(plt.Rectangle((0, 0), box.lx, box.ly, fill=False, ec="0.4", lw=0.8, ls="--"))
    ax[1, k].set_aspect("equal"); ax[1, k].grid(False)
    ax[1, k].set_xlim(-2, box.lx + 2); ax[1, k].set_ylim(-2, box.ly + 2)
    ax[1, k].set_xlabel("x (A)")
    ax[1, k].set_title(f"5 / 7 defects: {100 * np.mean(coord != 6):.1f}%", fontsize=9)

    ax[2, k].plot(h[f"snap_{T}_rdf_r"], h[f"snap_{T}_rdf_g"], color=COLORS[0], lw=1.1)
    ax[2, k].set_xlim(2, 17); ax[2, k].set_ylim(0, 5.2)
    ax[2, k].set_xlabel("$r$ (A)")
    if k == 0:
        annotate_shells(ax[2, k], radii, counts, y=5.0, max_shells=7)
ax[0, 0].set_ylabel("y (A)")
ax[1, 0].set_ylabel("y (A)")
ax[2, 0].set_ylabel("$g(r)$")
for a in ax[1, 1:]:
    a.set_ylabel("")
fig.suptitle("configurations coloured by hexatic order (top) and Voronoi coordination "
             "(middle); pair correlation (bottom)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT / "fig4_structure.png")
print("wrote fig4_structure.png")

# =========================================================== figure 5 ===== #
fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.5))

# (a) structure factor: Bragg spots -> ring
for k, T in enumerate([400, 2600]):
    a = ax[k]
    kx, ky, s = h[f"snap_{T}_kx"], h[f"snap_{T}_ky"], h[f"snap_{T}_sk"]
    im = a.pcolormesh(kx, ky, np.clip(s.T, 0.05, None), norm=LogNorm(vmin=0.05, vmax=N),
                      cmap="magma", shading="auto")
    a.set_aspect("equal"); a.grid(False)
    a.set_xlabel(r"$k_x$ (A$^{-1}$)"); a.set_ylabel(r"$k_y$ (A$^{-1}$)")
    a.set_title(f"(a{k + 1}) $S(\\mathbf{{k}})$, T = {T} K")
    fig.colorbar(im, ax=a, fraction=0.046, pad=0.03)

# (b) MSD
a = ax[2]
t = h["msd_time"]
temps = h["heat_T_set"]
sel = [1, 4, 7, 10, 11, 12, 14]
for k, idx in enumerate(sel):
    a.loglog(t[1:], h["heat_msd"][idx][1:], color=plt.cm.inferno(k / len(sel)),
             lw=1.2, label=f"{temps[idx]:.0f} K")
a.loglog(t[1:], 1.2 * t[1:], "k:", lw=0.9)
a.text(3.0, 6.0, r"$\propto t$", fontsize=8)
a.axhline(A0**2, color="0.6", lw=0.7, ls="--")
a.text(0.12, A0**2 * 1.15, r"$a^2$", fontsize=7, color="0.4")
a.set_xlabel("$t$ (ps)"); a.set_ylabel(r"MSD (A$^2$)")
a.set_title("(b) mean squared displacement")
a.legend(fontsize=6.5, ncol=2, loc="upper left")
fig.tight_layout()
fig.savefig(OUT / "fig5_dynamics.png")
print("wrote fig5_dynamics.png")

# ---- numbers for the README --------------------------------------------- #
print()
print("T_m (heating, steepest psi6 drop) =", round(tm), "K, T* =", round(TSTAR(tm), 3))
i_hot = np.argmin(np.abs(h["heat_T"] - 2400))
print("D jump across the transition: "
      f"{h['heat_diffusion'][i_hot - 1]:.4f} -> {h['heat_diffusion'][i_hot]:.4f} A^2/ps")
print("low-T Cv (heating, T<600K) =", np.round(h["heat_cv"][:3], 3), "k_B/atom")
print("residual defects after the quench (200 K, cooling) =",
      round(100 * float(c["cool_defects"][-1]), 1), "%")
