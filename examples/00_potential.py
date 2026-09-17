"""Figure 1 -- the potential itself, and the T = 0 lattice constant.

Left:   u(r) and f(r) under the three truncation schemes.
Middle: the static lattice energy per atom as a function of the lattice
        constant, whose minimum fixes a_0 and whose curvature gives the 2-D bulk
        modulus.
Right:  what the original code actually computed, for comparison.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

import trilattice as tl
from trilattice.plotting import use_style, COLORS
from trilattice.potentials import R_MIN_OVER_SIGMA

OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
use_style()

EPS, SIG, RC = 0.27, 2.88, 8.0
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))

# ---- (a) truncation schemes ---------------------------------------------- #
r = np.linspace(2.5, 9.0, 1600)
ax = axes[0]
pots = {m: tl.LennardJones(EPS, SIG, RC, mode=m) for m in ("cut", "shifted", "shifted-force")}
for k, (mode, p) in enumerate(pots.items()):
    ax.plot(r, p.energy(r), color=COLORS[k], label=mode)
ax.axvline(RC, color="0.6", lw=0.7, ls=":")
ax.axhline(0, color="0.6", lw=0.7)
ax.text(RC, 0.06, " $r_c$", fontsize=8, color="0.4")
ax.set_xlim(2.6, 9.0)
ax.set_ylim(-0.31, 0.12)
ax.set_xlabel(r"$r$ (A)")
ax.set_ylabel(r"$u(r)$ (eV)")
ax.set_title("(a) Lennard-Jones, three truncations")
ax.legend(loc="lower right", title="truncation")

# inset: the behaviour that actually differs, right at the cutoff
ins = ax.inset_axes([0.46, 0.50, 0.50, 0.44])
rz = np.linspace(7.2, 8.4, 800)
for k, (mode, p) in enumerate(pots.items()):
    ins.plot(rz, 1e3 * p.energy(rz), color=COLORS[k], lw=1.2)
    ins.plot(rz, 1e3 * p.force(rz), color=COLORS[k], lw=1.0, ls="--")
ins.axvline(RC, color="0.6", lw=0.7, ls=":")
ins.axhline(0, color="0.6", lw=0.6)
ins.set_xlim(7.2, 8.4)
ins.tick_params(labelsize=6)
ins.set_title(r"$10^3\,u$ (solid), $10^3 f$ (dashed)", fontsize=6.5)
ins.grid(alpha=0.2)

# ---- (b) cohesive curve --------------------------------------------------- #
ax = axes[1]
p = tl.LennardJones(EPS, SIG, RC, mode="cut")
a_grid = np.linspace(2.7, 4.2, 400)
e_lat = np.array([p.lattice_energy(a) for a in a_grid])
res = minimize_scalar(lambda a: p.lattice_energy(float(a)), bracket=(3.0, 3.2, 3.6))
a0, e0 = float(res.x), float(res.fun)

# 2-D bulk modulus  B = A d2E/dA2 ;  A = a^2 sqrt(3)/2 per atom
h = 1e-3
d2 = (p.lattice_energy(a0 + h) - 2 * e0 + p.lattice_energy(a0 - h)) / h**2
area_per_atom = a0**2 * np.sqrt(3) / 2
bulk = d2 * a0**2 / (4 * area_per_atom)  # eV/A^2

ax.plot(a_grid, e_lat, color=COLORS[0])
ax.plot([a0], [e0], "o", color=COLORS[1], ms=5)
ax.axvline(R_MIN_OVER_SIGMA * SIG, color="0.6", ls=":", lw=0.8)
ax.annotate(
    f"$a_0$ = {a0:.4f} A\n$E_0$ = {e0:.4f} eV/atom\n$B_{{2D}}$ = {bulk:.3f} eV/A$^2$",
    xy=(a0, e0), xytext=(3.45, e0 + 0.25), fontsize=8,
    arrowprops=dict(arrowstyle="->", lw=0.7, color="0.4"),
)
ax.text(R_MIN_OVER_SIGMA * SIG, e_lat.max() * 0.15,
        r"  $2^{1/6}\sigma$" + "\n  (dimer minimum)", fontsize=7, color="0.45")
ax.set_xlabel(r"lattice constant $a$ (A)")
ax.set_ylabel(r"$E_\mathrm{lattice}$ (eV/atom)")
ax.set_title("(b) static lattice sum")
ax.set_ylim(e0 - 0.15, e0 + 1.4)

# ---- (c) the original implementation -------------------------------------- #
ax = axes[2]
r = np.linspace(2.4, 9.0, 1200)
u_ok = 4 * EPS * ((SIG / r) ** 12 - (SIG / r) ** 6)
f_ok = 24 * EPS / r * (2 * (SIG / r) ** 12 - (SIG / r) ** 6)
u_bad = np.abs(u_ok)                                    # abs() in the original
f_bad = np.abs(24 / SIG * EPS * ((SIG / r) ** 13 - (SIG / r) ** 7))  # no 2, wrong prefactor
ax.plot(r, u_ok, color=COLORS[0], label="correct $u(r)$")
ax.plot(r, u_bad, color=COLORS[1], ls="--", label="original: $|u(r)|$")
ax.plot(r, f_ok, color=COLORS[2], label="correct $f(r)$")
ax.plot(r, f_bad, color=COLORS[3], ls="--", label="original $f(r)$")
ax.axhline(0, color="0.6", lw=0.7)
ax.set_xlim(2.4, 9.0)
ax.set_ylim(-0.4, 0.9)
ax.set_xlabel(r"$r$ (A)")
ax.set_ylabel("eV, eV/A")
ax.set_title("(c) what the original computed")
ax.legend(loc="upper right")

fig.tight_layout()
fig.savefig(OUT / "fig1_potential.png")
print(f"a0 = {a0:.6f} A   E0 = {e0:.6f} eV/atom   B_2D = {bulk:.4f} eV/A^2")
print(f"a0/sigma = {a0/SIG:.4f}  (dimer minimum is at {R_MIN_OVER_SIGMA:.4f})")
print(f"input deck used a = 3.2 A, i.e. {(3.2/a0-1)*100:+.2f}% strain")
print("wrote", OUT / "fig1_potential.png")
