"""Two-dimensional melting of the triangular lattice: the production run.

Sweeps temperature up from the perfect crystal and back down from the melt,
which brackets the transition and exposes its hysteresis.  At every temperature
it measures

  * potential energy per atom and, from its fluctuations, the heat capacity,
  * the hexatic order parameter |<psi_6>|,
  * the fraction of particles whose Voronoi coordination is not six
    (i.e. the density of topological defects),
  * the modified Lindemann parameter,
  * the self-diffusion constant from the Einstein relation,
  * g(r) and S(k) at a few selected temperatures.

Results are written to ``data/melting.npz`` so that the figures can be redrawn
without repeating the simulation.
"""
from pathlib import Path
import sys
import time
import numpy as np

import trilattice as tl
from trilattice import observables as obs
from trilattice.analysis import estimate

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

# ---------------------------------------------------------------- settings --
NX, NY = 20, 12
RHO_STAR = 0.85                 # reduced density: the 2-D LJ triple point sits near here
A0 = 2.88 * (2.0 / (RHO_STAR * 3 ** 0.5)) ** 0.5   # = 3.3567 A
EPS, SIG, RC = 0.27, 2.88, 8.0
DT = 0.002                      # ps
EQUIL, PROD = 8000, 20000       # steps  (16 ps + 40 ps)
SAMPLE, DEFECT_EVERY = 25, 250
PSI6_CUT = 1.35 * A0
TEMPS = np.arange(200.0, 3001.0, 200.0)
SNAPSHOT_AT = [400.0, 1400.0, 1800.0, 2600.0]
BRANCH = sys.argv[1] if len(sys.argv) > 1 else "heat"

pot = tl.LennardJones(EPS, SIG, RC, mode="shifted-force")
base = tl.triangular_lattice(NX, NY, A0)
N = base.n_particles
print(f"N = {N}, box = {base.box}, rho = {base.density:.5f} A^-2")
print(f"equil {EQUIL} + prod {PROD} steps at dt = {DT*1000:.0f} fs "
      f"({(EQUIL+PROD)*DT:.0f} ps) per temperature, {len(TEMPS)} temperatures, 2 branches\n")


def sweep(temperatures, seed, start_config=None, label=""):
    rows = []
    snaps = {}
    sim = None
    t_wall = time.perf_counter()
    for k, T in enumerate(temperatures):
        th = tl.Bussi(temperature=float(T), tau=1.0)
        cfg = base if start_config is None else start_config
        if sim is None:
            sim = tl.Simulation(cfg, pot, timestep=DT, thermostat=th, seed=seed, skin=1.2)
            sim.set_temperature(float(T))
        else:
            sim.thermostat = th
        sim.run(EQUIL, log_every=0)
        sim.clear_log()
        sim.reset_origin()

        traj, p6, defects = [], [], []
        def sampler(s, _k=[0]):
            _k[0] += 1
            traj.append(s.unwrapped.copy())
            p6.append(obs.global_psi6(s.positions, s.box, PSI6_CUT))
            if _k[0] % (DEFECT_EVERY // SAMPLE) == 0:
                defects.append(obs.defect_fraction(s.positions, s.box))

        log = sim.run(PROD, log_every=20, callback=sampler, callback_every=SAMPLE)

        traj = np.asarray(traj)
        times = np.arange(traj.shape[0]) * SAMPLE * DT
        msd = obs.mean_squared_displacement(traj)
        d_coef = obs.diffusion_coefficient(times, msd, fit_from=0.5)

        e_pot = np.array(log.e_pot)
        e_tot = np.array(log.e_tot)
        t_mean = float(np.mean(log.temperature))
        cv = obs.heat_capacity_nvt(e_tot, t_mean, N)

        rows.append(dict(
            T_set=float(T), T=t_mean,
            e_pot=estimate(e_pot / N).mean, e_pot_err=estimate(e_pot / N).error,
            pressure=float(np.mean(log.pressure)),
            psi6=float(np.mean(p6)), psi6_err=float(np.std(p6) / np.sqrt(max(1, len(p6)))),
            cv=cv,
            defects=float(np.mean(defects)) if defects else np.nan,
            lindemann=obs.lindemann_2d(traj, sim.box, A0),
            diffusion=d_coef,
            msd=msd, msd_time=times,
        ))
        if float(T) in SNAPSHOT_AT and label == "heat":
            rdf = obs.RDFAccumulator(sim.box, r_max=18.0, n_bins=360)
            for fr in traj[::4]:
                rdf.accumulate(sim.box.wrap(fr))
            kx, ky, sk = obs.structure_factor(sim.positions, sim.box, n_max=22)
            snaps[float(T)] = dict(
                positions=sim.positions.copy(),
                psi6=obs.psi6(sim.positions, sim.box, PSI6_CUT),
                coord=obs.coordination_by_delaunay(sim.positions, sim.box),
                rdf_r=rdf.result()[0], rdf_g=rdf.result()[1],
                kx=kx, ky=ky, sk=sk,
            )
        print(f"  [{label}] T = {T:6.0f} K -> {t_mean:7.1f}  psi6 = {rows[-1]['psi6']:.3f}  "
              f"E = {rows[-1]['e_pot']:+.4f} eV/at  def = {rows[-1]['defects']:.3f}  "
              f"D = {d_coef:8.4f} A^2/ps", flush=True)
    print(f"  [{label}] branch wall time {time.perf_counter()-t_wall:.0f} s\n")
    return rows, snaps, sim


if BRANCH == "heat":
    rows, snaps, sim_end = sweep(TEMPS, seed=101, label="heat")
    np.savez_compressed(DATA / "molten_seed.npz", positions=sim_end.positions,
                        masses=sim_end.masses, box=np.array([sim_end.box.lx, sim_end.box.ly]))
else:
    d = np.load(DATA / "molten_seed.npz")
    molten = tl.Configuration(d["positions"], tl.Box(*d["box"]), d["masses"])
    rows, snaps, _ = sweep(TEMPS[::-1], seed=202, start_config=molten, label="cool")


def pack(rows, keys):
    return {k: np.array([r[k] for r in rows]) for k in keys}


scalar_keys = ["T_set", "T", "e_pot", "e_pot_err", "pressure", "psi6", "psi6_err",
               "cv", "defects", "lindemann", "diffusion"]
out = {f"{BRANCH}_{k}": v for k, v in pack(rows, scalar_keys).items()}
out["msd_time"] = rows[0]["msd_time"]
out[f"{BRANCH}_msd"] = np.array([r["msd"] for r in rows])
out["box"] = np.array([base.box.lx, base.box.ly])
out["N"] = N
out["a0"] = A0
for T, dd in snaps.items():
    for key, val in dd.items():
        out[f"snap_{int(T)}_{key}"] = val

np.savez_compressed(DATA / f"melting_{BRANCH}.npz", **out)
print("wrote", DATA / f"melting_{BRANCH}.npz")
