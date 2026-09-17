# trilattice

2-D molecular dynamics of a Lennard-Jones triangular lattice — simulate, watch, and measure the melting of a 2-D crystal in real time.

![dashboard](figures/dashboard.gif)

## Install

```bash
git clone <this repo> && cd trilattice
pip install -e ".[dev,fast]"      # fast = numba (5x speedup), dev = pytest
```

Requires Python ≥ 3.11, NumPy, SciPy, Matplotlib. Numba is optional.

## Quick start

```bash
trilattice init -o config.toml          # write a commented input deck
trilattice animate config.toml -T 600   # open the live dashboard
```

That's it — drag the temperature slider up and watch the crystal melt.

## Command line

| command | what it does |
|---|---|
| `trilattice init` | write a commented `config.toml` |
| `trilattice info config.toml` | derived quantities (density, timescales, sanity checks) |
| `trilattice run config.toml -T 300 -n 20000` | run one simulation, print statistics |
| `trilattice scan config.toml --t-min 200 --t-max 3000` | sweep temperature, print a table |
| `trilattice animate config.toml -T 600` | **interactive dashboard** (below) |
| `trilattice animate config.toml --minimal` | small keyboard-driven viewer |
| `trilattice bench` | throughput vs. system size |

Add `--xyz traj.xyz` to `run` to export a trajectory for OVITO/VMD, with per-atom hexatic order and coordination as extra columns.

## The dashboard

```bash
trilattice animate examples/config.toml -T 600
```

| control | effect |
|---|---|
| **target T** slider | drives the thermostat, live |
| **rate** slider | simulation speed — down to slow motion (fractions of a step per frame), *never* changes the timestep |
| **time** slider | scrub back through the last 500 frames |
| `resume here` | restart the dynamics from the scrubbed frame, discarding the future |
| `quench` / `reheat` | relax to the nearest local minimum / redraw thermal velocities |
| `+vacancy` / `+interstitial` | inject a point defect into the running crystal |
| `reset` | back to the starting configuration |
| colour: `psi6` `coord` `speed` `energy` `displacement` | what the atoms are shaded by |
| overlays: bonds, trails, box, thermostat on/off | toggle on the fly |
| analysis panel: `g(r)` `S(k)` `MSD` `speeds` `coordination` `T–ψ6` | live measurement, updating as you watch |

Headless machines get a compact recorder instead of a window:

```bash
trilattice animate config.toml -T 2800 --save melting.gif --frames 150
```

## Python API

```python
import trilattice as tl

config    = tl.triangular_lattice(nx=20, ny=12, a=3.2, mass=24.305)
potential = tl.LennardJones(epsilon=0.27, sigma=2.88, cutoff=8.0)
sim = tl.Simulation(config, potential, timestep=0.002,
                    thermostat=tl.Bussi(temperature=1500.0, tau=0.2), seed=1)
sim.set_temperature(1500.0)
sim.run(10000, log_every=0)          # equilibrate
sim.clear_log(); sim.reset_origin()
log = sim.run(40000, log_every=50)   # produce

from trilattice import observables as obs
print(obs.global_psi6(sim.positions, sim.box, 1.35 * 3.2))   # hexatic order
print(obs.defect_fraction(sim.positions, sim.box))           # % not six-coordinated
```

## What's measured

![melting](figures/fig3_melting.png)

| quantity | function |
|---|---|
| hexatic order \|⟨ψ₆⟩\| | `observables.global_psi6` |
| topological defects (Voronoi ≠ 6) | `observables.defect_fraction` |
| pair correlation g(r) | `observables.RDFAccumulator` |
| structure factor S(k) | `observables.structure_factor` |
| self-diffusion D | `observables.diffusion_coefficient` |
| Lindemann parameter | `observables.lindemann_2d` |
| pressure, heat capacity | `observables.pressure_2d`, `heat_capacity_nvt` |

## Package layout

| module | responsible for |
|---|---|
| `lattice.py` | building the triangular lattice, periodic box, vacancies/interstitials |
| `potentials.py` | the Lennard-Jones pair potential (three truncation schemes) |
| `neighbors.py` | O(N) neighbour search (cell list + Verlet list) |
| `forces.py` | force/energy evaluation (NumPy and Numba kernels) |
| `simulation.py` | the velocity-Verlet integrator and FIRE relaxation |
| `thermostats.py` | Berendsen, Bussi, Langevin, Nosé–Hoover |
| `observables.py` | every physical quantity listed above |
| `dashboard.py` / `animate.py` | the interactive viewers |
| `config.py` | reading/writing `config.toml` |
| `cli.py` | the `trilattice` command |

## Tests

```bash
pytest -q
```