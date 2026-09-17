"""trilattice -- 2-D molecular dynamics of a triangular Lennard-Jones lattice.

Quick start::

    from trilattice import triangular_lattice, LennardJones, Simulation, Bussi

    config = triangular_lattice(nx=16, ny=9, a=3.2, mass=24.305)
    potential = LennardJones(epsilon=0.27, sigma=2.88, cutoff=8.0)
    sim = Simulation(config, potential, timestep=0.002,
                     thermostat=Bussi(temperature=300.0, tau=0.2), seed=1)
    sim.set_temperature(300.0)
    sim.run(20000)
"""

from .animate import LiveViewer, ViewerConfig
from .dashboard import Dashboard, DashboardConfig
from .analysis import Estimate, block_average, estimate, jackknife
from .config import Settings, read_legacy_start_txt
from .forces import ForceField, HAVE_NUMBA
from .lattice import (
    Box,
    Configuration,
    add_defects,
    interstitial_sites,
    make_binary_mixture,
    perfect_lattice_neighbour_shells,
    replicate,
    required_replication,
    triangular_lattice,
)
from .neighbors import NeighborList
from .observables import (
    RDFAccumulator,
    coordination_by_delaunay,
    defect_fraction,
    diffusion_coefficient,
    global_psi6,
    lindemann_2d,
    mean_squared_displacement,
    pressure_2d,
    psi6,
    psi6_correlation,
    structure_factor,
)
from .potentials import LennardJones, kob_andersen_like
from .simulation import RunLog, Simulation
from .thermostats import (
    Berendsen,
    Bussi,
    Langevin,
    NoseHooverChain,
    NoThermostat,
    make_thermostat,
)
from .trajectory import Trajectory
from .units import KB, FTM2V, MVV2E, ReducedUnits

__version__ = "1.0.0"

__all__ = [
    "__version__",
    "Box",
    "Configuration",
    "triangular_lattice",
    "add_defects",
    "interstitial_sites",
    "make_binary_mixture",
    "perfect_lattice_neighbour_shells",
    "replicate",
    "required_replication",
    "LennardJones",
    "kob_andersen_like",
    "NeighborList",
    "ForceField",
    "HAVE_NUMBA",
    "Simulation",
    "RunLog",
    "Trajectory",
    "Settings",
    "read_legacy_start_txt",
    "NoThermostat",
    "Berendsen",
    "Bussi",
    "Langevin",
    "NoseHooverChain",
    "make_thermostat",
    "RDFAccumulator",
    "psi6",
    "global_psi6",
    "psi6_correlation",
    "structure_factor",
    "coordination_by_delaunay",
    "defect_fraction",
    "mean_squared_displacement",
    "diffusion_coefficient",
    "lindemann_2d",
    "pressure_2d",
    "estimate",
    "block_average",
    "jackknife",
    "Estimate",
    "LiveViewer",
    "ViewerConfig",
    "Dashboard",
    "DashboardConfig",
    "KB",
    "MVV2E",
    "FTM2V",
    "ReducedUnits",
]
