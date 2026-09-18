"""Input handling: a TOML schema plus a reader for the original ``start.txt``.

The legacy format was a bare list of values whose meaning came from their line
number, decoded with a hard-coded ``l = 28``.  Adding one comment line to the
file silently shifted every parameter by one.  :func:`read_legacy_start_txt`
parses it by *key* instead, so it is robust, and it is kept so that existing
input decks still run unchanged.
"""

from __future__ import annotations

import tomllib
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .lattice import Configuration, add_defects, make_binary_mixture
from .lattices import LATTICE_NAMES, build_lattice, emptiest_points, lattice_spec
from .potentials import LennardJones, TruncationMode
from .thermostats import Thermostat, make_thermostat


@dataclass
class SystemSpec:
    lattice: str = "triangular"   #: triangular | square | honeycomb | kagome
    nx: int = 16
    ny: int = 9
    a: float = 3.2  #: A, nearest-neighbour distance
    mass: float = 24.305  #: u (magnesium)
    vacancies: int = 0
    interstitials: int = 0
    binary_fraction: float = 0.0


@dataclass
class PotentialSpec:
    epsilon: float = 0.27  #: eV
    sigma: float = 2.88  #: A
    cutoff: float = 8.0  #: A
    mode: TruncationMode = "shifted-force"
    #: strength of the three-body angular term in eV; 0 disables it.  Anything
    #: from ~0.25 eV holds the open lattices together (see trilattice.threebody).
    three_body: float = 0.0


@dataclass
class RunSpec:
    temperature: float = 300.0  #: K
    timestep: float = 0.002  #: ps
    steps: int = 20000
    equilibration: int = 5000
    thermostat: str = "bussi"
    tau: float = 0.2  #: ps
    skin: float = 1.0  #: A
    seed: int = 20240917
    log_every: int = 50
    sample_every: int = 100


@dataclass
class Settings:
    system: SystemSpec = field(default_factory=SystemSpec)
    potential: PotentialSpec = field(default_factory=PotentialSpec)
    run: RunSpec = field(default_factory=RunSpec)

    # -- construction ---------------------------------------------------- #
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Settings":
        return cls(
            system=SystemSpec(**data.get("system", {})),
            potential=PotentialSpec(**data.get("potential", {})),
            run=RunSpec(**data.get("run", {})),
        )

    @classmethod
    def from_toml(cls, path: str | Path) -> "Settings":
        with open(path, "rb") as fh:
            return cls.from_dict(tomllib.load(fh))

    @classmethod
    def load(cls, path: str | Path) -> "Settings":
        """Read either a TOML file or a legacy ``start.txt``."""
        path = Path(path)
        if path.suffix.lower() == ".toml":
            return cls.from_toml(path)
        return read_legacy_start_txt(path)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    # -- materialisation -------------------------------------------------- #
    def build_configuration(self, rng: np.random.Generator | None = None) -> Configuration:
        s = self.system
        rng = np.random.default_rng(self.run.seed) if rng is None else rng
        cfg = build_lattice(s.lattice, s.nx, s.ny, s.a, s.mass)
        if s.vacancies:
            cfg = add_defects(cfg, n_vacancies=s.vacancies, rng=rng)
        if s.interstitials:
            sites = emptiest_points(cfg.positions, cfg.box, n=s.interstitials)
            cfg = add_defects(cfg, n_interstitials=s.interstitials, a=s.a, rng=rng,
                              candidate_sites=sites, min_separation=0.0)
        if s.binary_fraction > 0.0:
            cfg = make_binary_mixture(cfg, s.binary_fraction, rng=rng)
        return cfg

    @property
    def lattice(self):
        """The :class:`~trilattice.lattices.LatticeSpec` for this deck."""
        return lattice_spec(self.system.lattice)

    def build_potential(self) -> LennardJones:
        p = self.potential
        if self.system.binary_fraction > 0.0:
            from .potentials import kob_andersen_like

            return kob_andersen_like(p.epsilon, p.sigma, p.cutoff, mode=p.mode)
        return LennardJones(p.epsilon, p.sigma, p.cutoff, mode=p.mode)

    def build_three_body(self):
        """The angular term for this deck's lattice, or None when disabled."""
        if self.potential.three_body <= 0.0:
            return None
        from .threebody import ThreeBodyAngular

        return ThreeBodyAngular.for_lattice(self.lattice, self.system.a,
                                            strength=self.potential.three_body)

    def build_thermostat(self, temperature: float | None = None) -> Thermostat:
        t = self.run.temperature if temperature is None else temperature
        kind = self.run.thermostat
        if kind in ("none", "nve"):
            return make_thermostat(kind, t)
        if kind == "langevin":
            return make_thermostat(kind, t, friction=1.0 / self.run.tau)
        return make_thermostat(kind, t, tau=self.run.tau)


# --------------------------------------------------------------------------- #
# Legacy reader
# --------------------------------------------------------------------------- #

_LEGACY_KEYS = {
    "number_x_entered": "nx",
    "number_y_entered": "ny",
    "sigma angstrom": "sigma",
    "sigma_cutoff angstrom": "cutoff",
    "a angstrom": "a",
    "m uam": "mass",
    "epsilon": "epsilon",
    "additional_particles": "interstitials",
    "vacancy": "vacancies",
    "termo": "temperature",
    "iter": "steps",
}


def read_legacy_start_txt(path: str | Path) -> Settings:
    """Parse the original ``start.txt`` (``#key`` / value line pairs)."""
    lines = [ln.strip() for ln in Path(path).read_text(encoding="utf-8").splitlines()]
    lines = [ln for ln in lines if ln]

    raw: dict[str, str] = {}
    key: str | None = None
    for ln in lines:
        if ln.startswith("#"):
            key = ln.lstrip("#").strip().lower()
        elif key is not None:
            raw[key] = ln
            key = None

    def get(name: str, default):
        for legacy, canonical in _LEGACY_KEYS.items():
            if canonical == name and legacy in raw:
                return type(default)(float(raw[legacy]))
        return default

    system = SystemSpec(
        nx=int(get("nx", 16)),
        ny=int(get("ny", 9)),
        a=get("a", 3.2),
        mass=get("mass", 24.305),
        vacancies=int(get("vacancies", 0)),
        interstitials=int(get("interstitials", 0)),
    )
    potential = PotentialSpec(
        epsilon=get("epsilon", 0.27),
        sigma=get("sigma", 2.88),
        cutoff=get("cutoff", 8.0),
    )
    run = RunSpec(
        temperature=get("temperature", 300.0),
        steps=int(get("steps", 20000)),
    )
    return Settings(system=system, potential=potential, run=run)


DEFAULT_TOML = """\
# trilattice input deck -- all quantities in metal units
# (angstrom, eV, atomic mass units, picoseconds, kelvin)

[system]
lattice = "triangular"   # triangular | square | honeycomb | kagome
                         # only "triangular" is a stable ground state of an
                         # isotropic pair potential; the others are metastable
                         # and collapse towards it when heated -- which is the
                         # point of being able to switch between them
nx = 16            # cells along x
ny = 9             # cells along y
a = 3.2            # nearest-neighbour distance, A
mass = 24.305      # u -- magnesium
vacancies = 0
interstitials = 0
binary_fraction = 0.0   # >0 turns on a size-asymmetric glass-forming mixture

[potential]
epsilon = 0.27     # eV
sigma = 2.88       # A
cutoff = 8.0       # A  (2.78 sigma)
mode = "shifted-force"   # cut | shifted | shifted-force
three_body = 0.0         # eV; >0 adds a Stillinger-Weber-style angular term that
                         # holds open lattices (square, honeycomb, kagome) together.
                         # 0.5 is a good value; 0 leaves a pure pair potential,
                         # in which only the triangular lattice is stable.

[run]
temperature = 300.0
timestep = 0.002   # ps = 2 fs
equilibration = 5000
steps = 20000
thermostat = "bussi"     # none | berendsen | bussi | langevin | nose-hoover
tau = 0.2          # ps
skin = 1.0         # A
seed = 20240917
log_every = 50
sample_every = 100
"""

__all__ = [
    "Settings",
    "SystemSpec",
    "PotentialSpec",
    "RunSpec",
    "read_legacy_start_txt",
    "DEFAULT_TOML",
    "LATTICE_NAMES",
]
