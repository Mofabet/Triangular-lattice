"""Unit system for the simulation.

The code works throughout in *metal* units (the LAMMPS convention), because that
is what the input parameters of the original project were naturally expressed in
(angstroms, electronvolts, atomic mass units) and because it avoids the
catastrophic dynamic range of SI floating point arithmetic in this problem:

===========  ================================
quantity     unit
===========  ================================
length       angstrom (A)
energy       electronvolt (eV)
mass         unified atomic mass unit (u)
time         picosecond (ps)
temperature  kelvin (K)
velocity     A/ps
force        eV/A
pressure     eV/A^2   (2-D: a force per unit length, i.e. a surface tension)
===========  ================================

The only non-trivial conversion is the one that ties mass, length and time to
energy.  With ``E = 1/2 m v^2``::

    1 u * (A/ps)^2 = 1.66053906660e-27 kg * 1e-20 m^2 / 1e-24 s^2
                   = 1.66053906660e-23 J
                   = 1.0364269...e-4 eV

so kinetic energies must be multiplied by :data:`MVV2E`, and Newton's second law
reads ``a[A/ps^2] = FTM2V * F[eV/A] / m[u]``.

Getting this factor wrong is the single most common bug in a hand-written MD
code: the trajectory still *looks* plausible, but the temperature, the diffusion
constant and the melting point are all off by orders of magnitude.

All constants are CODATA-2018 / SI-2019 exact values where such a value exists.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# --------------------------------------------------------------------------- #
# Fundamental constants (SI)
# --------------------------------------------------------------------------- #

#: Boltzmann constant, exact by SI definition [J/K].
KB_SI = 1.380649e-23
#: Elementary charge, exact by SI definition [C]; numerically equal to 1 eV in J.
EV_SI = 1.602176634e-19
#: Unified atomic mass unit [kg] (CODATA 2018).
AMU_SI = 1.66053906660e-27
#: Avogadro constant, exact by SI definition [1/mol].
NA = 6.02214076e23
#: Planck constant, exact by SI definition [J s].
H_SI = 6.62607015e-34

# --------------------------------------------------------------------------- #
# Metal units
# --------------------------------------------------------------------------- #

#: Boltzmann constant in metal units [eV/K].
KB = KB_SI / EV_SI  # 8.617333262e-05

#: Conversion from ``u * (A/ps)^2`` to eV.  Multiply ``1/2 m v^2`` by this.
MVV2E = AMU_SI * 1.0e-20 / 1.0e-24 / EV_SI  # 1.0364269e-04

#: Conversion from ``eV/A / u`` to ``A/ps^2``.  ``a = FTM2V * f / m``.
FTM2V = 1.0 / MVV2E  # 9648.533...

#: 1 picosecond in femtoseconds (convenience).
FS_PER_PS = 1000.0

#: Angstrom in metres, for reporting only.
ANGSTROM_SI = 1.0e-10


def kelvin_to_ev(temperature_k: float) -> float:
    """k_B * T in eV."""
    return KB * temperature_k


def ev_to_kelvin(energy_ev: float) -> float:
    """Temperature in K whose thermal energy k_B*T equals ``energy_ev``."""
    return energy_ev / KB


def thermal_velocity(temperature_k: float, mass_u: float) -> float:
    """Most probable speed sqrt(2 k_B T / m) of the 2-D Maxwell distribution [A/ps]."""
    return math.sqrt(2.0 * KB * temperature_k / (mass_u * MVV2E))


def einstein_timescale(mass_u: float, epsilon_ev: float, sigma_a: float) -> float:
    """Natural Lennard-Jones time unit tau = sigma * sqrt(m/epsilon) [ps].

    A velocity-Verlet step of ``dt ~ tau/200`` conserves energy to ~1e-6 per
    nanosecond for this potential; it is the right yardstick for choosing ``dt``.
    """
    return sigma_a * math.sqrt(mass_u * MVV2E / epsilon_ev)


# --------------------------------------------------------------------------- #
# Reduced (Lennard-Jones) units
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ReducedUnits:
    """Translator between metal units and Lennard-Jones reduced units.

    Reduced units are what the 2-D melting literature quotes, so every physical
    result produced by this package is reported in both systems.

    ``T* = k_B T / eps``, ``rho* = rho sigma^2``, ``E* = E / eps``,
    ``t* = t / (sigma sqrt(m/eps))``, ``P* = P sigma^2 / eps``,
    ``D* = D / (sigma^2 / tau)``.
    """

    epsilon: float  #: eV
    sigma: float  #: A
    mass: float  #: u

    @property
    def tau(self) -> float:
        """LJ time unit in ps."""
        return einstein_timescale(self.mass, self.epsilon, self.sigma)

    # --- forward (physical -> reduced) ------------------------------------- #
    def temperature(self, t_kelvin: float) -> float:
        return KB * t_kelvin / self.epsilon

    def energy(self, e_ev: float) -> float:
        return e_ev / self.epsilon

    def length(self, r_angstrom: float) -> float:
        return r_angstrom / self.sigma

    def density(self, rho_per_a2: float) -> float:
        return rho_per_a2 * self.sigma**2

    def time(self, t_ps: float) -> float:
        return t_ps / self.tau

    def pressure(self, p_ev_per_a2: float) -> float:
        return p_ev_per_a2 * self.sigma**2 / self.epsilon

    def diffusion(self, d_a2_per_ps: float) -> float:
        return d_a2_per_ps * self.tau / self.sigma**2

    # --- backward (reduced -> physical) ------------------------------------ #
    def to_kelvin(self, t_star: float) -> float:
        return t_star * self.epsilon / KB

    def to_ev(self, e_star: float) -> float:
        return e_star * self.epsilon

    def to_angstrom(self, r_star: float) -> float:
        return r_star * self.sigma

    def to_ps(self, t_star: float) -> float:
        return t_star * self.tau


__all__ = [
    "KB",
    "KB_SI",
    "EV_SI",
    "AMU_SI",
    "NA",
    "H_SI",
    "MVV2E",
    "FTM2V",
    "FS_PER_PS",
    "ANGSTROM_SI",
    "ReducedUnits",
    "kelvin_to_ev",
    "ev_to_kelvin",
    "thermal_velocity",
    "einstein_timescale",
]
