"""Thermostats.

All of them are applied as a half-step operator on either side of the
velocity-Verlet core, i.e. the propagator is factorised as

    exp(iL dt) ~ T(dt/2) . B(dt/2) . A(dt) . B(dt/2) . T(dt/2)

which is the standard Trotter scheme for Nose-Hoover chains and makes the
Langevin case coincide with the BAOAB splitting (the one with the smallest
configurational sampling error of the common Langevin integrators).

Correctness hierarchy, worth knowing before trusting a number:

``berendsen``
    Not canonical.  It rescales velocities deterministically towards the target,
    which suppresses kinetic-energy fluctuations and therefore gives a *wrong*
    heat capacity; it also suffers from the "flying ice cube" artefact, draining
    energy out of internal modes into centre-of-mass motion.  Kept because it is
    what the original code used, it is robust during equilibration, and it is
    still the fastest way to bring a system to temperature.
``bussi``
    Bussi-Donadio-Parrinello stochastic velocity rescaling.  Same simplicity as
    Berendsen, but with the correct noise term, so it samples the canonical
    ensemble exactly and has a conserved quantity.  Default for production.
``langevin``
    Ornstein-Uhlenbeck friction + noise per degree of freedom.  Canonical,
    destroys hydrodynamics (so do not measure viscosity with it), but excellent
    for equilibrating and for escaping metastable states.
``nose-hoover``
    Deterministic, time-reversible Nose-Hoover chain.  Canonical provided the
    dynamics is ergodic, which the chain of length >= 2 essentially guarantees.
``none``
    Pure NVE.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .units import KB, MVV2E


def kinetic_energy(velocities: np.ndarray, masses: np.ndarray) -> float:
    """Kinetic energy in eV."""
    return float(MVV2E * 0.5 * np.sum(masses * np.sum(velocities**2, axis=1)))


def temperature_from_ke(ke: float, dof: int) -> float:
    """Instantaneous temperature in K from the kinetic energy."""
    return 2.0 * ke / (dof * KB) if dof > 0 else 0.0


class Thermostat:
    """Base class.  ``half_step`` mutates velocities in place."""

    name = "none"

    def half_step(self, velocities, masses, dt, dof, rng) -> None:
        return None

    @property
    def conserved_offset(self) -> float:
        """Energy pumped in/out, so that ``E_tot + offset`` is a constant."""
        return 0.0

    def reset(self) -> None:
        return None


class NoThermostat(Thermostat):
    name = "none"


@dataclass
class Berendsen(Thermostat):
    temperature: float
    tau: float = 0.1  #: ps
    max_scale: float = 1.25
    name: str = field(default="berendsen", init=False)

    _drained: float = field(default=0.0, init=False)

    def half_step(self, velocities, masses, dt, dof, rng) -> None:
        ke = kinetic_energy(velocities, masses)
        if ke <= 0.0:
            # Cold start: seed the system rather than divide by zero.
            velocities += 1e-8
            return
        t_now = temperature_from_ke(ke, dof)
        lam = np.sqrt(1.0 + (0.5 * dt / self.tau) * (self.temperature / t_now - 1.0))
        lam = float(np.clip(lam, 1.0 / self.max_scale, self.max_scale))
        velocities *= lam
        self._drained += ke * (lam**2 - 1.0)

    @property
    def conserved_offset(self) -> float:
        return -self._drained

    def reset(self) -> None:
        self._drained = 0.0


@dataclass
class Bussi(Thermostat):
    """Stochastic velocity rescaling, J. Chem. Phys. 126, 014101 (2007)."""

    temperature: float
    tau: float = 0.1  #: ps
    name: str = field(default="bussi", init=False)

    _drained: float = field(default=0.0, init=False)

    def half_step(self, velocities, masses, dt, dof, rng) -> None:
        ke = kinetic_energy(velocities, masses)
        if ke <= 0.0:
            sigma = np.sqrt(KB * self.temperature / (masses[:, None] * MVV2E))
            velocities += sigma * rng.standard_normal(velocities.shape)
            return
        ke_target = 0.5 * dof * KB * self.temperature
        c = np.exp(-0.5 * dt / self.tau)
        r1 = rng.standard_normal()
        # sum of (dof-1) squared standard normals ~ chi^2_{dof-1} ~ 2*Gamma((dof-1)/2)
        s = 2.0 * rng.gamma(0.5 * (dof - 1)) if dof > 1 else 0.0
        ke_new = (
            ke
            + (1.0 - c) * (ke_target * (r1 * r1 + s) / dof - ke)
            + 2.0 * r1 * np.sqrt(ke * ke_target / dof * (1.0 - c) * c)
        )
        ke_new = max(ke_new, 1e-30)
        alpha = np.sqrt(ke_new / ke)
        velocities *= alpha
        self._drained += ke_new - ke

    @property
    def conserved_offset(self) -> float:
        return -self._drained

    def reset(self) -> None:
        self._drained = 0.0


@dataclass
class Langevin(Thermostat):
    """Ornstein-Uhlenbeck bath; ``friction`` is gamma in 1/ps."""

    temperature: float
    friction: float = 1.0
    name: str = field(default="langevin", init=False)

    _drained: float = field(default=0.0, init=False)

    def half_step(self, velocities, masses, dt, dof, rng) -> None:
        ke_before = kinetic_energy(velocities, masses)
        c1 = np.exp(-self.friction * 0.5 * dt)
        c2 = np.sqrt(1.0 - c1 * c1)
        sigma = np.sqrt(KB * self.temperature / (masses[:, None] * MVV2E))
        velocities *= c1
        velocities += c2 * sigma * rng.standard_normal(velocities.shape)
        self._drained += kinetic_energy(velocities, masses) - ke_before

    @property
    def conserved_offset(self) -> float:
        return -self._drained

    def reset(self) -> None:
        self._drained = 0.0


@dataclass
class NoseHooverChain(Thermostat):
    """Martyna-Tuckerman-Klein chain thermostat."""

    temperature: float
    tau: float = 0.1  #: ps
    length: int = 3
    name: str = field(default="nose-hoover", init=False)

    _xi: np.ndarray = field(default=None, init=False, repr=False)
    _v_xi: np.ndarray = field(default=None, init=False, repr=False)
    _dof: int = field(default=0, init=False)

    def _ensure(self, dof: int) -> None:
        if self._xi is None or self._dof != dof:
            self._xi = np.zeros(self.length)
            self._v_xi = np.zeros(self.length)
            self._dof = dof

    def half_step(self, velocities, masses, dt, dof, rng) -> None:
        self._ensure(dof)
        kt = KB * self.temperature
        q = np.full(self.length, kt * self.tau**2)
        q[0] = dof * kt * self.tau**2

        ke2 = 2.0 * kinetic_energy(velocities, masses)  # = sum m v^2
        dt2, dt4, dt8 = 0.5 * dt, 0.25 * dt, 0.125 * dt

        g = (ke2 - dof * kt) / q[0]
        self._v_xi[-1] += dt4 * (self._v_xi[-2] ** 2 * q[-2] - kt) / q[-1] if self.length > 1 else 0.0
        for m in range(self.length - 2, -1, -1):
            gm = g if m == 0 else (q[m - 1] * self._v_xi[m - 1] ** 2 - kt) / q[m]
            self._v_xi[m] = (
                self._v_xi[m] * np.exp(-dt8 * self._v_xi[m + 1])
                + dt4 * gm
            ) * np.exp(-dt8 * self._v_xi[m + 1])

        scale = float(np.exp(-dt2 * self._v_xi[0]))
        velocities *= scale
        ke2 *= scale * scale
        self._xi += dt2 * self._v_xi

        g = (ke2 - dof * kt) / q[0]
        for m in range(self.length - 1):
            gm = g if m == 0 else (q[m - 1] * self._v_xi[m - 1] ** 2 - kt) / q[m]
            self._v_xi[m] = (
                self._v_xi[m] * np.exp(-dt8 * self._v_xi[m + 1]) + dt4 * gm
            ) * np.exp(-dt8 * self._v_xi[m + 1])
        if self.length > 1:
            self._v_xi[-1] += dt4 * (q[-2] * self._v_xi[-2] ** 2 - kt) / q[-1]

    @property
    def conserved_offset(self) -> float:
        if self._xi is None:
            return 0.0
        kt = KB * self.temperature
        q = np.full(self.length, kt * self.tau**2)
        q[0] = self._dof * kt * self.tau**2
        return float(
            0.5 * np.sum(q * self._v_xi**2)
            + self._dof * kt * self._xi[0]
            + kt * np.sum(self._xi[1:])
        )

    def reset(self) -> None:
        self._xi = None
        self._v_xi = None


def make_thermostat(kind: str, temperature: float, **kw) -> Thermostat:
    kind = (kind or "none").lower()
    table = {
        "none": lambda: NoThermostat(),
        "nve": lambda: NoThermostat(),
        "berendsen": lambda: Berendsen(temperature, **kw),
        "bussi": lambda: Bussi(temperature, **kw),
        "csvr": lambda: Bussi(temperature, **kw),
        "langevin": lambda: Langevin(temperature, **kw),
        "nose-hoover": lambda: NoseHooverChain(temperature, **kw),
        "nhc": lambda: NoseHooverChain(temperature, **kw),
    }
    if kind not in table:
        raise ValueError(f"unknown thermostat {kind!r}; choose from {sorted(table)}")
    return table[kind]()


def maxwell_boltzmann_velocities(masses, temperature, rng, remove_drift=True):
    """Draw 2-D Maxwell-Boltzmann velocities and (optionally) zero the momentum."""
    sigma = np.sqrt(KB * temperature / (np.asarray(masses)[:, None] * MVV2E))
    v = sigma * rng.standard_normal((len(masses), 2))
    if remove_drift:
        v -= np.sum(masses[:, None] * v, axis=0) / np.sum(masses)
    return v


__all__ = [
    "Thermostat",
    "NoThermostat",
    "Berendsen",
    "Bussi",
    "Langevin",
    "NoseHooverChain",
    "make_thermostat",
    "kinetic_energy",
    "temperature_from_ke",
    "maxwell_boltzmann_velocities",
]
