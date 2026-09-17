"""The MD engine: velocity Verlet, optionally thermostatted, with logging.

Velocity Verlet in its kick-drift-kick form::

    v(t + dt/2) = v(t)        + (dt/2) f(t)/m
    r(t + dt)   = r(t)        + dt  v(t + dt/2)
    f(t + dt)   = F(r(t + dt))
    v(t + dt)   = v(t + dt/2) + (dt/2) f(t + dt)/m

Three properties make it, and not the Euler-like scheme of the original, the
right choice: it is second-order accurate, it is symplectic (so the energy error
is *bounded* rather than growing linearly with time), and it is exactly
time-reversible.  The original's ``r += v dt + a dt^2`` is neither second-order
(the factor 1/2 is missing) nor symplectic, so it heats up without limit.

A thermostat, when present, is applied as a half step on each side, giving the
Trotter factorisation described in :mod:`trilattice.thermostats`.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from .forces import ForceField
from .lattice import Box, Configuration
from .neighbors import NeighborList
from .observables import pressure_2d
from .potentials import LennardJones
from .thermostats import (
    NoThermostat,
    Thermostat,
    kinetic_energy,
    maxwell_boltzmann_velocities,
    temperature_from_ke,
)
from .units import FTM2V, KB


@dataclass
class RunLog:
    """Time series collected during a run."""

    step: list[int] = field(default_factory=list)
    time: list[float] = field(default_factory=list)
    temperature: list[float] = field(default_factory=list)
    e_pot: list[float] = field(default_factory=list)
    e_kin: list[float] = field(default_factory=list)
    e_tot: list[float] = field(default_factory=list)
    e_conserved: list[float] = field(default_factory=list)
    pressure: list[float] = field(default_factory=list)
    momentum: list[float] = field(default_factory=list)

    def as_dict(self) -> dict[str, np.ndarray]:
        return {k: np.asarray(v) for k, v in self.__dict__.items()}

    def __len__(self) -> int:
        return len(self.step)


class Simulation:
    """A single MD system: configuration + potential + integrator."""

    def __init__(
        self,
        configuration: Configuration,
        potential: LennardJones,
        *,
        timestep: float = 0.001,  # ps
        thermostat: Thermostat | None = None,
        skin: float = 1.0,
        seed: int | None = None,
        use_numba: bool | None = None,
        remove_drift: bool = True,
    ):
        self.config = configuration.copy()
        self.potential = potential
        self.dt = float(timestep)
        self.thermostat = thermostat or NoThermostat()
        self.rng = np.random.default_rng(seed)
        self.remove_drift = remove_drift

        self.box: Box = self.config.box
        self.positions = self.box.wrap(self.config.positions.copy())
        self.unwrapped = self.positions.copy()
        self.velocities = np.zeros_like(self.positions)
        self.masses = self.config.masses
        self.types = self.config.types

        self.neighbors = NeighborList(self.box, potential.cutoff, skin=skin)
        self.forcefield = ForceField(potential, self.box, use_numba=use_numba)

        self.step_count = 0
        self.elapsed = 0.0  # ps
        self.log = RunLog()
        self._wall_time = 0.0

        self.neighbors.build(self.positions)
        self._evaluate()

    # ------------------------------------------------------------------ #
    @property
    def n_particles(self) -> int:
        return self.positions.shape[0]

    @property
    def dof(self) -> int:
        """2N, minus the two momentum components if the drift is constrained."""
        return 2 * self.n_particles - (2 if self.remove_drift else 0)

    @property
    def kinetic(self) -> float:
        return kinetic_energy(self.velocities, self.masses)

    @property
    def temperature(self) -> float:
        return temperature_from_ke(self.kinetic, self.dof)

    @property
    def total_energy(self) -> float:
        return self.kinetic + self.e_pot

    @property
    def conserved(self) -> float:
        """Quantity that a correct integrator keeps constant, thermostat included."""
        return self.total_energy + self.thermostat.conserved_offset

    @property
    def pressure(self) -> float:
        return pressure_2d(self.temperature, self.virial, self.n_particles, self.box.area)

    # ------------------------------------------------------------------ #
    def set_temperature(self, temperature: float) -> None:
        """Assign fresh Maxwell-Boltzmann velocities."""
        self.velocities = maxwell_boltzmann_velocities(
            self.masses, temperature, self.rng, remove_drift=self.remove_drift
        )

    def zero_momentum(self) -> None:
        p = np.sum(self.masses[:, None] * self.velocities, axis=0)
        self.velocities -= p / np.sum(self.masses)

    # ------------------------------------------------------------------ #
    def _evaluate(self) -> None:
        i, j = self.neighbors.pairs
        res = self.forcefield(self.positions, i, j, self.types)
        self.forces = res.forces
        self.e_pot = res.energy
        self.virial = res.virial
        self.per_atom_energy = res.per_atom_energy

    def _accelerations(self) -> np.ndarray:
        return FTM2V * self.forces / self.masses[:, None]

    def step(self) -> None:
        """One velocity-Verlet step (with thermostat half-steps if any)."""
        dt = self.dt
        self.thermostat.half_step(self.velocities, self.masses, dt, self.dof, self.rng)

        self.velocities += 0.5 * dt * self._accelerations()
        displacement = dt * self.velocities
        self.positions += displacement
        self.unwrapped += displacement
        self.box.wrap_inplace(self.positions)

        self.neighbors.update(self.positions)
        self._evaluate()

        self.velocities += 0.5 * dt * self._accelerations()
        self.thermostat.half_step(self.velocities, self.masses, dt, self.dof, self.rng)

        self.step_count += 1
        self.elapsed += dt

    # ------------------------------------------------------------------ #
    def run(
        self,
        n_steps: int,
        *,
        log_every: int = 100,
        callback=None,
        callback_every: int = 100,
        progress: bool = False,
    ) -> RunLog:
        t0 = time.perf_counter()
        if self.remove_drift:
            self.zero_momentum()
        for k in range(n_steps):
            self.step()
            if log_every and (k + 1) % log_every == 0:
                self._record()
            if callback is not None and (k + 1) % callback_every == 0:
                callback(self)
            if progress and n_steps >= 10 and (k + 1) % max(1, n_steps // 10) == 0:
                print(
                    f"  {100 * (k + 1) // n_steps:3d}%  T = {self.temperature:8.2f} K  "
                    f"E = {self.total_energy / self.n_particles:+.6f} eV/atom",
                    flush=True,
                )
        self._wall_time += time.perf_counter() - t0
        return self.log

    def _record(self) -> None:
        ke = self.kinetic
        t = temperature_from_ke(ke, self.dof)
        self.log.step.append(self.step_count)
        self.log.time.append(self.elapsed)
        self.log.temperature.append(t)
        self.log.e_pot.append(self.e_pot)
        self.log.e_kin.append(ke)
        self.log.e_tot.append(ke + self.e_pot)
        self.log.e_conserved.append(ke + self.e_pot + self.thermostat.conserved_offset)
        self.log.pressure.append(
            pressure_2d(t, self.virial, self.n_particles, self.box.area)
        )
        self.log.momentum.append(
            float(np.linalg.norm(np.sum(self.masses[:, None] * self.velocities, axis=0)))
        )

    def clear_log(self) -> None:
        self.log = RunLog()

    def reset_origin(self) -> None:
        """Restart the unwrapped coordinates -- call after equilibration."""
        self.unwrapped = self.positions.copy()

    # ------------------------------------------------------------------ #
    # ------------------------------------------------------------------ #
    # Changing the particle count mid-run
    # ------------------------------------------------------------------ #
    def remove_particles(self, indices) -> int:
        """Delete particles by index, keeping every array consistent.

        Returns the new particle count.  The neighbour list is rebuilt and the
        forces re-evaluated, so the caller can continue integrating immediately.
        """
        idx = np.atleast_1d(np.asarray(indices, dtype=np.int64))
        if idx.size == 0:
            return self.n_particles
        if idx.size >= self.n_particles:
            raise ValueError("cannot remove every particle")
        keep = np.setdiff1d(np.arange(self.n_particles), idx)
        self.positions = np.ascontiguousarray(self.positions[keep])
        self.velocities = np.ascontiguousarray(self.velocities[keep])
        self.unwrapped = np.ascontiguousarray(self.unwrapped[keep])
        self.masses = np.ascontiguousarray(self.masses[keep])
        self.types = np.ascontiguousarray(self.types[keep])
        self._refresh()
        return self.n_particles

    def insert_particles(self, positions, mass=None, type_id: int = 0,
                         temperature: float | None = None) -> int:
        """Add particles at the given positions, with thermal velocities.

        A newly inserted atom generally overlaps its neighbours -- in a dense
        2-D crystal there is nowhere that does not -- so the caller should follow
        this with a short displacement-capped :meth:`minimize`.
        """
        from .thermostats import maxwell_boltzmann_velocities

        new = np.atleast_2d(np.asarray(positions, dtype=np.float64))
        n_new = new.shape[0]
        m = float(self.masses[0]) if mass is None else float(mass)
        t = self.temperature if temperature is None else float(temperature)
        masses_new = np.full(n_new, m)
        v_new = maxwell_boltzmann_velocities(masses_new, max(t, 1.0), self.rng,
                                             remove_drift=False)
        wrapped_new = self.box.wrap(new)
        self.positions = np.ascontiguousarray(np.vstack([self.positions, wrapped_new]))
        self.velocities = np.ascontiguousarray(np.vstack([self.velocities, v_new]))
        self.unwrapped = np.ascontiguousarray(np.vstack([self.unwrapped, wrapped_new]))
        self.masses = np.ascontiguousarray(np.concatenate([self.masses, masses_new]))
        self.types = np.ascontiguousarray(
            np.concatenate([self.types, np.full(n_new, type_id, dtype=self.types.dtype)]))
        self._refresh()
        return self.n_particles

    def set_state(self, positions, velocities, elapsed=None, unwrapped=None) -> None:
        """Overwrite the dynamical state, e.g. when rewinding to a stored frame."""
        self.positions = np.ascontiguousarray(np.array(positions, dtype=np.float64))
        self.velocities = np.ascontiguousarray(np.array(velocities, dtype=np.float64))
        self.unwrapped = (self.positions.copy() if unwrapped is None
                          else np.ascontiguousarray(np.array(unwrapped, dtype=np.float64)))
        n = self.positions.shape[0]
        if self.masses.shape[0] != n:
            self.masses = np.full(n, float(self.masses[0]))
            self.types = np.zeros(n, dtype=np.int32)
        if elapsed is not None:
            self.elapsed = float(elapsed)
        self._refresh()

    def _refresh(self) -> None:
        """Re-sync the neighbour list and forces after the state changed."""
        self.box.wrap_inplace(self.positions)
        self.neighbors.update(self.positions, force=True)
        self._evaluate()

    def minimize(
        self,
        max_steps: int = 2000,
        f_tol: float = 1e-6,
        dt_max: float | None = None,
        max_move: float = 0.1,
    ):
        """FIRE relaxation to the nearest local minimum of the potential energy.

        Fast Inertial Relaxation Engine (Bitzek et al., PRL 97, 170201).  Used to
        remove the residual stress of a constructed configuration before
        assigning velocities, and to compute T = 0 properties.

        ``max_move`` caps the displacement of the fastest atom in a single
        iteration (angstrom).  Without it a freshly inserted interstitial, which
        starts a/sqrt(3) from three neighbours where u ~ 200 eV, produces a force
        large enough that the first FIRE step throws it across the box and the
        energy overflows.  Capping the step costs nothing on well-behaved
        starting configurations and makes the minimiser usable on bad ones.
        """
        dt = 0.1 * self.dt
        dt_max = 10.0 * self.dt if dt_max is None else dt_max
        alpha, alpha_start = 0.1, 0.1
        f_inc, f_dec, f_alpha, n_min = 1.1, 0.5, 0.99, 5
        steps_since_negative = 0
        self.velocities[:] = 0.0

        for it in range(max_steps):
            self.velocities += dt * self._accelerations()
            power = float(np.sum(self.velocities * self.forces))
            if power > 0.0:
                steps_since_negative += 1
                vnorm = np.linalg.norm(self.velocities)
                fnorm = np.linalg.norm(self.forces)
                if fnorm > 0:
                    self.velocities = (1 - alpha) * self.velocities + alpha * (
                        self.forces / fnorm
                    ) * vnorm
                if steps_since_negative > n_min:
                    dt = min(dt * f_inc, dt_max)
                    alpha *= f_alpha
            else:
                steps_since_negative = 0
                dt *= f_dec
                alpha = alpha_start
                self.velocities[:] = 0.0

            disp = dt * self.velocities
            largest = float(np.max(np.abs(disp))) if disp.size else 0.0
            if largest > max_move:
                disp *= max_move / largest
                self.velocities *= max_move / largest
            self.positions += disp
            self.unwrapped += disp
            self.box.wrap_inplace(self.positions)
            self.neighbors.update(self.positions)
            self._evaluate()

            fmax = float(np.max(np.abs(self.forces))) if self.forces.size else 0.0
            if fmax < f_tol:
                break
        self.velocities[:] = 0.0
        return {"iterations": it + 1, "f_max": fmax, "e_pot": self.e_pot}

    # ------------------------------------------------------------------ #
    def performance(self) -> dict:
        steps = max(1, self.step_count)
        return {
            "wall_time_s": self._wall_time,
            "steps": self.step_count,
            "us_per_step": 1e6 * self._wall_time / steps,
            "ns_per_step_per_atom": 1e9 * self._wall_time / steps / self.n_particles,
            "backend": "numba" if self.forcefield.use_numba else "numpy",
            **self.neighbors.stats(),
        }

    def summary(self) -> str:
        return (
            f"N = {self.n_particles}, {self.box!r}, rho = {self.config.density:.4f} A^-2\n"
            f"{self.potential.summary()}\n"
            f"dt = {self.dt * 1000:.2f} fs, thermostat = {self.thermostat.name}, "
            f"dof = {self.dof}"
        )


__all__ = ["Simulation", "RunLog"]
