"""Trajectory recording and export.

Extended XYZ is written because it is the lingua franca of atomistic
visualisation -- OVITO, VMD and ASE all read it directly, including the
per-frame ``Lattice`` and the extra per-atom columns used here (``psi6`` and the
Voronoi coordination), so a run can be inspected without writing any plotting
code at all.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .lattice import Box


@dataclass
class Trajectory:
    """In-memory frames plus metadata."""

    box: Box
    positions: list[np.ndarray] = field(default_factory=list)
    unwrapped: list[np.ndarray] = field(default_factory=list)
    velocities: list[np.ndarray] = field(default_factory=list)
    times: list[float] = field(default_factory=list)
    types: np.ndarray | None = None

    def append(self, sim) -> None:
        self.positions.append(sim.positions.copy())
        self.unwrapped.append(sim.unwrapped.copy())
        self.velocities.append(sim.velocities.copy())
        self.times.append(sim.elapsed)
        if self.types is None:
            self.types = sim.types.copy()

    def __len__(self) -> int:
        return len(self.positions)

    # -- array views ------------------------------------------------------ #
    @property
    def r(self) -> np.ndarray:
        return np.asarray(self.positions)

    @property
    def u(self) -> np.ndarray:
        return np.asarray(self.unwrapped)

    @property
    def v(self) -> np.ndarray:
        return np.asarray(self.velocities)

    @property
    def t(self) -> np.ndarray:
        return np.asarray(self.times)

    # -- IO ---------------------------------------------------------------- #
    def save_npz(self, path: str | Path) -> None:
        np.savez_compressed(
            path,
            positions=self.r,
            unwrapped=self.u,
            velocities=self.v,
            times=self.t,
            types=self.types,
            box=np.array([self.box.lx, self.box.ly]),
        )

    @classmethod
    def load_npz(cls, path: str | Path) -> "Trajectory":
        d = np.load(path)
        box = Box(float(d["box"][0]), float(d["box"][1]))
        tr = cls(box)
        tr.positions = list(d["positions"])
        tr.unwrapped = list(d["unwrapped"])
        tr.velocities = list(d["velocities"])
        tr.times = list(d["times"])
        tr.types = d["types"]
        return tr

    def write_extxyz(
        self,
        path: str | Path,
        symbols: dict[int, str] | None = None,
        extra: dict[str, list[np.ndarray]] | None = None,
    ) -> None:
        """Write an extended-XYZ trajectory (z = 0 for this 2-D system)."""
        symbols = symbols or {0: "Mg", 1: "Ca"}
        extra = extra or {}
        lattice = f'Lattice="{self.box.lx} 0.0 0.0 0.0 {self.box.ly} 0.0 0.0 0.0 10.0"'
        with open(path, "w", encoding="utf-8") as fh:
            for k, pos in enumerate(self.positions):
                n = pos.shape[0]
                cols = "species:S:1:pos:R:3"
                for name in extra:
                    cols += f":{name}:R:1"
                fh.write(f"{n}\n")
                fh.write(f'{lattice} Properties={cols} Time={self.times[k]:.6f} pbc="T T F"\n')
                vals = [extra[name][k] for name in extra]
                for i in range(n):
                    sym = symbols.get(int(self.types[i]), "X")
                    row = f"{sym} {pos[i, 0]:.6f} {pos[i, 1]:.6f} 0.000000"
                    for v in vals:
                        row += f" {float(v[i]):.6f}"
                    fh.write(row + "\n")


__all__ = ["Trajectory"]
