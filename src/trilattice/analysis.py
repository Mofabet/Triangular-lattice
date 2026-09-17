"""Statistics for correlated time series.

An MD time series is correlated, so the naive standard error ``s/sqrt(n)``
underestimates the true uncertainty by ``sqrt(2 tau_int)`` -- often a factor of
five or more.  Reporting error bars without this correction is the most common
way to publish a number that disagrees with everyone else's at "8 sigma".
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Estimate:
    mean: float
    error: float
    tau_int: float
    n_effective: float

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return f"{self.mean:.6g} +/- {self.error:.2g}"


def autocorrelation(x: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Normalised autocorrelation function, computed via FFT."""
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    n = x.size
    if n < 2:
        return np.ones(1)
    size = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(x, size)
    acf = np.fft.irfft(f * np.conj(f), size)[:n].real
    acf /= acf[0] if acf[0] != 0 else 1.0
    return acf if max_lag is None else acf[: max_lag + 1]


def integrated_autocorrelation_time(x: np.ndarray, c: float = 6.0) -> float:
    """``tau_int`` with Sokal's automatic windowing."""
    acf = autocorrelation(x)
    taus = 2.0 * np.cumsum(acf) - 1.0
    window = np.arange(acf.size)
    ok = window >= c * taus
    idx = int(np.argmax(ok)) if ok.any() else acf.size - 1
    return float(max(taus[idx], 0.5))


def estimate(x: np.ndarray) -> Estimate:
    """Mean with a correlation-corrected standard error."""
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    if n < 2:
        return Estimate(float(x.mean()) if n else float("nan"), float("nan"), float("nan"), n)
    tau = integrated_autocorrelation_time(x)
    n_eff = n / (2.0 * tau)
    return Estimate(float(x.mean()), float(np.std(x, ddof=1) / np.sqrt(n_eff)), tau, n_eff)


def block_average(x: np.ndarray, n_blocks: int = 8) -> Estimate:
    """Blocking estimator -- a robust cross-check on :func:`estimate`."""
    x = np.asarray(x, dtype=np.float64)
    n_blocks = max(2, min(n_blocks, x.size // 2))
    blocks = np.array_split(x[: (x.size // n_blocks) * n_blocks], n_blocks)
    means = np.array([b.mean() for b in blocks])
    return Estimate(
        float(means.mean()),
        float(means.std(ddof=1) / np.sqrt(n_blocks)),
        float("nan"),
        n_blocks,
    )


def jackknife(x: np.ndarray, func, n_blocks: int = 16) -> Estimate:
    """Jackknife error for a nonlinear estimator (heat capacity, for example)."""
    x = np.asarray(x, dtype=np.float64)
    n_blocks = max(2, min(n_blocks, x.size // 2))
    blocks = np.array_split(x[: (x.size // n_blocks) * n_blocks], n_blocks)
    full = func(np.concatenate(blocks))
    partial = np.array(
        [func(np.concatenate([b for m, b in enumerate(blocks) if m != k])) for k in range(n_blocks)]
    )
    err = np.sqrt((n_blocks - 1) / n_blocks * np.sum((partial - partial.mean()) ** 2))
    return Estimate(float(full), float(err), float("nan"), n_blocks)


def find_transition(x: np.ndarray, y: np.ndarray) -> float:
    """Location of the steepest descent of ``y(x)`` -- a crude T_m estimator."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 3:
        return float("nan")
    d = np.gradient(y, x)
    return float(x[int(np.argmin(d))])


__all__ = [
    "Estimate",
    "autocorrelation",
    "integrated_autocorrelation_time",
    "estimate",
    "block_average",
    "jackknife",
    "find_transition",
]
