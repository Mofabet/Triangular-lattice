import numpy as np
import pytest

from trilattice import LennardJones
from trilattice.potentials import R_MIN_OVER_SIGMA


@pytest.mark.parametrize("mode", ["cut", "shifted", "shifted-force"])
def test_force_is_minus_the_derivative_of_the_energy(mode):
    p = LennardJones(0.27, 2.88, 8.0, mode=mode)
    r = np.linspace(2.8, 7.9, 400)
    h = 1e-6
    numeric = -(p.energy(r + h) - p.energy(r - h)) / (2 * h)
    assert np.allclose(p.force(r), numeric, rtol=1e-5, atol=1e-7)


def test_minimum_sits_at_the_textbook_position():
    p = LennardJones(1.0, 1.0, 100.0, mode="cut")
    r = R_MIN_OVER_SIGMA
    assert float(p.force(r)) == pytest.approx(0.0, abs=1e-10)
    assert float(p.energy(r)) == pytest.approx(-1.0, abs=1e-12)


def test_repulsive_term_carries_the_factor_of_two():
    # f(r) = 24 eps/r [2 (sig/r)^12 - (sig/r)^6]; dropping the 2 moves the zero
    # of the force from 2^(1/6) sigma to sigma, which is the original bug.
    p = LennardJones(1.0, 1.0, 100.0, mode="cut")
    assert float(p.force(1.0)) > 0.0
    assert float(p.force(R_MIN_OVER_SIGMA)) == pytest.approx(0.0, abs=1e-10)


def test_force_is_attractive_beyond_the_minimum():
    p = LennardJones(0.27, 2.88, 8.0, mode="cut")
    assert np.all(p.force(np.linspace(3.5, 7.5, 50)) < 0.0)


def test_truncation_continuity():
    rc = 8.0
    eps_r = 1e-9
    cut = LennardJones(0.27, 2.88, rc, mode="cut")
    shifted = LennardJones(0.27, 2.88, rc, mode="shifted")
    sf = LennardJones(0.27, 2.88, rc, mode="shifted-force")
    assert abs(float(cut.energy(rc - eps_r))) > 1e-4          # energy jumps
    assert abs(float(shifted.energy(rc - eps_r))) < 1e-9      # energy continuous
    assert abs(float(shifted.force(rc - eps_r))) > 1e-4       # force still jumps
    assert abs(float(sf.energy(rc - eps_r))) < 1e-9
    assert abs(float(sf.force(rc - eps_r))) < 1e-9


def test_lattice_sum_reproduces_the_known_2d_lj_ground_state():
    # For the 2-D triangular LJ crystal the T = 0 spacing is 1.1096 sigma and the
    # cohesive energy is -3.3820 eps (Hoover & co.).  Reproduce to 3 decimals.
    from scipy.optimize import minimize_scalar

    p = LennardJones(1.0, 1.0, 1000.0, mode="cut")
    res = minimize_scalar(lambda a: p.lattice_energy(float(a), n_shells=60),
                          bracket=(1.05, 1.11, 1.20))
    assert float(res.x) == pytest.approx(1.1096, abs=2e-3)
    assert float(res.fun) == pytest.approx(-3.3820, abs=5e-3)


def test_mixing_rules():
    p = LennardJones([1.0, 4.0], [1.0, 2.0], 10.0)
    assert p.eps_table[0, 1] == pytest.approx(2.0)   # geometric
    assert p.sig_table[0, 1] == pytest.approx(1.5)   # arithmetic


def test_tail_corrections_match_numerical_integration():
    from scipy.integrate import quad

    p = LennardJones(0.27, 2.88, 8.0, mode="cut")
    rho, n = 0.1128, 100
    exact = 0.5 * n * rho * quad(lambda r: p.bare_energy(r) * 2 * np.pi * r, 8.0, np.inf)[0]
    assert p.tail_energy(rho, n) == pytest.approx(exact, rel=1e-8)
