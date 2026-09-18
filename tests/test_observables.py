import numpy as np
import pytest

from trilattice import Box, triangular_lattice
from trilattice import observables as obs


def test_psi6_is_one_for_the_perfect_lattice():
    cfg = triangular_lattice(8, 5, 3.2)
    assert obs.global_psi6(cfg.positions, cfg.box, 1.35 * 3.2) == pytest.approx(1.0, abs=1e-12)


def test_psi6_is_small_for_a_random_gas():
    rng = np.random.default_rng(0)
    box = Box(60.0, 60.0)
    p = rng.uniform(0, 60, size=(400, 2))
    assert obs.global_psi6(p, box, 1.35 * 3.2) < 0.25


def test_every_site_is_six_coordinated_in_the_perfect_lattice():
    cfg = triangular_lattice(8, 5, 3.2)
    assert np.all(obs.coordination_by_delaunay(cfg.positions, cfg.box) == 6)
    assert obs.defect_fraction(cfg.positions, cfg.box) == 0.0


def test_a_vacancy_creates_defects():
    from trilattice import add_defects

    cfg = triangular_lattice(10, 6, 3.2)
    holed = add_defects(cfg, n_vacancies=1, rng=np.random.default_rng(0))
    assert obs.defect_fraction(holed.positions, holed.box) > 0.0


def test_rdf_of_an_ideal_gas_is_one():
    rng = np.random.default_rng(1)
    box = Box(120.0, 120.0)
    acc = obs.RDFAccumulator(box, r_max=20.0, n_bins=40)
    for _ in range(40):
        acc.accumulate(rng.uniform(0, 120, size=(1500, 2)))
    r, g = acc.result()
    m = r > 4.0
    assert np.allclose(g[m], 1.0, atol=0.06)


def test_rdf_peaks_sit_on_the_lattice_shells():
    from trilattice import perfect_lattice_neighbour_shells

    a = 3.2
    cfg = triangular_lattice(14, 8, a)
    acc = obs.RDFAccumulator(cfg.box, r_max=11.0, n_bins=550)
    acc.accumulate(cfg.positions)
    r, g = acc.result()
    peaks = r[g > 0]
    radii, _ = perfect_lattice_neighbour_shells(a, 11.0)
    for shell in radii[:4]:
        assert np.min(np.abs(peaks - shell)) < 0.05


def test_coordination_number_of_the_first_shell_is_six():
    a = 3.2
    cfg = triangular_lattice(14, 8, a)
    acc = obs.RDFAccumulator(cfg.box, r_max=11.0, n_bins=1100)
    acc.accumulate(cfg.positions)
    r, g = acc.result()
    z = obs.coordination_number(g, r, cfg.density, 1.3 * a)
    assert z == pytest.approx(6.0, rel=0.02)


def test_structure_factor_has_bragg_peaks_of_height_n():
    cfg = triangular_lattice(8, 5, 3.2)
    _, _, s = obs.structure_factor(cfg.positions, cfg.box, n_max=12)
    assert s.max() == pytest.approx(cfg.n_particles, rel=1e-8)


def test_diffusion_of_ballistic_motion():
    # r(t) = v t  ->  MSD = v^2 t^2 ; the linear fit is a sanity check of units
    t = np.linspace(0, 10, 200)
    msd = 4.0 * 0.5 * t          # D = 0.5 by construction
    assert obs.diffusion_coefficient(t, msd) == pytest.approx(0.5, rel=1e-6)


def test_pressure_of_an_ideal_gas():
    from trilattice.units import KB

    p = obs.pressure_2d(temperature=300.0, virial=0.0, n_particles=100, area=1000.0)
    assert p == pytest.approx(100 * KB * 300.0 / 1000.0)


def test_lindemann_is_zero_for_a_frozen_crystal():
    cfg = triangular_lattice(8, 5, 3.2)
    traj = np.repeat(cfg.positions[None], 10, axis=0)
    assert obs.lindemann_2d(traj, cfg.box, 3.2) == pytest.approx(0.0)


def test_heat_capacity_of_a_harmonic_2d_solid_is_two_kb():
    """Dulong-Petit in 2-D.  Catches the missing k_B^2 in the normalisation."""
    from trilattice.units import KB

    rng = np.random.default_rng(0)
    n, T = 480, 300.0
    # energies distributed with the variance a 2 k_B/atom solid must have
    var_target = 2.0 * n * KB**2 * T**2
    e = rng.normal(0.0, np.sqrt(var_target), size=400000)
    assert obs.heat_capacity_nvt(e, T, n) == pytest.approx(2.0, rel=0.02)


def test_structure_factor_k_max_gives_a_square_window_on_any_box():
    from trilattice import Box

    rng = np.random.default_rng(0)
    box = Box(40.0, 90.0)                      # deliberately anisotropic
    pos = rng.uniform(0, 1, size=(200, 2)) * box.lengths
    kx, ky, s = obs.structure_factor(pos, box, k_max=3.0)
    assert kx.max() == pytest.approx(3.0, abs=2 * np.pi / box.lx)
    assert ky.max() == pytest.approx(3.0, abs=2 * np.pi / box.ly)
    assert s.shape == (kx.size, ky.size)
    # n_max instead makes the window depend on the box: that is the old bug
    kx2, ky2, _ = obs.structure_factor(pos, box, n_max=16)
    assert kx2.max() > 2 * ky2.max()
