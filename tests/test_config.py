import numpy as np
import pytest

from trilattice import Settings, read_legacy_start_txt
from trilattice.config import DEFAULT_TOML

LEGACY = """#number_x_entered
5
#number_y_entered
5
#sigma Angstrom
2.88
#sigma_cutoff Angstrom
8
#r
0.4
#a Angstrom
3.2
#m uam
24.305
#bx Angstrom
5
#by Angstrom
5
#epsilon
0.27
#additional_particles
0
#vacancy
0
#termo
300
#iter
1000
"""


def test_legacy_deck_is_parsed_by_key(tmp_path):
    p = tmp_path / "start.txt"
    p.write_text(LEGACY, encoding="utf-8")
    s = read_legacy_start_txt(p)
    assert (s.system.nx, s.system.ny) == (5, 5)
    assert s.system.a == pytest.approx(3.2)
    assert s.system.mass == pytest.approx(24.305)
    assert s.potential.sigma == pytest.approx(2.88)
    assert s.potential.cutoff == pytest.approx(8.0)
    assert s.potential.epsilon == pytest.approx(0.27)
    assert s.run.temperature == pytest.approx(300.0)
    assert s.run.steps == 1000


def test_legacy_parsing_survives_extra_blank_lines(tmp_path):
    """The original decoded by line number and broke on any edit; this must not."""
    p = tmp_path / "start.txt"
    p.write_text("\n\n" + LEGACY.replace("#termo", "\n#termo"), encoding="utf-8")
    assert read_legacy_start_txt(p).run.temperature == pytest.approx(300.0)


def test_toml_roundtrip(tmp_path):
    p = tmp_path / "c.toml"
    p.write_text(DEFAULT_TOML, encoding="utf-8")
    s = Settings.load(p)
    assert s.system.nx == 16
    assert s.potential.mode == "shifted-force"
    assert s.run.thermostat == "bussi"


def test_settings_build_a_consistent_system():
    s = Settings()
    s.system.nx, s.system.ny, s.system.vacancies = 8, 5, 3
    cfg = s.build_configuration()
    assert cfg.n_particles == 2 * 8 * 5 - 3
    pot = s.build_potential()
    assert pot.cutoff < cfg.box.max_cutoff
