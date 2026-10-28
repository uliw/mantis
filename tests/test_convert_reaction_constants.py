"""Unit tests for convert_reaction_constants in mantis.diff_lib."""

import pytest
from mantis.diff_lib import convert_reaction_constants


def test_convert_reaction_constants_basic():
    velde = {
        "POC_fast": "10 1/year",
        "POC_slow": "0.1 1/year",
        "K_Fe3": "1e-3 dimensionless",
    }
    k_values = convert_reaction_constants(velde)
    assert pytest.approx(k_values["POC_fast"], rel=1e-6) == 10.0 / (365.25 * 86400) or pytest.approx(
        k_values["POC_fast"], rel=1e-4
    ) == 3.1688e-7
    assert pytest.approx(k_values["POC_slow"], rel=1e-4) == 3.1688e-9
    assert k_values["K_Fe3"] == 1e-3


def test_convert_reaction_constants_numeric():
    raw = {
        "k1": 1.5,
        "k2": 0,
    }
    k_values = convert_reaction_constants(raw)
    assert k_values["k1"] == 1.5
    assert k_values["k2"] == 0.0


def test_convert_reaction_constants_error_line_number(tmp_path):
    bad_file = tmp_path / "bad_constants.py"
    bad_file.write_text(
        "def get_reaction_constants(pH, phi):\n"
        "    velde = {\n"
        "        'valid_k': '10 1/year',\n"
        "        'bad_k': '10 invalid_unit_xyz',\n"
        "    }\n"
        "    return velde\n"
    )

    import importlib.util

    spec = importlib.util.spec_from_file_location("bad_mod", bad_file)
    bad_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bad_mod)

    velde = bad_mod.get_reaction_constants(7.5, 0.8)

    with pytest.raises(ValueError) as exc_info:
        convert_reaction_constants(velde, filename=bad_file)

    msg = str(exc_info.value)
    assert "bad_k" in msg
    assert "line 4" in msg
    assert str(bad_file) in msg


def test_convert_reaction_constants_non_string_non_number_error(tmp_path):
    bad_file = tmp_path / "bad_constants_type.py"
    bad_file.write_text(
        "def get_reaction_constants(pH, phi):\n"
        "    velde = {\n"
        "        'bad_list': [10, '1/year'],\n"
        "    }\n"
        "    return velde\n"
    )

    import importlib.util

    spec = importlib.util.spec_from_file_location("bad_mod_type", bad_file)
    bad_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bad_mod)

    velde = bad_mod.get_reaction_constants(7.5, 0.8)

    with pytest.raises(ValueError) as exc_info:
        convert_reaction_constants(velde, filename=bad_file)

    msg = str(exc_info.value)
    assert "bad_list" in msg
    assert "line 3" in msg
    assert "must be a string or number" in msg


def test_convert_reaction_constants_with_sed_age():
    import numpy as np
    from mantis.diff_lib import sed_age

    class MockMP:
        z = np.linspace(0, 0.5, 50)
        w = 0.1404 / 100.0  # m/yr
        phi = 0.85
        Db = None

    mp = MockMP()

    # Relative mode
    age = sed_age(mode="relative", a0=1.0, a=0.21, b=1.22)
    constants = {
        "k_POC": "age * 0.01 1/year",
        "K_SO4": "0.9 umol/cm^3",
    }
    k_vals = convert_reaction_constants(constants, mp=mp)
    assert isinstance(k_vals["k_POC"], np.ndarray)
    assert len(k_vals["k_POC"]) == 50
    # Surface age is 0 -> factor is 1.0 -> 0.01 1/year in SI
    assert pytest.approx(k_vals["k_POC"][0], rel=1e-5) == 0.01 / (365.25 * 86400)
    # Decreases with depth
    assert k_vals["k_POC"][-1] < k_vals["k_POC"][0]
    assert pytest.approx(k_vals["K_SO4"], rel=1e-6) == 0.9

    # Direct mode
    age_dir = sed_age(mode="direct", a0=1.0, a=0.21, b=1.22)
    constants_dir = {
        "k_POC": "age",
    }
    k_vals_dir = convert_reaction_constants(constants_dir, mp=mp)
    assert isinstance(k_vals_dir["k_POC"], np.ndarray)
    assert pytest.approx(k_vals_dir["k_POC"][0], rel=1e-5) == 0.21 / (365.25 * 86400)

