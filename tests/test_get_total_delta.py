"""Unit tests for get_total_delta in diff_lib.py."""

import numpy as np
import pytest

from mantis.diff_lib import data_container, get_delta, get_total_delta


def test_get_total_delta_basic():
    # Setup mock mp
    mp = data_container({
        "phi": 0.8,
        "isotope_standards": {
            "S": {"ratio": 0.044162589, "unit": "mUr"},
        },
    })

    # Setup mock concentrations at index -1
    # Dissolved: SO4, TS2 (phi = 0.8)
    # Solid: FeS (1-phi = 0.2)
    c = data_container({
        "SO4": np.array([28.0, 20.0]),
        "SO4_32": np.array([26.8, 19.1]),
        "TS2": np.array([1.0, 2.0]),
        "TS2_32": np.array([0.96, 1.92]),
        "FeS": np.array([5.0, 10.0]),
        "FeS_32": np.array([4.8, 9.6]),
    })

    # Expected at index -1 (second element):
    # phi = 0.8, f_s = 0.2
    # s = 0.8 * (20.0 + 2.0) + 0.2 * 10.0 = 0.8 * 22.0 + 2.0 = 17.6 + 2.0 = 19.6
    # s32 = 0.8 * (19.1 + 1.92) + 0.2 * 9.6 = 0.8 * 21.02 + 1.92 = 16.816 + 1.92 = 18.736
    expected_delta = get_delta(19.6, 18.736, 0.044162589)

    # Test with string species and isotope_standart="S" as requested by user
    res1 = get_total_delta(c, mp, species="SO4, TS2, FeS", isotope_standart="S")
    assert np.isclose(res1, expected_delta)

    # Test with isotope_standard
    res2 = get_total_delta(c, mp, species="SO4, TS2, FeS", isotope_standard="S")
    assert np.isclose(res2, expected_delta)

    # Test with direct float reference ratio
    res3 = get_total_delta(c, mp, species="SO4, TS2, FeS", isotope_standart=0.044162589)
    assert np.isclose(res3, expected_delta)

    # Test default call
    res_default = get_total_delta(c, mp)
    assert np.isclose(res_default, expected_delta)


def test_get_total_delta_stoichiometry_fes2():
    """Verify FeS2 stoichiometry: base gets 2*c while FeS2_32 gets 1*c."""
    mp = data_container({
        "phi": 0.8,
        "isotope_standards": {
            "S": {"ratio": 0.044162589, "unit": "mUr"},
        },
        "species_defs": {
            "SO4": {"phase": "liquid", "isotope": "SO4_32", "stoichiometry": 1.0},
            "FeS2": {"phase": "solid", "isotope": "FeS2_32", "stoichiometry": 2.0},
        },
    })

    c = data_container({
        "SO4": np.array([10.0]),
        "SO4_32": np.array([9.5]),
        "FeS2": np.array([5.0]),
        "FeS2_32": np.array([9.6]),
    })

    # phi = 0.8, f_s = 0.2
    # SO4: liquid -> 0.8 * 1.0 * 10.0 = 8.0, light -> 0.8 * 1.0 * 9.5 = 7.6
    # FeS2: solid -> 0.2 * 2.0 * 5.0 = 2.0, light -> 0.2 * 1.0 * 9.6 = 1.92
    # s = 8.0 + 2.0 = 10.0
    # s32 = 7.6 + 1.92 = 9.52
    expected = get_delta(10.0, 9.52, 0.044162589)

    res = get_total_delta(c, mp, species="SO4, FeS2", isotope_standart="S")
    assert np.isclose(res, expected)


def test_get_total_delta_list_input():
    mp = data_container({
        "phi": 0.5,
        "VCDT": 0.044162589,
    })
    c = data_container({
        "SO4": np.array([10.0]),
        "SO4_32": np.array([9.5]),
    })
    expected = get_delta(0.5 * 10.0, 0.5 * 9.5, 0.044162589)
    res = get_total_delta(c, mp, species=["SO4"], isotope_standard="VCDT")
    assert np.isclose(res, expected)


def test_get_total_delta_profile_none_index():
    mp = data_container({
        "phi": np.array([0.8, 0.7]),
        "VCDT": 0.044162589,
    })
    c = data_container({
        "SO4": np.array([10.0, 20.0]),
        "SO4_32": np.array([9.5, 19.0]),
    })
    res = get_total_delta(c, mp, species="SO4", isotope_standard="VCDT", index=None)
    assert isinstance(res, np.ndarray)
    assert len(res) == 2


def test_get_total_delta_missing_species_error():
    mp = data_container({
        "phi": 0.5,
        "VCDT": 0.044162589,
    })
    c = data_container({
        "SO4": np.array([10.0]),
        "SO4_32": np.array([9.5]),
    })
    with pytest.raises(KeyError, match="Species 'FeS2' not found"):
        get_total_delta(c, mp, species="SO4, FeS2", isotope_standard="VCDT")
