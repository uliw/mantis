"""test_equilibrium_speciation_isotopes.py

Tests declarative equilibrium isotope fractionation in speciation systems (e.g. HS- <-> H2S).
"""

import numpy as np
import pytest

from mantis.diff_lib import (
    data_container,
    get_delta,
    partition_equilibrium_isotope,
    update_all_species,
)
from mantis.species_parser import validate_and_normalize_speciation


def test_partition_equilibrium_isotope_mass_and_ratio():
    """Verify analytical 2-species isotope partitioning and isotopic mass conservation."""
    std_ratio = 0.044162589
    c_tot = 2.0  # mmol/L total sulfide
    # Suppose bulk TS2 has delta = +10 permil
    r_sample = std_ratio * (1.0 + 10.0 / 1000.0)
    c_32 = c_tot / (1.0 + r_sample)

    f_hs = 0.8
    f_h2s = 0.2
    alpha_eq = 0.998  # (34S/32S)_h2s / (34S/32S)_hs ~ 1.002

    hs_32 = partition_equilibrium_isotope(c_32, f_hs, f_h2s, alpha_eq)
    h2s_32 = c_32 - hs_32

    # 1. Mass conservation
    assert np.isclose(hs_32 + h2s_32, c_32, rtol=1e-14)

    # 2. Delta calculation
    c_hs = c_tot * f_hs
    c_h2s = c_tot * f_h2s
    d_hs = get_delta(c_hs, hs_32, std_ratio)
    d_h2s = get_delta(c_h2s, h2s_32, std_ratio)

    # R_sample definitions
    r_hs = (c_hs - hs_32) / hs_32
    r_h2s = (c_h2s - h2s_32) / h2s_32

    # Under partition_equilibrium_isotope(c_iso, f1, f2, alpha):
    # denom = f1 + alpha * f2
    # c1_iso = c_iso * f1 / (f1 + alpha * f2)
    # c2_iso = c_iso * alpha * f2 / (f1 + alpha * f2)
    # Since R = (c - c_iso) / c_iso = c / c_iso - 1, we have:
    # 1 + R1 = (c / c_iso) * (f1 + alpha * f2)
    # 1 + R2 = (c / c_iso) * (f1 + alpha * f2) / alpha
    # Therefore: (1 + R2) / (1 + R1) = 1.0 / alpha
    assert np.isclose((1.0 + r_h2s) / (1.0 + r_hs), 1.0 / alpha_eq, rtol=1e-12)
    assert d_h2s > d_hs  # H2S enriched in 34S relative to HS-


def test_update_all_species_speciation_isotope_fractionation():
    """Verify that update_all_species evaluates child isotope companions from speciation block."""
    raw_speciation = {
        "sulfide_acid_base": {
            "type": "acid_base",
            "primary": "TS2",
            "constants": {"pKa1": 7.0},
            "equilibrium_fractionation": {
                "standard": "d34S",
                "alpha": {"H2S_S32_to_HS_S32": 0.998},
            },
            "species": {
                "HS": {
                    "phase": "liquid",
                    "fraction": "Ka1 / (Hplus + Ka1)",
                    "isotopes": {"HS_S32": "S"},
                },
                "H2S": {
                    "phase": "liquid",
                    "fraction": "Hplus / (Hplus + Ka1)",
                    "isotopes": {"H2S_S32": "S"},
                },
            },
        }
    }

    species_defs = {
        "TS2": {
            "phase": "liquid",
            "include": True,
            "isotopes": {"TS2_S32": "S"},
        },
        "TS2_S32": {
            "phase": "liquid",
            "include": True,
            "isotope_of": "TS2",
        },
    }

    norm_spec = validate_and_normalize_speciation(raw_speciation, species_defs)

    mp = data_container({
        "phi": 0.8,
        "pH": 7.0,  # pH = pKa => 50% HS-, 50% H2S
        "species_defs": species_defs,
        "speciation": norm_spec,
        "isotopes": True,
        "isotope_standards": {
            "d34S": {"element": "S", "standard": "VCDT", "ratio": 0.044162589}
        },
    })

    n_cells = 10
    c_ts2 = np.full(n_cells, 1.5)
    std_r = 0.044162589
    c_ts2_32 = c_ts2 / (1.0 + std_r)

    class CellVar:
        def __init__(self, arr):
            self.value = np.array(arr, dtype=np.float64)

        def setValue(self, val):
            self.value[:] = val

        def __setitem__(self, s, val):
            self.value[s] = val

    c = data_container({
        "TS2": CellVar(c_ts2),
        "TS2_S32": CellVar(c_ts2_32),
        "HS": CellVar(np.zeros(n_cells)),
        "H2S": CellVar(np.zeros(n_cells)),
        "HS_S32": CellVar(np.zeros(n_cells)),
        "H2S_S32": CellVar(np.zeros(n_cells)),
    })

    update_all_species(c, mp)

    # Verify bulk concentrations
    assert np.allclose(c.HS.value, 0.75)
    assert np.allclose(c.H2S.value, 0.75)

    # Verify isotope companions are populated and sum to parent
    assert np.all(c.HS_S32.value > 0.0)
    assert np.all(c.H2S_S32.value > 0.0)
    assert np.allclose(c.HS_S32.value + c.H2S_S32.value, c.TS2_S32.value, rtol=1e-14)

    # Verify fractionation direction
    d_hs = get_delta(c.HS.value, c.HS_S32.value, std_r)
    d_h2s = get_delta(c.H2S.value, c.H2S_S32.value, std_r)
    assert np.all(d_h2s > d_hs)


def test_equilibrium_fractionation_invalid_species_raises_error():
    """Verify that specifying an unknown species key in equilibrium_fractionation exits with error."""
    raw_speciation = {
        "sulfide_acid_base": {
            "type": "acid_base",
            "primary": "TS2",
            "constants": {"pKa1": 7.0},
            "equilibrium_fractionation": {
                "standard": "d34S",
                "alpha": {"Unknown_to_HS_S32": 0.998},
            },
            "species": {
                "HS": {
                    "phase": "liquid",
                    "fraction": "Ka1 / (Hplus + Ka1)",
                    "isotopes": {"HS_S32": "S"},
                },
                "H2S": {
                    "phase": "liquid",
                    "fraction": "Hplus / (Hplus + Ka1)",
                    "isotopes": {"H2S_S32": "S"},
                },
            },
        }
    }
    species_defs = {
        "TS2": {"phase": "liquid", "include": True, "isotopes": {"TS2_S32": "S"}},
        "TS2_S32": {"phase": "liquid", "include": True, "isotope_of": "TS2"},
    }
    with pytest.raises(SystemExit):
        validate_and_normalize_speciation(raw_speciation, species_defs)

