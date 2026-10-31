"""test_carbonate_chemistry.py

Unit and integration tests for the marine carbonate equilibrium chemistry module.
"""

from pathlib import Path

import pytest
import numpy as np
from mantis.seawater_carbonate import (
    calc_seawater_density,
    calc_carbonate_constants,
    solve_carbonate_speciation_from_ta_dic,
    solve_carbonate_speciation_from_dic_ph,
    solve_carbonate_speciation_from_ta_ph,
    adjust_for_calcite_equilibrium,
)


def test_seawater_density_unesco():
    """Test UNESCO 1981 equation of state against standard benchmark."""
    # At 10 C, 35 PSU, density is approx 1026.9 kg/m^3
    rho = calc_seawater_density(temp_c=10.0, salinity=35.0)
    assert 1026.0 < rho < 1028.0


def test_carbonate_constants_benchmarks():
    """Test K1, K2, and Ksp calculations against oceanographic benchmarks on free and total scales."""
    # 1. Default Free pH scale
    const_free = calc_carbonate_constants(temp_c=10.0, salinity=34.0, pressure_bar=1.0)
    assert const_free.ph_scale == "free"
    assert 6.00 < const_free.pK1 < 6.12
    assert 9.20 < const_free.pK2 < 9.35
    assert 0.40 < const_free.Ksp_calcite < 0.48

    # 2. Total pH scale (Lueker et al. 2000: pK1 ~ 6.00, pK2 ~ 9.22)
    const_tot = calc_carbonate_constants(temp_c=10.0, salinity=34.0, pressure_bar=1.0, ph_scale="total")
    assert const_tot.ph_scale == "total"
    assert 5.95 < const_tot.pK1 < 6.05
    assert 9.15 < const_tot.pK2 < 9.30
    assert 0.40 < const_tot.Ksp_calcite < 0.48

    # Scale factor check: pK_free - pK_total = log10(scale_factor)
    assert const_free.pK1 == pytest.approx(const_tot.pK1 + np.log10(const_free.scale_factor), abs=1e-8)
    assert const_free.pK2 == pytest.approx(const_tot.pK2 + np.log10(const_free.scale_factor), abs=1e-8)


def test_adjust_for_calcite_equilibrium_combinations():
    """Test all 3 boundary combinations with adjust_omega = True."""
    const = calc_carbonate_constants(temp_c=10.0, salinity=34.0, pressure_bar=1.0)
    ca2_top = 10.3  # mol/m^3

    # 1. pH & DIC -> adjust TA
    spec1, msg1 = adjust_for_calcite_equilibrium(
        ca2_top=ca2_top, dic=2.2, ph=7.8, const=const
    )
    omega1 = (ca2_top * spec1["CO3"]) / const.Ksp_calcite
    assert pytest.approx(omega1, rel=1e-5) == 1.0
    assert "Adjusted TA" in msg1

    # 2. pH & TA -> adjust DIC
    spec2, msg2 = adjust_for_calcite_equilibrium(
        ca2_top=ca2_top, ta=2.3, ph=7.8, const=const
    )
    omega2 = (ca2_top * spec2["CO3"]) / const.Ksp_calcite
    assert pytest.approx(omega2, rel=1e-5) == 1.0
    assert "Adjusted DIC" in msg2

    # 3. DIC & TA -> adjust pH
    spec3, msg3 = adjust_for_calcite_equilibrium(
        ca2_top=ca2_top, dic=2.2, ta=2.4, const=const
    )
    omega3 = (ca2_top * spec3["CO3"]) / const.Ksp_calcite
    assert pytest.approx(omega3, rel=1e-5) == 1.0
    assert "Adjusted pH" in msg3


def test_missing_ca2_raises_error():
    """Test that missing or non-positive Ca2 raises a ValueError."""
    const = calc_carbonate_constants(temp_c=10.0, salinity=34.0, pressure_bar=1.0)
    with pytest.raises(ValueError, match="requires species 'Ca2'"):
        adjust_for_calcite_equilibrium(ca2_top=None, dic=2.2, ta=2.4, const=const)

    with pytest.raises(ValueError, match="requires species 'Ca2'"):
        adjust_for_calcite_equilibrium(ca2_top=-1.0, dic=2.2, ta=2.4, const=const)


FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "carbonate_test"


def test_model_carbonate_overrides_species_bc():
    """Test that carbonate_equilibrium_chemistry overrides any TA, DIC, pH in species.py."""
    from mantis import DiageneticModel

    eq_file = Path(__file__).resolve().parent / "test_carbonate_chemistry_equations.py"
    model = DiageneticModel(
        base_dir=FIXTURES_DIR,
        species="species.py",
        reactions="fes_reactions.py",
        constants="kinetic_constants_fast.py",
        equations=eq_file,
        layout="plot_layout.py",
    )
    model.configure(
        w="0.2 cm/yr",
        phi=0.8,
        carbonate_equilibrium_chemistry={
            "pressure": 1,
            "salinity": 34,
            "TA": 2.4,
            "DIC": 2.2,
            "adjust_omega": True,
        },
    )
    model.setup()
    assert model.mp.is_carbonate is True
    assert hasattr(model.mp, "carb_spec_top")
    assert model.mp.carb_spec_top["TA"] == pytest.approx(2.182, rel=1e-2)


def test_compute_carbonate_speciation_with_include_in_ta():
    """Test that declarative include_in_TA properly subtracts non-carbonate alkalinity."""
    from mantis.diff_lib import compute_carbonate_speciation

    class Container:
        pass

    class Params:
        pass

    # 1. Base case: TA = 2.4, DIC = 2.2, no non-carbonate species
    c1 = Container()
    c1.TA = np.array([2.4])
    c1.DIC = np.array([2.2])
    mp1 = Params()
    mp1.is_carbonate = True
    mp1.K1_carb = 10.0 ** (-6.0)
    mp1.K2_carb = 10.0 ** (-9.0)

    compute_carbonate_speciation(c1, mp1)
    ph_base = float(c1.pH[0])

    # 2. Case with TS2 and declarative include_in_TA: 1.0 on HS
    c2 = Container()
    c2.TA = np.array([2.4])
    c2.DIC = np.array([2.2])
    c2.TS2 = np.array([1.0])  # 1.0 mM sulfide
    mp2 = Params()
    mp2.is_carbonate = True
    mp2.K1_carb = 10.0 ** (-6.0)
    mp2.K2_carb = 10.0 ** (-9.0)
    mp2.species_defs = {
        "TS2": {"phase": "liquid", "include": True},
        "HS": {
            "phase": "liquid",
            "include": False,
            "equilibrium_of": "TS2",
            "fraction": "Ka1 / (Hplus + Ka1)",
            "constants": {"pKa1": 7.0},
            "include_in_TA": 1.0,
        },
    }

    compute_carbonate_speciation(c2, mp2)
    ph_with_hs = float(c2.pH[0])
    ta_carb = float(c2.ta_carb[0])
    ta_non_carb = float(c2.ta_non_carb[0])

    # HS- takes up alkalinity, so carbonate alkalinity is reduced: CA < TA
    assert ta_carb < 2.4
    assert ta_non_carb > 0.0
    assert pytest.approx(ta_carb + ta_non_carb, rel=1e-6) == 2.4
    # With lower CA, pH must be more acidic than without sulfide alkalinity
    assert ph_with_hs < ph_base


def test_compute_carbonate_speciation_negative_factor():
    """Test that a negative include_in_TA factor (proton donor) increases carbonate alkalinity."""
    from mantis.diff_lib import compute_carbonate_speciation

    class Container:
        pass

    class Params:
        pass

    c = Container()
    c.TA = np.array([2.4])
    c.DIC = np.array([2.2])
    c.Acid = np.array([0.2])
    mp = Params()
    mp.is_carbonate = True
    mp.K1_carb = 10.0 ** (-6.0)
    mp.K2_carb = 10.0 ** (-9.0)
    mp.species_defs = {
        "Acid": {
            "phase": "liquid",
            "include": True,
            "include_in_TA": -1.0,
        },
    }

    compute_carbonate_speciation(c, mp)
    ta_carb = float(c.ta_carb[0])
    ta_non_carb = float(c.ta_non_carb[0])

    assert pytest.approx(ta_non_carb, rel=1e-6) == -0.2
    assert pytest.approx(ta_carb, rel=1e-6) == 2.6  # 2.4 - (-0.2) = 2.6
    assert ta_carb > 2.4


def test_species_parser_rejects_non_numeric_include_in_ta():
    """Test that species parser rejects invalid non-numeric include_in_TA."""
    from mantis.species_parser import parse_and_validate_species_dict

    species = {
        "HS": {
            "phase": "liquid",
            "include": False,
            "include_in_TA": "invalid_string",
        }
    }
    with pytest.raises(SystemExit):
        parse_and_validate_species_dict(species, {})


def test_carbonate_configured_via_speciation_block():
    """Test that defining carbonate system under speciation in species.py automatically configures carbonate chemistry."""
    from mantis import DiageneticModel
    from pathlib import Path

    simplefe_dir = Path(__file__).resolve().parent.parent / "nbk" / "SimpleFe"
    if not (simplefe_dir / "species.py").exists():
        pytest.skip("SimpleFe directory not present.")

    model = DiageneticModel(
        base_dir=simplefe_dir,
        species="species.py",
        reactions="reactions.py",
        constants="kinetic_constants.py",
        equations="equations.py",
    )
    # Configure without passing carbonate_equilibrium_chemistry or iron speciation
    model.configure(
        w="0.14 cm/yr",
        phi=0.8,
    )
    model.setup()

    assert model.mp.is_carbonate is True
    assert model.mp.enable_iron_speciation is True
    assert hasattr(model.mp, "carb_spec_top")
    assert model.mp.carb_constants.pressure_bar == 50.0
    assert model.mp.carb_constants.salinity == 35.0
    assert model.mp.carb_spec_top["DIC"] == pytest.approx(2.345, rel=1e-3)
    assert model.mp.carb_spec_top["pH"] == pytest.approx(7.59, rel=1e-3)
    assert hasattr(model.mp, "beta_FeHS")
    assert model.mp.beta_FeHS == pytest.approx(251.1886, rel=1e-3)


