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
