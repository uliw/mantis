"""Unit tests for get_flux in diff_lib.py."""

import numpy as np
import pytest

from mantis.diff_lib import data_container, get_flux
from mantis.model import DiageneticModelResult


def test_get_flux_pure_advection():
    """Verify advective flux calculation: J = eff_phi * w * C."""
    z = np.linspace(0.0, 1.0, 100)
    c_arr = np.full(100, 10.0)  # 10 mol/m^3
    phi = 0.8
    w = 1e-9  # m/s

    c = data_container({"SO4": c_arr})
    mp = data_container({
        "phi": phi,
        "w": w,
        "species_defs": {"SO4": {"phase": "liquid"}},
    })

    # For zero gradient, diffusive flux is 0
    # Dissolved flux = phi * w * C = 0.8 * 1e-9 * 10 = 8e-9 mol/(m^2*s)
    flux_bottom = get_flux("SO4", c=c, mp=mp, z=z)
    assert pytest.approx(flux_bottom, rel=1e-5) == 8e-9

    # Check at depth 0.5 m
    flux_mid = get_flux("SO4", depth=0.5, c=c, mp=mp, z=z)
    assert pytest.approx(flux_mid, rel=1e-5) == 8e-9


def test_get_flux_diffusion_and_advection():
    """Verify flux with linear concentration gradient: J = phi*w*C - phi*D*(dC/dz)."""
    z = np.linspace(0.0, 1.0, 100)
    # C(z) = 2 + 3 * z -> dC/dz = 3
    c_arr = 2.0 + 3.0 * z
    phi = 0.5
    w = 2e-8
    D = 1e-8

    c = data_container({"SO4": c_arr})
    mp = data_container({
        "phi": phi,
        "w": w,
        "species_defs": {"SO4": {"phase": "liquid"}},
    })
    D_mol = data_container({"SO4": D, "D_bio": 0.0})

    # At z = 1.0 (bottom):
    # C = 5.0, dC/dz = 3.0
    # J = 0.5 * 2e-8 * 5.0 - 0.5 * 1e-8 * 3.0 = 5e-8 - 1.5e-8 = 3.5e-8
    flux_bot = get_flux("SO4", c=c, mp=mp, z=z, D_mol=D_mol)
    assert pytest.approx(flux_bot, rel=1e-3) == 3.5e-8

    # At z = 0.5:
    # C = 3.5, dC/dz = 3.0
    # J = 0.5 * 2e-8 * 3.5 - 0.5 * 1e-8 * 3.0 = 3.5e-8 - 1.5e-8 = 2.0e-8
    flux_mid = get_flux("SO4", depth="50 cm", c=c, mp=mp, z=z, D_mol=D_mol)
    assert pytest.approx(flux_mid, rel=1e-3) == 2.0e-8


def test_get_flux_solid_species():
    """Verify solid species flux scaling by (1 - phi)."""
    z = np.linspace(0.0, 1.0, 100)
    c_arr = np.full(100, 20.0)  # solid conc
    phi = 0.7  # 1 - phi = 0.3
    w = 1e-9

    c = data_container({"FeS": c_arr})
    mp = data_container({
        "phi": phi,
        "w": w,
        "species_defs": {"FeS": {"phase": "solid"}},
    })

    # J = (1 - phi) * w * C = 0.3 * 1e-9 * 20 = 6e-9 mol/(m^2*s)
    flux_bot = get_flux("FeS", c=c, mp=mp, z=z)
    assert pytest.approx(flux_bot, rel=1e-5) == 6e-9


def test_get_flux_via_result_and_unit_conversion():
    """Verify passing DiageneticModelResult and using unit conversion."""
    z = np.linspace(0.0, 1.0, 100)
    c_arr = np.full(100, 10.0)
    phi = 0.8
    w = 1e-9  # m/s

    c = data_container({"SO4": c_arr})
    mp = data_container({
        "phi": phi,
        "w": w,
        "species_defs": {"SO4": {"phase": "liquid"}},
    })
    D_mol = data_container({"SO4": 0.0, "D_bio": 0.0})

    result = DiageneticModelResult(
        mp=mp,
        c=c,
        k=data_container(),
        species_list=["SO4"],
        z=z,
        D_mol=D_mol,
        diagenetic_reactions=None,
        converged="True",
        step=1,
        total_time=100.0,
    )

    # Internal flux = 8e-9 mol/(m^2*s)
    flux_default = get_flux("SO4", result=result)
    assert pytest.approx(flux_default, rel=1e-5) == 8e-9

    # ModelResult method
    flux_method = result.get_flux("SO4")
    assert pytest.approx(flux_method, rel=1e-5) == 8e-9

    # Unit conversion to mmol/(m^2*s)
    flux_mmol = result.get_flux("SO4", units="mmol/(m^2*s)")
    assert pytest.approx(flux_mmol, rel=1e-5) == 8e-6
