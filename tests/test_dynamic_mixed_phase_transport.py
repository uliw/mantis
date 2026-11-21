"""Tests for dynamic state-dependent mixed-phase transport and isotope retardation."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pytest

from mantis import DiageneticModel
from mantis.diff_lib import data_container, get_flux, is_mixed_species_dynamic
from mantis.monolithic_solver import build_two_phase_mixed_stencil


def _write_model_files(
    target_dir: Path,
    species_code: str,
    reactions_code: str,
    constants_code: str,
) -> None:
    """Write declarative model definition files into a temporary directory."""
    (target_dir / "species.py").write_text(species_code, encoding="utf-8")
    (target_dir / "reactions.py").write_text(reactions_code, encoding="utf-8")
    (target_dir / "kinetic_constants.py").write_text(constants_code, encoding="utf-8")


def test_two_phase_stencil_conservation_and_equilibrium_invariance():
    """Verify exact mass conservation and C_liq = const equilibrium invariance for varying f_dis(z)."""
    N = 50
    L = 0.20
    dx = np.full(N, L / N)
    z = (np.arange(N) + 0.5) * (L / N)
    phi = np.full(N, 0.8)
    f_dis = 0.8 - 0.6 * (z / L)
    D_mol_raw = np.full(N, 1e-9)
    D_bio = np.full(N, 1e-8)

    # 1. Equilibrium invariance when C_liq = f_dis * C_tot = 1.0 (with w=0, D_bio=0)
    bc_dir = {
        "type": "Dirichlet",
        "top": 1.0 / float(f_dis[0]),
        "phase": "mixed",
        "basis": "porewater",
        "bot_type": "Neumann",
        "bot": 0.0,
    }
    ab_eq, b_eq = build_two_phase_mixed_stencil(
        z=z,
        dx=dx,
        phi=phi,
        D_mol_raw=D_mol_raw,
        D_bio=0.0,
        w_pw=0.0,
        w_solid=0.0,
        f_dis=f_dis,
        bc_props=bc_dir,
    )
    c_eq = 1.0 / f_dis
    Tu_eq = ab_eq[1] * c_eq
    Tu_eq[1:] += ab_eq[2, :-1] * c_eq[:-1]
    Tu_eq[:-1] += ab_eq[0, 1:] * c_eq[1:]
    assert np.max(np.abs(Tu_eq - b_eq)) <= 1e-18

    # 2. Exact discrete mass conservation (zero-flux Neumann BCs with D_mol, D_bio > 0)
    bc_nf = {
        "type": "Neumann",
        "top": 0.0,
        "phase": "mixed",
        "basis": "porewater",
        "bot_type": "Neumann",
        "bot": 0.0,
    }
    ab_cons, b_cons = build_two_phase_mixed_stencil(
        z=z,
        dx=dx,
        phi=phi,
        D_mol_raw=D_mol_raw,
        D_bio=D_bio,
        w_pw=0.0,
        w_solid=0.0,
        f_dis=f_dis,
        bc_props=bc_nf,
    )
    assert np.max(np.abs(b_cons)) == 0.0
    col_sums = ab_cons[1].copy()
    col_sums[:-1] += ab_cons[2, :-1]
    col_sums[1:] += ab_cons[0, 1:]
    assert np.max(np.abs(col_sums)) <= 1e-18

    # 3. Verify get_flux returns zero at exact C_liq = const equilibrium
    c_obj = data_container({"Fe2_total": c_eq})
    mp_obj = data_container({
        "phi": phi,
        "w": 0.0,
        "advection": 0.0,
        "Fe2_total_dis_frac": f_dis,
        "species_defs": {"Fe2_total": {"phase": "mixed", "basis": "porewater"}},
    })
    d_mol_obj = data_container({
        "Fe2_total": D_mol_raw,
        "Fe2_total_diss_frac": 1.0,
        "D_bio": 0.0,
    })
    for depth_check in (0.02, 0.10, 0.18):
        flux = get_flux("Fe2_total", depth=depth_check, c=c_obj, mp=mp_obj, z=z, D_mol=d_mol_obj)
        assert abs(flux) <= 1e-18


def _make_sorbent_species_code(k_ads: float) -> str:
    return f"""
isotope_standards = {{
    "d56Fe": {{
        "element": "Fe",
        "standard": "IRMM-014",
        "ratio": 15.68376068,
        "unit": "mUr",
    }},
}}

species = {{
    "Sorbent": {{
        "phase": "solid",
        "include": True,
        "formula": "FeOOH",
        "diffusion": None,
        "bc_top": {{"Dirichlet": 1000.0}},
        "prefill": True,
    }},
    "Fe2_total": {{
        "phase": "mixed",
        "basis": "porewater",
        "include": True,
        "formula": "Fe",
        "diffusion": {{"m0": 3.31, "m1": 0.150}},
        "isotopes": {{"Fe2_total_54": "Fe"}},
        "bc_top": {{"Dirichlet": 1.0, "d56Fe": 0.0}},
        "prefill": False,
    }},
    "CO2": {{"phase": "liquid", "include": False, "formula": "CO2"}},
}}

speciation = {{
    "fe2_partitioning": {{
        "type": "partition",
        "primary": "Fe2_total",
        "constants": {{"K_ads": {k_ads}}},
        "fractions": {{
            "dissolved": "1.0 / (1.0 + K_ads * Sorbent * dstopw)",
            "sorbed": "(K_ads * Sorbent * dstopw) / (1.0 + K_ads * Sorbent * dstopw)",
        }},
        "equilibrium_fractionation": {{
            "standard": "d56Fe",
            "alpha": {{"Fe2_sorbed_54_to_Fe2_liq_54": 0.998}},
        }},
        "species": {{
            "Fe2_liq": {{
                "phase": "liquid",
                "alias_of": "Fe2_total.dissolved",
                "isotopes": {{"Fe2_liq_54": "Fe"}},
            }},
            "Fe2_sorbed": {{
                "phase": "solid",
                "alias_of": "Fe2_total.sorbed",
                "isotopes": {{"Fe2_sorbed_54": "Fe"}},
            }},
        }},
    }},
}}
"""


def test_dynamic_sorbent_mixed_transport_and_chromatographic_isotopes(tmp_path: Path, monkeypatch):
    """Verify dynamic sorbent decay updates f_dis(z, t) and accelerates iron penetration."""
    dyn_dir = tmp_path / "dynamic_sorbent"
    stat_dir = tmp_path / "static_sorbent"
    dyn_dir.mkdir()
    stat_dir.mkdir()

    k_ads = 0.008  # At Sorbent=1000, dstopw=0.25 -> K_ads*Sorbent*dstopw = 2.0 -> f_dis(z=0) = 1/3
    species_code = _make_sorbent_species_code(k_ads)
    reactions_decay = """
reactions = [
    {
        "reaction_name": "sorbent_decay",
        "reaction": "Sorbent -> CO2",
        "k_value_name": "k_decay",
    },
]
"""
    constants_decay = 'kinetic_constants = {"k_decay": "2.0 1/yr"}\n'

    reactions_static = "reactions = []\n"
    constants_static = "kinetic_constants = {}\n"

    _write_model_files(dyn_dir, species_code, reactions_decay, constants_decay)
    _write_model_files(stat_dir, species_code, reactions_static, constants_static)

    monkeypatch.chdir(dyn_dir)
    m_dyn = DiageneticModel(
        base_dir=dyn_dir,
        species="species.py",
        reactions="reactions.py",
        constants="kinetic_constants.py",
        experiment="dyn_sorbent",
    )
    m_dyn.configure(
        grid_setup=[["0 cm", "20 cm", "2 mm", "2 mm"]],
        t_end="1.0 yr",
        dt_init="1 day",
        dt_max="4 day",
        temp=[10.0, 10.0],
        w="1 cm/yr",
        advection="-4 cm/yr",
        phi=0.8,
        isotopes=True,
        delta_threshold_liquid=1e-4,
        delta_threshold_solid=1e-4,
        bioturbation=None,
        bioirrigation=None,
    )
    res_dyn = m_dyn.run()

    # Check live f_dis and f_dis_54 on mp match exact analytical formula on final Sorbent(z)
    sorb_final = np.asarray(res_dyn.c.Sorbent.value)
    assert sorb_final[-1] < 0.3 * sorb_final[0], "Sorbent should decay substantially with depth"
    dstopw = (1.0 - 0.8) / 0.8
    expected_f = 1.0 / (1.0 + k_ads * sorb_final * dstopw)
    expected_f_54 = expected_f / (expected_f + 0.998 * (1.0 - expected_f))
    np.testing.assert_allclose(res_dyn.mp.Fe2_total_dis_frac, expected_f, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(res_dyn.mp.Fe2_total_54_dis_frac, expected_f_54, rtol=1e-10, atol=1e-12)

    # Check chromatographic isotope signature in dynamic run
    df_dyn, _ = res_dyn.save_csv()
    assert abs(df_dyn["d56Fe_Fe2_sorbed"].iloc[0] - df_dyn["d56Fe_Fe2_liq"].iloc[0] - 2.13) < 0.05
    assert df_dyn["d56Fe_Fe2_total"].min() < -1.0

    # Compare penetration depth against static high-sorbent control run
    monkeypatch.chdir(stat_dir)
    m_stat = DiageneticModel(
        base_dir=stat_dir,
        species="species.py",
        reactions="reactions.py",
        constants="kinetic_constants.py",
        experiment="stat_sorbent",
    )
    m_stat.configure(
        grid_setup=[["0 cm", "20 cm", "2 mm", "2 mm"]],
        t_end="1.0 yr",
        dt_init="1 day",
        dt_max="4 day",
        temp=[10.0, 10.0],
        w="1 cm/yr",
        advection="-4 cm/yr",
        phi=0.8,
        isotopes=True,
        delta_threshold_liquid=1e-4,
        delta_threshold_solid=1e-4,
        bioturbation=None,
        bioirrigation=None,
    )
    res_stat = m_stat.run()

    fe_dyn_tot = np.asarray(res_dyn.c.Fe2_total.value)
    fe_stat_tot = np.asarray(res_stat.c.Fe2_total.value)
    z_grid = res_dyn.z
    assert np.sum(fe_dyn_tot * res_dyn.mp.Fe2_total_dis_frac) > np.sum(
        fe_stat_tot * res_stat.mp.Fe2_total_dis_frac
    ), "Dynamic sorbent decay increases f_dis at depth and must yield deeper dissolved Fe2+ penetration"
    # Specifically, in the deep half (> 8 cm) total mass is higher when sorbent decays and advection carries Fe2 faster
    deep_mask = z_grid > 0.08
    assert np.sum(fe_dyn_tot[deep_mask] * res_dyn.mp.Fe2_total_dis_frac[deep_mask]) > np.sum(
        fe_stat_tot[deep_mask] * res_stat.mp.Fe2_total_dis_frac[deep_mask]
    )


def test_dynamic_nonlinear_self_sorption_isotherm(tmp_path: Path, monkeypatch):
    """Verify nonlinear Langmuir self-sorption f_dis(Fe2_total) converges cleanly with isotopes."""
    run_dir = tmp_path / "langmuir"
    run_dir.mkdir()

    k_ads = 3.0
    b_lang = 4.0
    species_code = f"""
isotope_standards = {{
    "d56Fe": {{
        "element": "Fe",
        "standard": "IRMM-014",
        "ratio": 15.68376068,
        "unit": "mUr",
    }},
}}

species = {{
    "Fe2_total": {{
        "phase": "mixed",
        "basis": "porewater",
        "include": True,
        "formula": "Fe",
        "diffusion": {{"m0": 3.31, "m1": 0.150}},
        "isotopes": {{"Fe2_total_54": "Fe"}},
        "bc_top": {{"Dirichlet": 1.0, "d56Fe": 0.0}},
        "prefill": False,
    }},
}}

speciation = {{
    "fe2_partitioning": {{
        "type": "partition",
        "primary": "Fe2_total",
        "constants": {{"K_ads": {k_ads}, "b_lang": {b_lang}}},
        "fractions": {{
            "dissolved": "1.0 / (1.0 + K_ads / (1.0 + b_lang * np.maximum(Fe2_total, 0.0)))",
            "sorbed": "1.0 - (1.0 / (1.0 + K_ads / (1.0 + b_lang * np.maximum(Fe2_total, 0.0))))",
        }},
        "equilibrium_fractionation": {{
            "standard": "d56Fe",
            "alpha": {{"Fe2_sorbed_54_to_Fe2_liq_54": 0.998}},
        }},
        "species": {{
            "Fe2_liq": {{
                "phase": "liquid",
                "alias_of": "Fe2_total.dissolved",
                "isotopes": {{"Fe2_liq_54": "Fe"}},
            }},
            "Fe2_sorbed": {{
                "phase": "solid",
                "alias_of": "Fe2_total.sorbed",
                "isotopes": {{"Fe2_sorbed_54": "Fe"}},
            }},
        }},
    }},
}}
"""
    _write_model_files(run_dir, species_code, "reactions = []\n", "kinetic_constants = {}\n")

    monkeypatch.chdir(run_dir)
    model = DiageneticModel(
        base_dir=run_dir,
        species="species.py",
        reactions="reactions.py",
        constants="kinetic_constants.py",
        experiment="langmuir",
    )
    model.configure(
        grid_setup=[["0 cm", "20 cm", "2 mm", "2 mm"]],
        t_end="0.5 yr",
        dt_init="0.5 day",
        dt_max="2 day",
        temp=[10.0, 10.0],
        w="0 cm/yr",
        advection="-5 cm/yr",
        phi=0.8,
        isotopes=True,
        delta_threshold_liquid=1e-4,
        delta_threshold_solid=1e-4,
        bioturbation=None,
        bioirrigation=None,
    )
    res = model.run()
    df, _ = res.save_csv()

    c_tot = np.asarray(res.c.Fe2_total.value)
    c_54 = np.asarray(res.c.Fe2_total_54.value)
    expected_f = 1.0 / (1.0 + k_ads / (1.0 + b_lang * np.maximum(c_tot, 0.0)))
    expected_f54 = expected_f / (expected_f + 0.998 * (1.0 - expected_f))

    np.testing.assert_allclose(res.mp.Fe2_total_dis_frac, expected_f, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(res.mp.Fe2_total_54_dis_frac, expected_f54, rtol=1e-10, atol=1e-12)

    # Check mass conservation between parent and sub-species
    np.testing.assert_allclose(df["c_Fe2_liq"] + df["c_Fe2_sorbed"], c_tot, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(df["c_Fe2_liq_54"] + df["c_Fe2_sorbed_54"], c_54, rtol=1e-12, atol=1e-12)

    # Non-constant f_dis: higher at top (c=1 -> 1/(1+0.6)=0.625) than at leading edge (c~0 -> 0.25)
    assert res.mp.Fe2_total_dis_frac[0] > 0.55
    assert res.mp.Fe2_total_dis_frac[-1] < 0.30
    assert df["d56Fe_Fe2_total"].min() < -1.0


def test_dynamic_complexation_with_isotopes(tmp_path: Path, monkeypatch):
    """Verify 3-sub-species complexation block with dynamic POC adsorption and isotopes."""
    run_dir = tmp_path / "complexation_iso"
    run_dir.mkdir()

    k_ads = 0.01
    beta_cl = 0.5
    species_code = f"""
isotope_standards = {{
    "d56Fe": {{
        "element": "Fe",
        "standard": "IRMM-014",
        "ratio": 15.68376068,
        "unit": "mUr",
    }},
}}

species = {{
    "POC": {{
        "phase": "solid",
        "include": True,
        "formula": "CH2O",
        "diffusion": None,
        "bc_top": {{"Dirichlet": 1000.0}},
        "prefill": True,
    }},
    "Fe2_total": {{
        "phase": "mixed",
        "basis": "porewater",
        "include": True,
        "formula": "Fe",
        "diffusion": {{"m0": 3.31, "m1": 0.150}},
        "isotopes": {{"Fe2_total_54": "Fe"}},
        "bc_top": {{"Dirichlet": 1.0, "d56Fe": 0.0}},
        "prefill": False,
    }},
    "CO2": {{"phase": "liquid", "include": False, "formula": "CO2"}},
}}

speciation = {{
    "fe2_speciation": {{
        "type": "complexation",
        "primary": "Fe2_total",
        "constants": {{"K_ads": {k_ads}, "beta_Cl": {beta_cl}}},
        "equilibrium_fractionation": {{
            "standard": "d56Fe",
            "alpha": {{"Fe2_sorbed_54_to_Fe2_dis_54": 0.998}},
        }},
        "species": {{
            "Fe2_dis": {{
                "phase": "liquid",
                "fraction": "1.0 / (1.0 + K_ads * POC * dstopw)",
                "isotopes": {{"Fe2_dis_54": "Fe"}},
            }},
            "Fe2_liq": {{
                "phase": "liquid",
                "fraction": "1.0 / ((1.0 + K_ads * POC * dstopw) * (1.0 + beta_Cl))",
            }},
            "Fe2_sorbed": {{
                "phase": "solid",
                "alias_of": "Fe2_total.sorbed",
                "isotopes": {{"Fe2_sorbed_54": "Fe"}},
            }},
        }},
    }},
}}
"""
    reactions_code = """
reactions = [
    {
        "reaction_name": "poc_decay",
        "reaction": "POC -> CO2",
        "k_value_name": "k_poc",
    },
]
"""
    constants_code = 'kinetic_constants = {"k_poc": "2.0 1/yr"}\n'
    _write_model_files(run_dir, species_code, reactions_code, constants_code)

    monkeypatch.chdir(run_dir)
    model = DiageneticModel(
        base_dir=run_dir,
        species="species.py",
        reactions="reactions.py",
        constants="kinetic_constants.py",
        experiment="comp_iso",
    )
    model.configure(
        grid_setup=[["0 cm", "20 cm", "2 mm", "2 mm"]],
        t_end="1.0 yr",
        dt_init="1 day",
        dt_max="4 day",
        temp=[10.0, 10.0],
        w="1 cm/yr",
        advection="-4 cm/yr",
        phi=0.8,
        isotopes=True,
        delta_threshold_liquid=1e-4,
        delta_threshold_solid=1e-4,
        bioturbation=None,
        bioirrigation=None,
    )
    res = model.run()
    df, _ = res.save_csv()

    poc_final = np.asarray(res.c.POC.value)
    dstopw = (1.0 - 0.8) / 0.8
    expected_f = 1.0 / (1.0 + k_ads * poc_final * dstopw)
    expected_f54 = expected_f / (expected_f + 0.998 * (1.0 - expected_f))

    np.testing.assert_allclose(res.mp.Fe2_total_dis_frac, expected_f, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(res.mp.Fe2_total_54_dis_frac, expected_f54, rtol=1e-10, atol=1e-12)

    # Verify all sub-species and isotopes obey exact mass balance
    np.testing.assert_allclose(df["c_Fe2_dis"] + df["c_Fe2_sorbed"], df["c_Fe2_total"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(df["c_Fe2_liq"] * (1.0 + beta_cl), df["c_Fe2_dis"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(df["c_Fe2_dis_54"] + df["c_Fe2_sorbed_54"], df["c_Fe2_total_54"], rtol=1e-12, atol=1e-12)

    # Verify isotope offset & negative chromatographic excursion at leading front
    assert abs(df["d56Fe_Fe2_sorbed"].iloc[0] - df["d56Fe_Fe2_dis"].iloc[0] - 2.13) < 0.05
    assert df["d56Fe_Fe2_total"].min() < -1.0


def test_static_vs_dynamic_mixed_classification():
    """Verify AST classification distinguishes static constants from state-dependent f_dis."""
    base_species = {
        "Fe2_total": {"phase": "mixed"},
        "Fe2_total_54": {"phase": "mixed", "isotope_of": "Fe2_total"},
        "POC": {"phase": "solid"},
    }

    assert not is_mixed_species_dynamic(
        "Fe2_total",
        {**base_species, "Fe2_total": {"phase": "mixed", "partition": {"dissolved": "1.0 / (1.0 + k.Fe2_p_eq)"}}},
        {},
    )
    assert not is_mixed_species_dynamic(
        "Fe2_total_54",
        {**base_species, "Fe2_total": {"phase": "mixed", "partition": {"dissolved": "1.0 / (1.0 + K_d)"}}},
        {},
    )
    assert is_mixed_species_dynamic(
        "Fe2_total",
        {**base_species, "Fe2_total": {"phase": "mixed", "partition": {"dissolved": "1.0 / (1.0 + K_ads * POC * dstopw)"}}},
        {},
    )
    assert is_mixed_species_dynamic(
        "Fe2_total_54",
        {**base_species, "Fe2_total": {"phase": "mixed", "partition": {"dissolved": "1.0 / (1.0 + K_ads * c.POC * dstopw)"}}},
        {},
    )
    carb_species = {
        **base_species,
        "TA": {"phase": "liquid"},
        "DIC": {"phase": "liquid"},
        "Fe2_total": {"phase": "mixed", "partition": {"dissolved": "1.0 / (1.0 + K_ads * 10**(-pH))"}},
    }
    assert is_mixed_species_dynamic("Fe2_total", carb_species, {})
    assert not is_mixed_species_dynamic(
        "Fe2_total",
        {**base_species, "Fe2_total": {"phase": "mixed", "partition": {"dissolved": "1.0 / (1.0 + K_ads * 10**(-pH))"}}},
        {},
        mp=data_container({"is_carbonate": False}),
    )
