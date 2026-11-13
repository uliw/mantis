"""Unit tests for the MANTIS simulation report and metadata generator."""

from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from mantis.diff_lib import data_container
from mantis.model import DiageneticModelResult
from mantis.report_generator import (
    formula_to_elements,
    generate_simulation_report,
    stoichiometry_to_mhchem,
)


def test_stoichiometry_to_mhchem():
    """Verify chemical equation to mhchem formatting."""
    eq1 = "POC_fast + O2 -> CO2"
    tex1 = stoichiometry_to_mhchem(eq1)
    assert r"$\ce{" in tex1
    assert "(CH2O)_f" in tex1
    assert "O2" in tex1
    assert "->" in tex1

    eq2 = "POC_fast + 4 Fe3 -> 4 Fe2_total"
    tex2 = stoichiometry_to_mhchem(eq2)
    assert "4 FeOOH" in tex2
    assert "4 Fe^{2+}_{tot}" in tex2


def test_formula_to_elements():
    """Verify elemental parsing for common geochemical species."""
    assert formula_to_elements("FeS2") == {"Fe": 1.0, "S": 2.0}
    assert formula_to_elements("SO4") == {"S": 1.0, "O": 4.0}
    assert formula_to_elements("HCO3") == {"H": 1.0, "C": 1.0, "O": 3.0}
    assert formula_to_elements("FeOOH") == {"Fe": 1.0, "O": 2.0, "H": 1.0}
    assert formula_to_elements("Fe2", name="Fe2_total") == {"Fe": 1.0}


def test_generate_simulation_report(tmp_path):
    """Verify report generation produces all required CSV tables, Org, and text files."""
    z = np.linspace(0.0, 0.5, 50)
    c_arr = np.full(50, 10.0)

    c = data_container({"SO4": c_arr, "FeS": c_arr})
    mp = data_container({
        "plot_name": "test_sim",
        "phi": 0.8,
        "w": 1e-9,
        "pH": 7.5,
        "species_defs": {
            "SO4": {"phase": "liquid", "formula": "SO4"},
            "FeS": {"phase": "solid", "formula": "FeS"},
        },
        "bc_map": {
            "SO4": {"type": "Dirichlet", "top": 10.0, "bot_type": "Neumann", "bot": 0.0},
            "FeS": {"type": "Robin", "top": 1e-8, "bot_type": "Neumann", "bot": 0.0},
        },
        "reactions": [
            {
                "reaction_name": "sulfate_reduction",
                "reaction": "2 POC_fast + SO4 -> FeS",
                "k_value_name": "k_sr",
                "limiters": ["SO4"],
            }
        ],
        "limiters": {"SO4": "c.SO4 / (c.SO4 + 1.0)"},
        "raw_reaction_constants": {"k_sr": "1.0 1/year"},
    })
    k = data_container({"k_sr": 3.17e-8})
    D_mol = data_container({"SO4": 1e-9, "D_bio": 1e-11})

    result = DiageneticModelResult(
        mp=mp,
        c=c,
        k=k,
        species_list=["SO4", "FeS"],
        z=z,
        D_mol=D_mol,
        diagenetic_reactions=None,
        converged="Yes",
        step=10,
        total_time=1000.0,
        cpu_time=0.5,
        wall_time=0.6,
    )

    out_dir = tmp_path / "logs"
    rep = result.save_report(output_dir=out_dir)

    # Verify all expected files exist
    assert Path(rep["org_path"]).exists()
    assert Path(rep["txt_path"]).exists()

    expected_csvs = [
        "01_simulation_metadata.csv",
        "02_parameters.csv",
        "03_reactions.csv",
        "04_kinetic_rate_laws.csv",
        "05_boundary_conditions_and_fluxes.csv",
        "06_kinetic_constants.csv",
        "07_elemental_mass_balance.csv",
    ]
    for csv_name in expected_csvs:
        csv_file = out_dir / csv_name
        assert csv_file.exists(), f"Missing CSV table: {csv_name}"

    # Verify Org-mode content
    org_text = Path(rep["org_path"]).read_text()
    assert "#+TITLE: MANTIS Diagenetic Simulation Report: test_sim" in org_text
    assert "#+NAME: meta" in org_text
    assert "#+NAME: parameters" in org_text
    assert "#+NAME: reactions" in org_text
    assert "#+NAME: constants" in org_text
    assert "#+NAME: mass_balance" not in org_text


def test_report_generator_with_isotopes(tmp_path):
    z = np.linspace(0, 0.5, 51)
    c = data_container()
    mp = data_container()
    mp.plot_name = "iso_sim"
    mp.z = z
    mp.phi = 0.8
    mp.w = 0.002
    mp.pH = 7.5
    mp.isotopes = True
    mp.reactions = [
        {
            "reaction_name": "sulfate_reduction_fast",
            "reaction": "2 POC_fast + SO4 -> TS2",
            "k_value_name": "POC_fast",
            "isotope_fractionation": {"SO4_S32 to TS2_S32": ["msr_alpha", "SO4_alpha_explicit"]},
        }
    ]
    mp.reaction_parameters = {"msr_alpha": 1.07}
    mp.speciation = {
        "sulfide_acid_base": {
            "type": "acid_base",
            "primary": "TS2",
            "constants": {"pKa1": 7.0},
            "equilibrium_fractionation": {
                "standard": "d34S",
                "alpha": {"H2S_S32_to_HS_S32": 0.99991542},
            },
        }
    }
    mp.species_defs = {
        "SO4": {"phase": "liquid", "diffusion": 1e-9, "bc_top": {"Robin": 28.2}},
        "TS2": {"phase": "liquid", "diffusion": 1e-9, "bc_top": {"Robin": 0.0}},
    }
    for sp in mp.species_defs:
        setattr(c, sp, np.ones_like(z))

    result = DiageneticModelResult(
        mp=mp,
        c=c,
        k=data_container({"k_fast": 1.0}),
        species_list=["SO4", "TS2"],
        z=z,
        D_mol=data_container(),
        diagenetic_reactions=None,
        converged="Yes",
        step=5,
        total_time=100.0,
        cpu_time=0.2,
        wall_time=0.3,
    )

    out_dir = tmp_path / "logs"
    rep = result.save_report(output_dir=out_dir)

    iso_csv = out_dir / "08_isotope_fractionation.csv"
    assert iso_csv.exists()
    df_iso = pd.read_csv(iso_csv)
    assert list(df_iso.columns) == ["Process", "Enrichment factor", "Type"]
    assert "sulfate reduction" in df_iso["Process"].values
    assert "-70 mUr" in df_iso["Enrichment factor"].values
    assert "sulfide speciation" in df_iso["Process"].values
    assert "+2 mUr" in df_iso["Enrichment factor"].values

    txt_text = Path(rep["txt_path"]).read_text()
    assert "ISOTOPE FRACTIONATION & EQUILIBRIUM EFFECTS" in txt_text
    assert "sulfate reduction" in txt_text
    assert "-70 mUr" in txt_text
    assert "sulfide speciation" in txt_text
    assert "+2 mUr" in txt_text

    org_text = Path(rep["org_path"]).read_text()
    assert "#+NAME: isotopes" in org_text
    assert "08_isotope_fractionation.csv" in org_text
