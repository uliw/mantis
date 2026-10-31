"""species.py

Declarative Chemical Species Definitions for the SimpleFe Benchmark Model
(Du 2023, Geosci. Model Dev., 16, 5865-5894).
==========================================================================

Overview:
---------
This module defines all chemical species modeled in the SimpleFe diagenesis benchmark.
Species are categorized into three physical phases:
  1. Particulate Solids (`phase="solid"`):
     Bound within the sediment matrix; concentrations are expressed per unit volume of
     solid sediment (mol/m^3_solid). They move solely by sediment accumulation (burial)
     and macrofaunal bioturbation mixing.
  2. Dissolved Porewater Solutes (`phase="liquid"`):
     Dissolved in the pore fluid; concentrations are expressed per unit volume of
     porewater (mol/m^3_pw = mM). They move by molecular diffusion, advection, and
     bio-irrigation.
  3. Mixed Pools (`phase="mixed"`):
     Lumped species like total ferrous iron (Fe2_total) that partition into a dissolved
     fraction and a surface-adsorbed fraction.

Non-transported Diagnostic & Speciation Species:
-----------------------------------------------
Species with `include: False` are not integrated as primary PDE state variables by the
solver. Instead, they represent dynamic equilibrium speciation products (e.g. Fe2_liq,
HS, H2S, CO3, HCO3, Hplus) that are updated algebraically from their parent pools
during reaction and Jacobian evaluation and exported to the output DataFrame and state NPZ.

Unit Conversion Reference:
--------------------------
  • Solids: SedTrace stores solids in mmol/cm^3_solid.
    Conversion to Mantis SI: 1 mmol/cm^3_solid = 1000 mol/m^3_solid.
    Conversion to dry wt%: wt% = C_solid [mol/m^3] * MW [g/mol] / (10 * rho_s [g/cm^3] * 1000).
  • Dissolved: SedTrace stores solutes in mmol/cm^3_pw.
    Conversion to Mantis SI: 1 mmol/cm^3_pw = 1000 mol/m^3_pw = 1000 mM = 10^6 uM.
"""

species = {
    # -------------------------------------------------------------------------
    # Particulate Solid Phases (mol/m^3 of solid sediment)
    # -------------------------------------------------------------------------
    "POC": {
        "phase": "solid",
        "include": True,
        "formula": "CH2O",
        "diffusion": None,
        # Surface flux: FPOC0 = 0.31 mmol/(cm^2 * yr) = 3.1 mol/(m^2 * yr)
        # Translates to initial steady-state surface concentration: J_solid / w_top ~ 7660 mol/m^3 (~3.54 wt%)
        "bc_top": {"Robin": "0.31 mmol/(cm^2 * yr)"},
        "prefill": True,
    },
    "FeOOH": {
        "phase": "solid",
        "include": True,
        "formula": "FeOOH",
        "diffusion": None,
        # Surface flux: FFeOOH0 = 0.023 mmol/(cm^2 * yr) = 0.23 mol/(m^2 * yr)
        # Translates to initial surface concentration: J_solid / w_top ~ 568 mol/m^3 (~0.10 wt% Fe)
        "bc_top": {"Robin": "0.023 mmol/(cm^2 * yr)"},
        "prefill": True,
    },
    "FeS": {
        "phase": "solid",
        "include": True,
        "formula": "FeS",
        "diffusion": None,
        # Zero depositional flux at sediment surface; generated authigenically by precipitation
        "bc_top": {"Robin": "0.0 mmol/(cm^2 * yr)"},
        "prefill": False,
    },

    # -------------------------------------------------------------------------
    # Dissolved Solutes (Porewater Transport, mol/m^3 = mM)
    # -------------------------------------------------------------------------
    "SO4": {
        "phase": "liquid",
        "include": True,
        "formula": "SO4",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 4.99, "m1": 0.231},  # Boudreau (1997) temperature parameters
        # In Du (2023), SO4 has a Robin boundary condition with a 0.05 cm diffusive boundary layer (delta = 0.05 cm)
        "bc_top": {"Dirichlet": 28, "dbl": "0.05 cm"},
        "prefill": True,
    },
    "TS2": {
        "phase": "liquid",
        "include": True,
        "formula": "H2S",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 10.4, "m1": 0.273},
        "bc_top": {"Dirichlet": 0},
        "prefill": False,
    },
    "Fe2_total": {
        "phase": "mixed",
        "basis": "porewater",
        "include": True,
        "formula": "Fe2",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 3.31, "m1": 0.150},
        # Surface boundary condition: in Du (2023), TFe_dis0 = 2.95e-5 mM and FTFe_ads0 = 0.
        # The total lumped iron at the SWI is TFe0 = TFe_dis0 * (1 + K_ads * POC0 * dstopw) ~ 5.46e-5 mM.
        "bc_top": {"Dirichlet": "5.4574e-5 mmol/l"},
        "prefill": True,
    },

    # -------------------------------------------------------------------------
    # Carbonate System & Dynamic Porewater pH
    # -------------------------------------------------------------------------
    "DIC": {
        "phase": "liquid",
        "include": True,
        "formula": "CO2",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 5.06, "m1": 0.275},
        "bc_top": {"Dirichlet": 2.345},
        "prefill": True,
    },
    "TA": {
        "phase": "liquid",
        "include": True,
        "formula": "HCO3",
        "diff_coeff_function": "diff_coeff",
        "diffusion": "DIC",                          # Inherit DIC diffusion profile
        "bc_top": {"Dirichlet": 2.329},
        "prefill": True,
    },
    "Hplus": {
        "phase": "liquid",
        "include": False,
        "formula": "H",
        "bc_top": {"Dirichlet": 7.59},
    },

    # -------------------------------------------------------------------------
    # Diagnostic Solutes
    # --------------------------------------------------------------------------
    "Ca2": {
        "phase": "liquid",
        "include": False,
        "formula": "Ca",
        "bc_top": {"Dirichlet": 10.3},    # 10.3 mM typical seawater Ca2+
    },
}


# =============================================================================
# Unified Speciation & Phase Partitioning
# =============================================================================
  
speciation = {
    # -------------------------------------------------------------------------
    # Marine Carbonate System Speciation (DIC + TA -> CO3, HCO3, Hplus, pH)
    # -------------------------------------------------------------------------
    "carbonate_system": {
        "type": "carbonate_system",
        "primary": {"dic": "DIC", "ta": "TA"},
        "species": ["Hplus", "CO3", "HCO3"],
        "pressure": 50.0,       # 500 m water depth (~50 bar hydrostatic pressure)
        "salinity": 35.0,       # 35 PSU bottom water salinity
        "adjust_omega": False,  # Retain exact boundary values without forcing calcite equilibrium
        "ph_scale": "free",
        "controls": ["DIC", "pH"],
    },

    # -------------------------------------------------------------------------
    # Ferrous Iron Mixed-Phase Partitioning & Ligand Speciation
    # -------------------------------------------------------------------------
    "fe2_speciation": {
        "type": "complexation",
        "primary": "Fe2_total",
        "constants": {
            "KFe_ads": 0.001,      # (mol/m^3)^-1 = 1.0 (mmol/cm^3)^-1 for POC adsorption
            "beta_Cl": 0.4292,     # 10^-0.12 * 0.5658 M Cl- inorganic complex
            "beta_SO4": 0.00912,   # 10^0.96 / 1000 (mol/m^3)^-1 for FeSO4 complex
            "beta_CO3": 4.4668,    # 10^3.65 / 1000 (mol/m^3)^-1 for FeCO3 complex
            "beta_FeHS": 251.1886, # 10^5.4 / 1000 (mol/m^3)^-1 for soluble FeHS+ complex
        },
        "species": {
            "Fe2_dis": {
                "phase": "liquid",
                "fraction": "1.0 / (1.0 + KFe_ads * POC * dstopw)",
            },
            "Fe2_liq": {
                "phase": "liquid",
                "fraction": "1.0 / ((1.0 + KFe_ads * POC * dstopw) * (1.0 + beta_Cl + beta_SO4 * SO4 + beta_CO3 * CO3 + beta_FeHS * HS))",
            },
        },
    },

    # -------------------------------------------------------------------------
    # Dissolved Sulfide Acid-Base Speciation (TS2 -> HS- + H2S)
    # -------------------------------------------------------------------------
    "sulfide_acid_base": {
        "type": "acid_base",
        "primary": "TS2",
        "constants": {"pKa1": 6.8912},
        "species": {
            "HS": {
                "fraction": "Ka1 / (Hplus + Ka1)",
                "include_in_TA": 1.0,
            },
            "H2S": {
                "fraction": "Hplus / (Hplus + Ka1)",
            },
        },
    },
}
