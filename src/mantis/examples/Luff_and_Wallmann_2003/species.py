"""species.py

Declarative definition of all chemical species for the Luff & Wallmann (2003)
Hydrate Ridge cold vent diagenetic model.

Reference:
----------
Luff R. and Wallmann K. (2003) Fluid flow, methane fluxes, carbonate precipitation
and biogeochemical turnover in gas hydrate-bearing sediments at Hydrate Ridge,
Cascadia Margin: Numerical modeling and mass balances.
Geochimica et Cosmochimica Acta 67(18), 3403-3421.
"""

import numpy as np

species = {
    # -------------------------------------------------------------------------
    # Dissolved Liquid Species (Porewater Transport with Upward Advection)
    # -------------------------------------------------------------------------
    "SO4": {
        "phase": "liquid",
        "include": True,
        "formula": "SO4",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 4.88, "m1": 0.232},
        "bc_top": {"Dirichlet": 27.8},  # mM from salinity
        "bc_bot": {"Dirichlet": 0.0},
        "prefill": True,
    },
    "CH4": {
        "phase": "liquid",
        "include": True,
        "formula": "CH4",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 8.0, "m1": 0.39},  # Boudreau 1996 for CH4
        "bc_top": {"Dirichlet": 0.0},
        "bc_bot": {"Dirichlet": 68.0},  # Hydrate equilibrium at 78 atm, 4.2 C
        "prefill": False,
    },
    "TS2": {
        "phase": "liquid",
        "include": True,
        "formula": "H2S",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 10.4, "m1": 0.273},
        "bc_top": {
            "Dirichlet": 3.5
        },  # Upper sediment/mat boundary (Table 1: H2S 1.0 + HS- 2.5 mM)
        "bc_bot": {
            "Dirichlet": 20.0
        },  # Bottom boundary (Table 1: H2S 4.1 + HS- 15.9 mM)
        "prefill": True,
    },
    "Ca2": {
        "phase": "liquid",
        "include": True,
        "formula": "Ca",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 4.19, "m1": 0.231},
        "bc_top": {"Dirichlet": 10.1},  # mM from salinity
        "bc_bot": {"Dirichlet": 3.5},
        "prefill": True,
    },
    "DIC": {
        "phase": "liquid",
        "include": True,
        "formula": "CO2",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 5.06, "m1": 0.278},
        "bc_top": {"Dirichlet": 2.53},  # CO2(0.3) + HCO3-(2.2) + CO32-(0.03) = 2.53 mM
        "bc_bot": {
            "Dirichlet": 16.55
        },  # CO2(0.75) + HCO3-(15.6) + CO32-(0.2) = 16.55 mM
        "prefill": True,
    },
    "TA": {
        "phase": "liquid",
        "include": True,
        "formula": "HCO3",
        "diff_coeff_function": "diff_coeff",
        "diffusion": "DIC",
        "bc_top": {"Dirichlet": 4.76},  # Table 1: TA_carb (2.26) + HS- (2.5) = 4.76 mM
        "bc_bot": {
            "Dirichlet": 29.394
        },  # Calcite equilibrium at bottom: CA (15.83) + HS- (13.57) = 29.4 mM (Table 1: 29.4 mM)
        "prefill": True,
    },
    "TPO4": {
        "phase": "liquid",
        "include": True,
        "formula": "H3PO4",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 3.54, "m1": 0.187},
        "bc_top": {"Dirichlet": 0.047},  # 47.3 uM
        "bc_bot": {"Dirichlet": 0.012},  # 11.7 uM
        "prefill": True,
    },
    "TNH4": {
        "phase": "liquid",
        "include": True,
        "formula": "NH4",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 10.7, "m1": 0.365},
        "bc_top": {"Dirichlet": 0.117},  # 117 uM
        "bc_bot": {"Dirichlet": 0.210},  # ~210 uM
        "prefill": True,
    },
    # -------------------------------------------------------------------------
    # Diagnostic Speciation Solutes (Not transported dynamically)
    # -------------------------------------------------------------------------
    "Hplus": {
        "phase": "liquid",
        "include": False,
        "formula": "H",
        "bc_top": {"Dirichlet": 7.5},
        "prefill": True,
    },
    "CO2": {"phase": "liquid", "include": False, "formula": "CO2"},
    "H2O": {"phase": "liquid", "include": False, "formula": "H2O"},
    # -------------------------------------------------------------------------
    # Solid Mineral & Particulate Species (Burial Only)
    # -------------------------------------------------------------------------
    "POC_fast": {
        "phase": "solid",
        "include": True,
        "formula": "CH2O6",
        "diffusion": None,
        "bc_top": {"Robin": "55 umol/(cm^2 * yr)"},
        # Steady-state exponential decay: C0 * exp(-(k/w)*z)
        # C0 = 55 / (w * (1 - phi)) = 7692.3 mol/m^3_solid, k/w = 0.2 / 0.0275 = 727.27 m^-1
        "init": lambda z: 7692.3 * np.exp(-727.27 * z),
        "prefill": False,
    },
    "POC_slow": {
        "phase": "solid",
        "include": True,
        "formula": "CH2O6",
        "diffusion": None,
        "bc_top": {"Robin": "26 umol/(cm^2 * yr)"},
        # C0 = 26 / (w * (1 - phi)) = 3636.4 mol/m^3_solid, k/w = 3e-4 / 0.0275 = 1.0909 m^-1
        "init": lambda z: 3636.4 * np.exp(-1.0909 * z),
        "prefill": False,
    },
    "Calcite": {
        "phase": "solid",
        "include": True,
        "formula": "CaCO3",
        "diffusion": None,
        "bc_top": {"Robin": "1.0 umol/(cm^2 * yr)"},
        # Background sediment carbonate (~18 wt% = 4500 mol/m^3_solid)
        "init": lambda z: np.full_like(z, 4500.0),
        "prefill": False,
    },
    "Aragonite": {
        "phase": "solid",
        "include": True,
        "formula": "CaCO3",
        "diffusion": None,
        "bc_top": {"Robin": "1.0 umol/(cm^2 * yr)"},
        # Pure authigenic mineral phase, initialized to 0.0
        "init": 0.0,
        "prefill": False,
    },
}


# =============================================================================
# Unified Speciation & Phase Partitioning
# =============================================================================
speciation = {
    "carbonate_system": {
        "type": "carbonate_system",
        "primary": {"dic": "DIC", "ta": "TA"},
        "species": ["Hplus", "CO3", "HCO3"],
        "pressure": 79.0,  # 500 m water depth (~50 bar hydrostatic pressure)
        "salinity": 34.3,  # 35 PSU bottom water salinity
        "adjust_omega": False,  # Retain exact boundary values without forcing calcite equilibrium
        "ph_scale": "free",
        "controls": ["DIC", "TA"],
        "DIC": 2.53,
        "TA": 4.76,
    },
    "sulfide_acid_base": {
        "type": "acid_base",
        "primary": "TS2",
        "constants": {"pKa1": 7.0},
        "species": {
            "HS": {
                "phase": "liquid",
                "fraction": "Ka1 / (Hplus + Ka1)",
                "include_in_TA": 1.0,
            },
            "H2S": {
                "phase": "liquid",
                "fraction": "Hplus / (Hplus + Ka1)",
            },
        },
    },
}
