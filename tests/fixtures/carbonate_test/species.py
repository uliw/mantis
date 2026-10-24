"""species.py

Declarative definition of all chemical species in the diagenetic system with carbonate chemistry.
"""

import numpy as np


# -----------------------------------------------------------------------------
# Isotope Reference Standards
# -----------------------------------------------------------------------------
isotope_standards: dict = {
    "d34S": {
        "element": "S",
        "standard": "VCDT",
        "ratio": 0.044162589,  # 34S / 32S
        "unit": "mUr",
    },
}


def boudreau_o2(T_profile, phi):
    """Temperature porosity corrected O2 molecular diffusion."""
    phi_val = phi.value if hasattr(phi, "value") else phi
    return (
        (0.2604 + 0.006363 * (T_profile + 273.15))
        * 1e-9
        / (1.0 - np.log(phi_val**2))
    )


species = {
    # -------------------------------------------------------------------------
    # Dissolved Liquid Species (Porewater Transport)
    # -------------------------------------------------------------------------
    "SO4": {
        "phase": "liquid",
        "include": True,
        "formula": "SO4",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 4.88, "m1": 0.232},
        "isotopes": {"SO4_S32": "S"},
        "bc_top": {"Dirichlet": 28.2, "d34S": 21.0},
        "prefill": True,
    },
    "TS2": {
        "phase": "liquid",
        "include": True,
        "formula": "H2S",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 10.4, "m1": 0.273},
        "isotopes": {"TS2_S32": "S"},
        "bc_top": {"Dirichlet": 0.0, "d34S": 0.0},
    },
    "O2": {
        "phase": "liquid",
        "include": True,
        "formula": "O2",
        "diff_coeff_function": "boudreau_o2",
        "diffusion": None,
        "bc_top": {"Dirichlet": 0.28},
        "prefill": True,
    },

    # -------------------------------------------------------------------------
    # Diagnostic Speciation Solutes
    # -------------------------------------------------------------------------
    "Hplus": {
        "phase": "liquid",
        "include": False,
        "formula": "H",
        "bc_top": {"Dirichlet": None},  # ph is derived from TA and DIC and set to be in Equilibrium with CaCO3
    },

    # -------------------------------------------------------------------------
    # Mixed / Lumped Species (Porewater + Sorbed Solid)
    # -------------------------------------------------------------------------
    "Fe2_total": {
        "phase": "mixed",
        "include": True,
        "formula": "Fe2",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 27.7, "m1": 1.0},
        "bc_top": {"Dirichlet": 0.0},
    },

    # -------------------------------------------------------------------------
    # Solid Mineral & Particulate Species (Burial + Bioturbation)
    # -------------------------------------------------------------------------
    "POC_fast": {
        "phase": "solid",
        "include": True,
        "formula": "CH2O6",
        "diffusion": None,
        "bc_top": {"Robin": "365 umol/(cm^2 * yr)"},
        "prefill": False,
    },
    "POC_slow": {
        "phase": "solid",
        "include": True,
        "formula": "CH2O6",
        "diffusion": None,
        "bc_top": {"Robin": "183 umol/(cm^2 * yr)"},
        "prefill": False,
    },
    "Fe3": {
        "phase": "solid",
        "include": True,
        "formula": "FeOOH",
        "diffusion": None,
        "bc_top": {"Robin": "12 umol/(cm^2 * yr)"},
        "prefill": False,
    },
    "FeS": {
        "phase": "solid",
        "include": True,
        "formula": "FeS",
        "diffusion": None,
        "isotopes": {"FeS_S32": "S"},
        "bc_top": {"Robin": 0.0, "d34S": 0.0},
    },
    "FeS2": {
        "phase": "solid",
        "include": True,
        "formula": "FeS2",
        "diffusion": None,
        "isotopes": {"FeS2_S32": "S"},
        "stoichiometry": 2.0,
        "bc_top": {"Robin": 0.0, "d34S": 0.0},
    },
    "S0": {
        "phase": "solid",
        "include": True,
        "formula": "S",
        "diffusion": None,
        "isotopes": {"S0_S32": "S"},
        "bc_top": {"Robin": 0.0, "d34S": 8.0},
    },

    # External / Reaction Balance Species (Excluded from transport)
    "CO2": {"phase": "liquid", "include": False, "formula": "CO2"},
    "H2O": {"phase": "liquid", "include": False, "formula": "H2O"},

    # -------------------------------------------------------------------------
    # Carbonate System Species
    # -------------------------------------------------------------------------
    "DIC": {
        "phase": "liquid",
        "include": True,
        "formula": "CO2",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 5.06, "m1": 0.278},
        "bc_top": {"Dirichlet": 2.2},
        "prefill": True,
    },
    "TA": {
        "phase": "liquid",
        "include": True,
        "formula": "HCO3",
        "diff_coeff_function": "diff_coeff",
        "diffusion": "DIC",
        "bc_top": {"Dirichlet": 2.4},
        "prefill": True,
    },
    "Ca2": {
        "phase": "liquid",
        "include": True,
        "formula": "Ca",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 4.19, "m1": 0.231},
        "bc_top": {"Dirichlet": 10.3},
        "prefill": True,
    },
    "CaCO3": {
        "phase": "solid",
        "include": True,
        "formula": "CaCO3",
        "diffusion": None,
        "bc_top": {"Robin": "20 umol/(cm^2 * yr)"},
        "prefill": True,
    },
}


# =============================================================================
# Unified Speciation & Phase Partitioning
# =============================================================================
speciation: Dict[str, Dict[str, Any]] = {
    "carbonate_system": {
        "type": "carbonate_system",
        "primary": {"dic": "DIC", "ta": "TA"},
        "species": ["Hplus", "CO3", "HCO3"],
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
                "isotopes": {"HS_S32": "S"},
            },
            "H2S": {
                "phase": "liquid",
                "fraction": "Hplus / (Hplus + Ka1)",
                "isotopes": {"H2S_S32": "S"},
            },
        },
    },
    "fe2_partitioning": {
        "type": "partition",
        "primary": "Fe2_total",
        "fractions": {
            "dissolved": "1.0 / (1.0 + k.Fe2_p_eq)",
            "sorbed": "k.Fe2_p_eq / (1.0 + k.Fe2_p_eq)",
        },
        "species": {
            "Fe2_liq": {"phase": "liquid", "alias_of": "Fe2_total.dissolved"},
            "Fe2_sorbed": {"phase": "solid", "alias_of": "Fe2_total.sorbed"},
        },
    },
}
