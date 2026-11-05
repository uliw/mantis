"""species.py

Declarative definition of all chemical species in the diagenetic system.

Template Format & Structure:
----------------------------
Each key in the `species` dictionary represents a chemical species name:
- phase: 'liquid', 'solid', or 'mixed'.
    * 'liquid': Dissolved species in porewater (scaled by porosity phi).
    * 'solid': Mineral or particulate species (scaled by 1 - phi).
    * 'mixed': Lumped pool (e.g. Fe2_total) partitioned into porewater and sorbed phase.
- include: bool.
    * True: Active transport PDE variable solved by FiPy / Newton solver.
    * False: Equilibrium, diagnostic, or fixed external species.
- formula: Chemical formula string (e.g., 'FeS2', 'SO4', 'CH2O6').
- equilibrium_of: (Optional) Parent species for fast dynamic speciation (e.g., 'TS2').
- fraction: (Optional) Analytical expression or constant defining the fraction of parent.
- partition: (Optional for 'mixed') Fractions partitioned into dissolved and sorbed phases.
- diffusion: (Optional) Recipe for molecular diffusion coefficient calculation:
    * dict with {'type': 'diff_coeff', 'm0': float, 'm1': float} (Boudreau 1996)
    * dict with {'type': 'boudreau_o2'} (temperature & Archie porosity corrected O2 diffusion)
    * str pointing to another species name to inherit its diffusion coefficient.
    * None for solid species (which only transport via sediment burial and bioturbation).
- bc: (Optional) Default boundary condition dict:
- bc: (Optional) Default boundary condition dict:
    * {'top/bottom': float or pint-string, 'type': 'Dirichlet', 'Robin' or 'Neumann'}
    * unless specified, bottom bc defaults to Neumann 
    * For Hplus: {'pH': 7.5} sets porewater pH as the single source of truth.
- isotopes: (Optional){"isotope_species_name": "element"}
- isotope_of: (Optional) Base species name if this entry is an isotope companion.
- delta_threshold: (Optional) Minimum concentration threshold for calculating delta values (mmol/L).
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


"""
Most species are covered by the diff_coeff function with the m0, m1 parameter listed in
Table 4.8 in Boudreau 1996.  If no parameters are availble, please provide any custom
functions here.  Note that diffusions coefficients are only needed for species that are
par of the reaction-transport martrix (i.e., where include = True).
"""
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
        "bc_top": {"Dirichlet": 28.2, "d34S": 21},
        "prefill": True,
    },
    "TS2": {
        "phase": "liquid",
        "include": True,
        "formula": "H2S",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 10.4, "m1": 0.273},
        "isotopes": {"TS2_S32": "S"},
        "bc_top": {"Neumann": 0.0},
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
        "bc_top": {"Dirichlet": 7.5},  # Single source of truth for porewater pH
        "prefill": True,
    },

    # -------------------------------------------------------------------------
    # Mixed / Lumped Species (Porewater + Sorbed Solid)
    # -------------------------------------------------------------------------
    "Fe2_total": {
        "phase": "mixed",
        "include": True,
        "formula": "Fe2",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 3.31, "m1": 0.15},
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
    },
    "POC_slow": {
        "phase": "solid",
        "include": True,
        "formula": "CH2O6",
        "diffusion": None,
        "bc_top": {"Robin": "183 umol/(cm^2 * yr)"},
    },
    "Fe3": {
        "phase": "solid",
        "include": True,
        "formula": "FeOOH",
        "diffusion": None,
        "bc_top": {"Robin": "12 umol/(cm^2 * yr)"},
    },
    "FeS": {
        "phase": "solid",
        "include": True,
        "formula": "FeS",
        "diffusion": None,
        "isotopes": {"FeS_S32": "S"},
        "bc_top": {"Robin": 0.0, "d34S": 0},
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
    }, # species we need for the carbonate system
    "DIC": {
        "phase": "liquid",
        "include": True,
        "formula": "CO2",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 5.06, "m1": 0.278},
        "bc_top": {"Dirichlet": 2.53},  # CO2(0.3) + HCO3-(2.2) + CO32-(0.03) = 2.53 mM
        "bc_bot": {"Dirichlet": 16.55}, # CO2(0.75) + HCO3-(15.6) + CO32-(0.2) = 16.55 mM
        "prefill": True,
    },
    "TA": {
        "phase": "liquid",
        "include": True,
        "formula": "HCO3",
        "diff_coeff_function": "diff_coeff",
        "diffusion": "DIC",
        "bc_top": {"Dirichlet": 4.76},  # Table 1: TA_carb (2.26) + HS- (2.5) = 4.76 mM
        "bc_bot": {"Dirichlet": 29.394},  # Calcite equilibrium at bottom: CA (15.83) + HS- (13.57) = 29.4 mM (Table 1: 29.4 mM)
        "prefill": True,
    },
     "Calcite": {
        "phase": "solid",
        "include": True,
        "formula": "CaCO3",
        "diffusion": None,
        "bc_top": {"Robin": "100.0 umol/(cm^2 * yr)"},
        # Background sediment carbonate (~18 wt% = 4500 mol/m^3_solid)
        #"init": lambda z: np.full_like(z, 4500.0),
        "prefill": True,
    },
    # External / Reaction Balance Species (Excluded from transport)
    "CO2": {"phase": "liquid", "include": False, "formula": "CO2"},
    "H2O": {"phase": "liquid", "include": False, "formula": "H2O"},
}


# =============================================================================
# Unified Speciation & Phase Partitioning
# =============================================================================
speciation = {
    "sulfide_acid_base": {
        "type": "acid_base",
        "primary": "TS2",
        "constants": {"pKa1": 7.0},
        "equilibrium_fractionation": {
            "standard": "d34S",
            "alpha": {"H2S_S32_to_HS_S32": 0.99991542},
        },
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
    }, # carbonate speciation
    "carbonate_system": {
        "type": "carbonate_system",
        "primary": {"dic": "DIC", "ta": "TA"},
        "species": ["Hplus", "CO3", "HCO3"],
    },
}



