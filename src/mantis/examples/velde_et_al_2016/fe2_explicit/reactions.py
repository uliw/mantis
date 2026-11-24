"""Mantis-RTM: A Python library for reaction-transport modeling

Copyright (C), 2026 - Ulrich G. Wortmann

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
"""fes_reactions.py

Declarative definition of diagenetic chemical reactions and kinetic limiters.

Template Structure & Design:
----------------------------
1. `limiters`:
   Dictionary mapping unique limiter names to their mathematical expressions.
   Moving limiters here keeps all kinetic terms and thresholds adjacent to the
   reactions they regulate, eliminating detached limiter scripts.

2. `reactions`:
   List of chemical reaction specifications. Each entry contains:
   - reaction: Stoichiometric equation string (e.g., '1 FeS + 1 HS -> 1 FeS2').
     Reactants and products are parsed and verified against species.py.
     Dynamic equilibrium species (like HS and Fe2_liq) are recognized automatically.
   - reaction_name: Unique identifier for reporting and logging.
   - k_value_name: Rate constant identifier defined in kinetic_constants_slow.py.
     Can be a string (e.g. 'POC_fast') or a dict for multi-branch reactions:
     {'precipitation': 'FeS_isp', 'dissolution': 'FeS_isd', 'equilibrium': 'FeS_sp'}.
   - limiters: Dict or list of limiters applied to this reaction.
   - isotope_species: (Optional) Explicit mapping between base and 32S companion species.
   - isotope_fractionation: (Optional) Isotope fractionation factor and limiter.
"""


limiters = {
    # Base implicit limiters: L_imp = 1.0 / (C + K)
    # The parser automatically derives explicit (C * L_imp) and inhibition (K * L_imp = 1 - C * L_imp)
    "O2": "1.0 / (c.O2 + K_O2)",
    "O2_H2S_tot": "1.0 / (c.O2 + K_O2_H2S_tot)",
    "SO4": "1.0 / (c.SO4 + K_SO4)",
    "Fe3": "1.0 / (np.maximum(c.Fe3, 0.0) + K_Fe3)",
    "Fe3_diss_red": "1.0 / (np.maximum(c.Fe3, 0.0) + K_Fe3_diss_red)",
}

reactions = [
    # -------------------------------------------------------------------------
    # Organic Carbon Primary Respiration
    # -------------------------------------------------------------------------
    {
        "reaction_name": "aerobic_respiration_fast",
        "reaction": "POC_fast + O2 -> CO2",
        "k_value_name": "POC_fast",
        "limiters": ["O2"],
    },
    {
        "reaction_name": "aerobic_respiration_slow",
        "reaction": "POC_slow + O2 -> CO2",
        "k_value_name": "POC_slow",
        "limiters": ["O2"],
    },
    {
        "reaction_name": "dissimilatory_iron_reduction_fast",
        "reaction": "POC_fast + 4 Fe3 -> 4 Fe2_liq",
        "k_value_name": "POC_fast",
        "limiters": ["Fe3_diss_red"],
        "inhibitors": ["O2"],
    },
    {
        "reaction_name": "dissimilatory_iron_reduction_slow",
        "reaction": "POC_slow + 4 Fe3 -> 4 Fe2_liq",
        "k_value_name": "POC_slow",
        "limiters": ["Fe3_diss_red"],
        "inhibitors": ["O2"],
    },
    {
        "reaction_name": "sulfate_reduction_fast",
        "reaction": "2 POC_fast + SO4 -> H2S_tot",
        "k_value_name": "POC_fast",
        "limiters": ["SO4"],
        "inhibitors": ["O2", "Fe3_diss_red"],
    },
    {
        "reaction_name": "sulfate_reduction_slow",
        "reaction": "2 POC_slow + SO4 -> H2S_tot",
        "k_value_name": "POC_slow",
        "limiters": ["SO4"],
        "inhibitors": ["O2", "Fe3_diss_red"],
        },

    # -------------------------------------------------------------------------
    # Secondary Redox & Sulfide Oxidation
    # -------------------------------------------------------------------------
    {
        "reaction_name": "hs_oxidation",
        "reaction": "1 HS + 2 O2 -> 1 SO4",
        "k_value_name": "H2S_tot_O2",
        # "limiters": ["O2_H2S_tot"],
    },
    {
        "reaction_name": "sulfide_mediated_iron_reduction_velde",
        "reaction": "HS + 8 Fe3 -> SO4 + 8 Fe2_liq",
        "k_value_name": "Fe3_hs",
        # "limiters": ["Fe3"],
        "inhibitors": ["O2"],
    },
    {
        "reaction_name": "Fe2_liq_oxidation",
        "reaction": "4 Fe2_liq + O2 -> 4 Fe3",
        "k_value_name": "Fe2_O2",
    },
    {
        "reaction_name": "Fe2_sorbed_oxidation",
        "reaction": "4 Fe2_sorbed + O2 -> 4 Fe3",
        "k_value_name": "Fe2_O2",
    },
    {
        "reaction_name": "FeS_oxidation",
        "reaction": "4 FeS + 9 O2 -> 4 Fe3 + 4 SO4",
        "k_value_name": "FeS_O2",
    },

    # -------------------------------------------------------------------------
    # Explicit Fast Fe2 Sorption / Desorption
    # -------------------------------------------------------------------------
    {
        "reaction_name": "Fe2_adsorption",
        "reaction": "Fe2_liq -> Fe2_sorbed",
        "k_value_name": "Fe2_ads",
    },
    {
        "reaction_name": "Fe2_desorption",
        "reaction": "Fe2_sorbed -> Fe2_liq",
        "k_value_name": "Fe2_des",
    },

    # -------------------------------------------------------------------------
    # Mineral Precipitation & Dissolution (Smooth Transition)
    # -------------------------------------------------------------------------
    {
        "reaction_name": "FeS_precipitation_dissolution",
        "reaction": "Fe2_liq + HS -> FeS",
        "k_value_name": {
            "precipitation": "FeS_isp",
            "dissolution": "FeS_isd",
            "equilibrium": "FeS_sp",
        },
    },
]
