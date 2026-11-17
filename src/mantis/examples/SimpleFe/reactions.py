"""simplefe_reactions.py

Chemical Reaction Network & Kinetic Limiters for the SimpleFe Benchmark Model
(Du 2023, Geosci. Model Dev., 16, 5865-5894).
=============================================================================

Reaction Network Overview:
--------------------------
The SimpleFe benchmark couples 4 primary biogeochemical reaction processes:
  1. Dissimilatory Iron Reduction (DIR; RFeOOHPOC):
     CH2O + 4 FeOOH + 8 H+ -> CO2 + 4 Fe2+ + 7 H2O
     In Mantis's DIC/TA conservation formulation:
     POC + 4 FeOOH -> DIC + 4 Fe2_total + 8 TA
     Rate: R1 = k_POC * [POC] * [FeOOH] / (K_FeOOH + [FeOOH])

  2. Organoclastic Sulfate Reduction (OSR; RSO4POC):
     CH2O + 0.5 SO4^2- + H+ -> CO2 + 0.5 H2S + H2O
     Normalized per 1.0 mole of SO4^2-:
     2 POC + SO4 -> 2 DIC + TS2 + 2 TA
     Rate: R2 = k_POC * [POC] * [SO4] / (K_SO4 + [SO4]) * K_FeOOH / (K_FeOOH + [FeOOH])
     (Note the Monod limitation by sulfate and inhibition by reactive iron oxide).

  3. Abiotic Sulfide Oxidation by Iron Oxide (RFeOOHH2S):
     2 FeOOH + H2S + 4 H+ -> 2 Fe2+ + S0 + 4 H2O
     In Mantis's DIC/TA conservation formulation:
     2 FeOOH + TS2 -> 2 Fe2_total + 4 TA
     Rate: R3 = k_FeOOHH2S * [FeOOH] * [TS2]

  4. Iron Monosulfide Precipitation & Reversible Dissolution (RFeS_pre):
     Fe2+ + HS- -> FeS + H+  (consumes 1 TA)
     Precipitation is governed by the saturation state:
         Omega = [Fe2_free] * [HS-] / ([H+] * Ksp)
     Rate: R4 = k_isp * max(Omega - 1, 0)
     Only free uncomplexed Fe2+ (Fe2_liq) drives precipitation, avoiding unphysical
     hyper-precipitation of complexed iron (e.g. FeHS+).

Alkalinity and Stoichiometric Mapping:
--------------------------------------
Mantis models porewater pH by tracking Total Dissolved Inorganic Carbon (DIC) and
Total Alkalinity (TA). Each reaction's proton production/consumption is mapped to TA:
  • Consuming 1 mole of H+ produces 1 mole of TA (e.g. +8 TA in DIR, +4 TA in sulfide oxidation).
  • Producing 1 mole of H+ consumes 1 mole of TA (e.g. -1 TA in FeS precipitation).
"""

from typing import Dict, List, Any


limiters = {
    # Base implicit limiter definitions: 1.0 / (C + K)
    # Mantis automatically derives:
    #   - Monod limitation term: C / (C + K) = C * (1.0 / (C + K))
    #   - Inhibition term: K / (C + K) = K * (1.0 / (C + K))
    "FeOOH": "1.0 / (np.maximum(c.FeOOH, 0.0) + K_FeOOH)",
    "SO4": "1.0 / (c.SO4 + K_SO4)",
}

reactions: List[Dict[str, Any]] = [
    # -------------------------------------------------------------------------
    # 1. Dissimilatory Iron Reduction (DIR; RFeOOHPOC)
    # CH2O + 4 FeOOH + 8 H+ -> CO2 + 4 Fe2+ + 7 H2O
    # -------------------------------------------------------------------------
    {
        "reaction_name": "iron_reduction_poc",
        "reaction": "POC + 4 FeOOH -> DIC + 4 Fe2_total + 8 TA",
        "k_value_name": "k_POC",
        "ref_species": "POC",
        "limiters": ["FeOOH"],
    },

    # -------------------------------------------------------------------------
    # 2. Organoclastic Sulfate Reduction (OSR; RSO4POC)
    # 2 CH2O + SO4^2- + 2 H+ -> 2 CO2 + H2S + 2 H2O
    # -------------------------------------------------------------------------
    {
        "reaction_name": "sulfate_reduction_poc",
        "reaction": "2 POC + SO4 -> 2 DIC + TS2 + 2 TA",
        "k_value_name": "k_POC",
        "ref_species": "POC",
        "limiters": ["SO4"],
        "inhibitors": ["FeOOH"],
    },

    # -------------------------------------------------------------------------
    # 3. Sulfide Oxidation by Iron Oxide (RFeOOHH2S)
    # 2 FeOOH + H2S + 4 H+ -> 2 Fe2+ + S0 + 4 H2O
    # -------------------------------------------------------------------------
    {
        "reaction_name": "iron_oxide_sulfide_oxidation",
        "reaction": "2 FeOOH + TS2 -> 2 Fe2_total + 4 TA",
        "k_value_name": "k_FeOOHH2S",
        "ref_species": "TS2",
    },

    # -------------------------------------------------------------------------
    # 4. Iron Monosulfide Precipitation & Dissolution (RFeS_pre)
    # Fe2+ + HS- <-> FeS + H+
    # Driven by the saturation state: Omega = [Fe2_liq][HS] / ([H+] * Ksp)
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
