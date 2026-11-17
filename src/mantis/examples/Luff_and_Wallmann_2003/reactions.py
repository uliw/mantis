"""reactions.py

Reaction network and kinetic limiters for the Luff & Wallmann (2003) Hydrate Ridge model.
"""

reaction_parameters = {
    # Redfield stoichiometry: C:N:P = 106:16:1
    "C_N_ratio": 106.0 / 16.0,
    "C_P_ratio": 106.0 / 1.0,
}
fractionation_factors = reaction_parameters

limiters = {
    # Sulfate limitation and inhibition
    "SO4": "1.0 / (c.SO4 + K_SO4)",

    # Calcite saturation and rate limiters (Equation 4: RPPT = kPPT * (Omega - 1))
    "calcite_precip_limiter": "np.maximum((c.Ca2 * c.CO3) / (mp.Ksp_calcite + 1e-30) - 1.0, 0.0)",
    "calcite_diss_limiter": "np.maximum(1.0 - (c.Ca2 * c.CO3) / (mp.Ksp_calcite + 1e-30), 0.0)",

    # Aragonite saturation and rate limiters
    # Ksp_arag = 0.8168 (mol/m^3)^2 at 4.2 C, 34.3 PSU, 78 atm
    "aragonite_precip_limiter": "np.maximum((c.Ca2 * c.CO3) / (0.8168 + 1e-30) - 1.0, 0.0)",
    "aragonite_diss_limiter": "np.maximum(1.0 - (c.Ca2 * c.CO3) / (0.8168 + 1e-30), 0.0)",
}

reactions = [
    # -------------------------------------------------------------------------
    # 1. Anaerobic Oxidation of Methane (AOM)
    # CH4 + SO4^2- -> HCO3^- + HS^- + H2O
    # In DIC / TA notation:
    # 1 CH4 + 1 SO4 -> 1 DIC + 2 TA + 1 TS2
    # Rate: k_AOM * [CH4] * [SO4] (already bimolecular by master/reactant)
    # -------------------------------------------------------------------------
    {
        "reaction_name": "AOM",
        "reaction": "CH4 + SO4 -> DIC + 2 TA + TS2",
        "k_value_name": "k_AOM",
    },

    # -------------------------------------------------------------------------
    # 2. Organic Carbon Degradation via Sulfate Reduction
    # 2 CH2O + SO4^2- -> 2 HCO3^- + HS^- + H+ (net: 2 DIC + 2 TA + 1 TS2)
    # Nutrients release: N: 0.151 TNH4, P: 0.0094 TPO4 per C
    # -------------------------------------------------------------------------
    {
        "reaction_name": "POC_fast_sulfate_reduction",
        "reaction": "2 POC_fast + SO4 -> 2 DIC + 2 TA + TS2 + 0.302 TNH4 + 0.0188 TPO4",
        "k_value_name": "POC_fast",
        "limiters": ["SO4"],
    },
    {
        "reaction_name": "POC_slow_sulfate_reduction",
        "reaction": "2 POC_slow + SO4 -> 2 DIC + 2 TA + TS2 + 0.302 TNH4 + 0.0188 TPO4",
        "k_value_name": "POC_slow",
        "limiters": ["SO4"],
    },

    # -------------------------------------------------------------------------
    # 3. Methanogenesis of Organic Matter (when SO4 depleted)
    # 2 CH2O -> CH4 + CO2  (1 DIC + 1 CH4)
    # -------------------------------------------------------------------------
    {
        "reaction_name": "POC_fast_methanogenesis",
        "reaction": "2 POC_fast -> CH4 + DIC + 0.302 TNH4 + 0.0188 TPO4",
        "k_value_name": "POC_fast",
        "inhibitors": ["SO4"],
    },
    {
        "reaction_name": "POC_slow_methanogenesis",
        "reaction": "2 POC_slow -> CH4 + DIC + 0.302 TNH4 + 0.0188 TPO4",
        "k_value_name": "POC_slow",
        "inhibitors": ["SO4"],
    },

    # -------------------------------------------------------------------------
    # 4. Authigenic Carbonate Mineral Dynamics
    # Precipitation: Ca2+ + CO3^2- -> CaCO3
    # Net in DIC / TA: Ca2 + DIC + 2 TA -> CaCO3
    # Dissolution: CaCO3 -> Ca2 + DIC + 2 TA
    # Note, we cannot use the coupled precipitation/dissolution scheme,
    # since bicarbonate relates to both, DIC and TA in a non-linear fashion.
    # so we need to write this as two distinct reactions 
    # -------------------------------------------------------------------------
    {
        "reaction_name": "calcite_precipitation",
        "reaction": "Ca2 + DIC + 2 TA -> Calcite",
        "k_value_name": "k_calcite_precip",
        "limiters": ["calcite_precip_limiter"],
        "reaction_order": 0,
    },
    {
        "reaction_name": "calcite_dissolution",
        "reaction": "Calcite -> Ca2 + DIC + 2 TA",
        "k_value_name": "k_calcite_diss",
        "limiters": ["calcite_diss_limiter"],
        "reaction_order": 1,
    },
    {
        "reaction_name": "aragonite_precipitation",
        "reaction": "Ca2 + DIC + 2 TA -> Aragonite",
        "k_value_name": "k_aragonite_precip",
        "limiters": ["aragonite_precip_limiter"],
        "reaction_order": 0,
    },
    {
        "reaction_name": "aragonite_dissolution",
        "reaction": "Aragonite -> Ca2 + DIC + 2 TA",
        "k_value_name": "k_aragonite_diss",
        "limiters": ["aragonite_diss_limiter"],
        "reaction_order": 1,
    },
]
