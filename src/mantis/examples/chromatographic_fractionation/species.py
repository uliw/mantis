"""species.py

Declarative Chemical Species Definitions for Chromatographic Isotope Fractionation.
===================================================================================

This benchmark simulates the 1D advective-dispersive transport of dissolved ferrous iron
invading an initial column devoid of iron, undergoing linear equilibrium sorption to sediment:
    Fe2_liq <-> Fe2_sorbed

With isotopic equilibrium fractionation between dissolved and sorbed iron (d56Fe standard IRMM-014):
    Fe2_sorbed_54 <-> Fe2_liq_54   (alpha_eq = 1.002)

Because the heavy isotope (56Fe) sorbs preferentially (alpha = 1.002), the light companion
(54Fe) is enriched in the dissolved phase, giving it a higher dissolved fraction f_dis,54 > f_dis.
Consequently, 54Fe experiences less retardation and migrates faster than bulk 56Fe,
producing chromatographic separation and a pronounced negative delta56Fe excursion
at the invasion front.
"""

# -----------------------------------------------------------------------------
# Isotope Reference Standards
# -----------------------------------------------------------------------------
isotope_standards: dict = {
    "d56Fe": {
        "element": "Fe",
        "standard": "IRMM-014",
        "ratio": 15.68376068,  # 56Fe / 54Fe standard atomic ratio (Mantis convention: heavy/light)
        "unit": "mUr",
    },
}

species = {
    # -------------------------------------------------------------------------
    # Mixed-phase Fe2+ Pool (Porewater + Sorbed Fraction)
    # -------------------------------------------------------------------------
    "Fe2_total": {
        "phase": "mixed",
        "basis": "porewater",
        "include": True,
        "formula": "Fe",
        "diffusion": {"m0": 3.31, "m1": 0.150},
        "isotopes": {"Fe2_total_54": "Fe"},
        "bc_top": {"Dirichlet": 1.0, "d56Fe": 0.0},
        "prefill": False,
    },
}

# =============================================================================
# Speciation & Phase Partitioning
# =============================================================================
speciation = {
    "fe2_partitioning": {
        "type": "partition",
        "primary": "Fe2_total",
        "fractions": {
            "dissolved": "0.5",
            "sorbed": "0.5",
        },
        "equilibrium_fractionation": {
            "standard": "d56Fe",
            "alpha": {"Fe2_sorbed_54_to_Fe2_liq_54": 0.998},
        },
        "species": {
            "Fe2_liq": {
                "phase": "liquid",
                "alias_of": "Fe2_total.dissolved",
                "isotopes": {"Fe2_liq_54": "Fe"},
            },
            "Fe2_sorbed": {
                "phase": "solid",
                "alias_of": "Fe2_total.sorbed",
                "isotopes": {"Fe2_sorbed_54": "Fe"},
            },
        },
    },
}
