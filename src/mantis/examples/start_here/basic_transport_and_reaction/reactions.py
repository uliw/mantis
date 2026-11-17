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


# Monod type reaction rate limiters 
limiters = {
    # Base implicit limiters: L_imp = 1.0 / (C + K)
    # The parser automatically derives explicit (C * L_imp) and inhibition (K * L_imp = 1 - C * L_imp)
    "SO4": "1.0 / (c.SO4 + K_SO4)",
}

reactions = [
    # -------------------------------------------------------------------------
    # Organic Carbon Primary Respiration
    # -------------------------------------------------------------------------
    {
        "reaction_name": "sulfate_reduction",
        "reaction": "2 POC + SO4 -> TS2",
        "k_value_name": "POC",
        "limiters": ["SO4"],
    },
]
   
