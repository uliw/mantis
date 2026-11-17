"""
plot_layout.py

Declarative layout configuration for diagenetic modeling results (Velde et al. 2016, no isotopes).
"""

# ----  optional comparison against existing data --- #
comparison_data = "14_with_irr.csv"

# Species mapping between mantis and external data.
# Format: "mantis_species": [external_column_name, source_unit, target_display_unit]
species_mapping = {
    "z": ["depth_cm", "cm", "m"],
    "HS": ["C_HS", "mmol/L_pw", "mmol/L"],
    "FeS": ["C_FeS", "mmol/L_solid", "wt%"],
    "Fe3": ["C_FeOOH", "mmol/L_solid", "wt%"],
    "Fe2_liq": ["C_Fe", "mmol/L_pw", "mmol/L"],
    "O2": ["C_O2", "mmol/L_pw", "mmol/L"],
}

# ----  color map and plot style setup  --- #
colormap = "default"
colormap_override = {
    "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True},
}

# ----  concentration data panel  --- #
liquid_species = {
    "scale": "log",
    "species": [
        "O2",
        "SO4",
        "HS",
        "Fe2_liq",
    ],
    "lim": (1e-6, 50),
    "units": "mmol/l",
}

solid_species = [
    {
        "scale": "lin",
        "units": ["wt"],
        "relative_to": "Fe",
        "species": ["FeS", "Fe3"],
        "lim": (0, 1),
    }
]

# ----  reaction rate data panel  --- #
reaction_rates = {
    "species": [
        "POC",
        "SO4",
        "HS",
        "Fe3",
        "FeS",
        "D_bio",
    ],
}

plot_panels = [
    [liquid_species, solid_species],
    [reaction_rates],
]

