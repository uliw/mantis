"""
plot_layout.py

Declarative layout configuration for diagenetic modeling results (Velde et al. 2016, no isotopes).
"""

# ----  optional comparison against existing data --- #
comparison_data = "velde_bt04_d6.csv"

# Species mapping between mantis and external data.
# The second list element is a conversion factor or function.
species_mapping = {
    "z": ["depth_cm", 1 / 100],
    "HS": ["C_HS", 1.0],
    "FeS": ["C_FeS", "wt"],
    "Fe3": ["C_FeOOH", "wt"],
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
