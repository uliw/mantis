"""
plot_layout.py

Declarative layout configuration for diagenetic modeling results (Carbonate Test Fixture).
"""

colormap = "default"
colormap_override = {
    "SO4": {"color": "cornflowerblue", "lw": 2, "fill_only": True, "zorder": -10},
    "FeS": {"color": "#ee3377", "lw": 4},
    "FeS2": {"color": "darkorange", "lw": 2, "fill_only": True, "zorder": -9},
    "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True},
}

liquid_species = {
    "scale": "log",
    "species": ["SO4", "O2", "TS2", "Fe2_liq"],
    "lim": (1e-6, 50),
    "units": "mmol/l",
}

solid_species = [
    {
        "scale": "lin",
        "units": ["wt"],
        "relative_to": "Fe",
        "species": ["FeS2", "FeS", "Fe3"],
        "lim": (0, 1),
    },
    {
        "scale": "lin",
        "units": ["wt"],
        "relative_to": "S",
        "species": ["S0"],
        "lim": (0, 0.5),
    },
]

reaction_rates = {
    "species": ["SO4", "FeS", "Fe3", "D_bio", "TS2", "DIC", "TA", "CaCO3"],
}

isotopes = {
    "species": ["d_SO4", "d_TS2", "d_FeS", "d_FeS2", "d_S0"],
    "lim": (-45, -37, 20, 30),
}
