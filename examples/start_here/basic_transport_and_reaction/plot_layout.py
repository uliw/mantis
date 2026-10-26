"""
plot_layout.py

Declarative layout configuration for diagenetic modeling results (Velde et al. 2016, no isotopes).
"""

# ----  optional comparison against existing data --- #
# comparison_data = ""

# Species mapping between mantis and external data.
# The second list element is a conversion factor or function.
# species_mapping = {
#     "z": ["depth_cm", 1 / 100],
#     "HS": ["C_HS", 1.0],
#     "FeS": ["C_FeS", "wt"],
#     "Fe3": ["C_FeOOH", "wt"],
# }

# ----  color map and plot style setup  --- #
colormap = "default"
# colormap_override = {
#     "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True},
# }

#  ---- z-axis (Depth) Scale and Limits ----
z_scale = "lin"  # Or "linear" / "lin" (default: "log")
z_lim = (0, 1)  # (z_min, z_max) in meters (default for log: (1e-4, display_length))

# ----  concentration data panel  --- #
liquid_species_1 = liquid_species({
    "scale": "lin",
    "species": ["SO4", "TS2"],
    "lim": (0, 30),
    "units": "mmol/l",
})

solid_species_1 = solid_species([
    {
        "scale": "lin",
        "units": ["wt"],
        "relative_to": "C",
        "species": ["POC"],
        # "lim": (0, 1),
    }
])

# ----  reaction rate data panel  --- #
reaction_rates_1 = reaction_rates({"species": ["POC"]})

# ---- create final panel arrangement ------ #
plot_panels = [
    [liquid_species_1, solid_species_1],
    [reaction_rates_1],
]
