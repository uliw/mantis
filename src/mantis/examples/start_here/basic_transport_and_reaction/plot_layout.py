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
    "species": ["SO4", "H2S_tot"],
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
