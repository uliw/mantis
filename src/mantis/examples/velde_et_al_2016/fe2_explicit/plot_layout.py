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
comparison_data = "velde_bt04_d6.csv"

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

