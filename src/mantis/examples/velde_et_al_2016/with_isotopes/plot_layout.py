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

Declarative layout configuration for diagenetic modeling results (Velde et al. 2016, with isotopes).
"""

# ----  color map and plot style setup  --- #
colormap = "default"
colormap_override = {
    "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True, "zorder": -10},
}

# ----  concentration data panel  --- #
liquid_species = {
    "scale": "log",
    "species": [
        "O2",
        "SO4",
        "TS2",
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
        "TS2",
        "Fe3",
        "FeS",
        "D_bio",
    ],
}

# ----  isotope data panel  --- #
isotopes = {
    "species": ["d34S_SO4", "d34S_TS2", "d34S_FeS"],
    "lim": (-45, -37, 20, 30),  # broken y-axis: exclude -37 to 20
}


plot_panels = [
    [liquid_species, solid_species],
    [reaction_rates],
    [isotopes],
]
