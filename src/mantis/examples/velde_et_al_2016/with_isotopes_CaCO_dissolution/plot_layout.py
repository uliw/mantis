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
from mantis import liquid_species, solid_species, isotopes, reaction_rates

# ----  color map and plot style setup  --- #
colormap = "default"
colormap_override = {
    "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True, "zorder": -10},
}

# ----  concentration data panel  --- #
liquid_species_1 = liquid_species({
    # "scale": "log",
    "species": [
        #"O2",
        # "SO4",
        # "TS2",
        # "Fe2_liq",
        "Ca2",
        "TA",
        "DIC",
    ],
    #"lim": (1e-6, 50),
    "units": "mmol/l",
})

solid_species_1 = solid_species([
    {
        "scale": "lin",
        "units": "wt%",
        "relative_to": "CaCO3",
        "species": ["CaCO3"],
        # "lim": (0, 5),
    },
])

pH = solid_species([
    {
        "scale": "lin",
        "units": "pH",
        "species": ["pH"],
        #"lim": (6.5, 9),
    },
    {
        "scale": "lin",
        "units": r"$\Omega$",
        "species": ["Omega_calcite"],
        #"lim": (6.5, 9),
    }
])

# ----  reaction rate data panel  --- #
reaction_rates_1 = reaction_rates(
{
    "species": [
        "CaCO3",
        "D_bio",
    ],
})

# ----  isotope data panel  --- #
isotopes_1 = isotopes({
    "species": ["d34S_SO4", "d34S_TS2", "d34S_FeS"],
    "lim": (-45, -37, 20, 30),  # broken y-axis: exclude -37 to 20
})

isotopes_2 = isotopes({
    "species": ["d13C_DIC", "d13C_CaCO3"],
    # "lim": (-45, -37, 20, 30),  # broken y-axis: exclude -37 to 20
})


plot_panels = [
    [liquid_species_1, solid_species_1],
    [reaction_rates_1, pH],
    [isotopes_1, isotopes_2],
]
