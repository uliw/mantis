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
"""Define reaction constants.

Note, Velde et al likely have a typo in the FeS equilibrium constant.
Rickard (2006) defines two reaction constants, one with Ksp = 3.5
and the other with Ksp = - 3.5. Here we use the latter since otherwise
we would never precipitate FeS unless H2S > 1 mmol/l
"""

kinetic_constants: dict = {
    "POC": "0.001 1/year",
    "K_SO4": "0.9 umol/cm^3",
}


if __name__ == "__main__":
    from mantis.diff_lib import convert_reaction_constants

    k_values = convert_reaction_constants(kinetic_constants)
    for k, v in k_values.items():
        print(f"{k} = {v:.2e}")
