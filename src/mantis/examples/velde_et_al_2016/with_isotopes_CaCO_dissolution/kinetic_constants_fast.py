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
    "POC_fast": "10 1/year",  # highly reactive OM
    "POC_slow": "0.1 1/year",  # reactive OM, Velde et al
    "Fe2_p_eq": "696 dimensionless",  # sorbed vs Fe2+ Fe2+ liq.
    # FeS + HS -> FeS2
    "FeS_TS2": "5e-8 m^3/(mol*second)",  # FeS + H2S -> FeS2, at 10C -> notes.org
    "Hplus": f"{10 ** (-7.5)} mol/l",
    # Oxidation reactions after Velde at eal 2016
    "TS2_O2": "1e7 cm^3/(umol*year)",
    "Fe3_hs": "494 cm^3/(umol*year)",  # check with halevy
    "FeS_O2": "1e7 cm^3/(umol*year)",
    "Fe2_O2": "1e7 cm^3/(umol*year)",
    # Iron sulfide reactions after Velde at eal 2016
    "FeS_isp": "1e4 umol/((cm^3*year))",  # FeS precipitation
    "FeS_isd": "3 1/year",  # FeS dissolution
    "FeS_sp": "10**-3.5 mol/l",  # FeS equilibrium constant
    # Monod half-saturation & inhibition thresholds
    "K_O2": "0.001 umol/cm^3",
    "K_O2_TS2": "0.001 umol/cm^3",
    "K_SO4": "0.9 umol/cm^3",
    "K_Fe3_diss_red": "10.4 umol/cm^3",
    "K_Fe3": "1e-3 dimensionless",
    "k_calcite_precip": "0.08 mol/(m^3*year)",
    "k_calcite_diss": "0.015, 1/year",
}


if __name__ == "__main__":
    from mantis.diff_lib import convert_reaction_constants

    k_values = convert_reaction_constants(kinetic_constants)
    for k, v in k_values.items():
        print(f"{k} = {v:.2e}")

    print(f"k.FeS_sp * k.Hplus: {k_values['FeS_sp'] * k_values['Hplus']:.2e}")
    print(f"k.Hplus: {k_values['Hplus']:.2e}")
