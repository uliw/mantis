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
"""Define reaction constants."""

kinetic_constants: dict = {
    "POC_fast": "0.2 1/year",
    "POC_slow": "3.0e-4 1/year",
    "k_AOM": "2500.0 cm^3/(mmol*year)",
    "K_SO4": "0.5 umol/cm^3",
    "k_calcite_precip": "0.08 mol/(m^3*year)",
    "k_calcite_diss": "0.015 1/year",
    "k_aragonite_precip": "0.8 mol/(m^3*year)",
    "k_aragonite_diss": "5.0 1/year",
}
