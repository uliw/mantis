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
"""kinetic_constants.py

Thermodynamic & Kinetic Constants for the SimpleFe Benchmark Model
(Du 2023, Geosci. Model Dev., 16, 5865-5894).
==================================================================

Overview & Unit Conversion Reference:
-------------------------------------
All constants declared here are translated from Du (2023) Table 6 into Mantis's internal
SI unit system (meters, seconds, moles):
  • Time: 1 year = 365.25 * 86400 seconds = 3.15576 * 10^7 s.
  • Concentration: 1 mmol/cm^3 = 1000 mol/m^3 = 1000 mM = 10^6 uM.
  • Degradation rate k_POC = 0.01 yr^-1 = 3.16888 * 10^-10 s^-1.
  • Half-saturation constant K_FeOOH = 3.0 mmol/cm^3_solid = 3000.0 mol/m^3_solid.
  • Half-saturation constant K_SO4 = 1.0 * 10^-3 mmol/cm^3_pw = 1.0 mol/m^3_pw.
  • Sulfide oxidation rate k_FeOOHH2S = 4000 (mmol/cm^3)^-1 yr^-1 = 4.0 m^3/(mol * yr).
  • FeS solubility product Ksp_FeS = 10^(-3.2) (mmol/cm^3)^-1 -> 1000 * 10^(-3.2) = 0.630957 mol/m^3.
  • FeS precipitation rate kFeSpre = 0.2 * ds_rho = 0.52 mmol/(cm^3 * yr) = 520 mol/(m^3 * yr).
"""

kinetic_constants: dict = {
    # Organic carbon remineralization rate constant (SimpleFe benchmark: constant 0.01 yr^-1)
    # Note: If age-dependent scaling is desired in custom models, use:
    #   from mantis.diff_lib import sed_age
    #   age = sed_age(mode="relative", a0=1.0, a=0.21, b=1.22)
    #   "k_POC": "age * 0.01 1/year"
    "k_POC": "0.01 1/year",

    # Half-saturation constants
    "K_FeOOH": "3000.0 mol/m^3",
    "K_SO4": "1.0 mol/m^3",

    # Sulfide oxidation by iron oxide (Du 2023: 4000 (mmol/cm^3)^-1 yr^-1 = 4.0 (mol/m^3)^-1 yr^-1)
    "k_FeOOHH2S": "4.0 m^3/(mol*year)",

    # Iron monosulfide (FeS) precipitation & dissolution
    # KspFeS = 10^(-3.2) (mmol/cm^3)^-1 in Du's scale -> 0.630957 mol/m^3
    "FeS_sp": f"{1000.0 * 10.0 ** (-3.2)} mol/m^3",
    # kFeSpre = 0.2 * ds_rho = 0.52 mmol/cm^3/yr = 520 mol/(m^3 * yr)
    "FeS_isp": "520.0 mol/(m^3*year)",
    # Reversible dissolution rate constant
    "FeS_isd": "0.01 1/year",

    # Dynamic proton concentration in mol/m^3 (pH = 7.5)
    "Hplus": f"{1000.0 * 10.0 ** (-7.5)} mol/m^3",

    # Partitioning coefficient for Fe2+ adsorption onto POC (Du 2023, Tables 4 & 6)
    "Fe2_p_eq": "0.8444 dimensionless",
}
