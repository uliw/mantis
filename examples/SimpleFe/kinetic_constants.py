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

from typing import Any


def get_reaction_constants(pH: float, phi: float | Any) -> dict:
    """Define reaction constants in initial units; diff_lib converts to model SI units."""
    constants: dict = {
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

        # Dynamic proton concentration in mol/m^3
        "Hplus": f"{1000.0 * 10.0 ** (-float(pH))} mol/m^3",

        # Partitioning coefficient for Fe2+ adsorption onto POC (Du 2023, Tables 4 & 6)
        "Fe2_p_eq": "0.8444 dimensionless",
    }

    return constants


if __name__ == "__main__":
    from mantis.diff_lib import convert_reaction_constants

    constants = get_reaction_constants(7.5, 0.8)
    print("Defined constants:", constants)
