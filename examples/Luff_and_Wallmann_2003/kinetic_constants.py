from typing import Any
import pint

def get_reaction_constants(pH: float, phi: float | Any, k_values: dict = {}) -> tuple[dict, dict]:
    ureg = pint.UnitRegistry()
    Q_ = ureg.Quantity
    constants_def: dict = {
        "POC_fast": [0.2, "1/year", "1/second"],
        "POC_slow": [3.0e-4, "1/year", "1/second"],
        "k_AOM": [2500.0, "cm^3/(mmol*year)", "m^3/(mol*second)"],
        "K_SO4": [0.5, "umol/cm^3", "mol/m^3"],
        "k_calcite_precip": [0.08, "mol/(m^3*year)", "mol/(m^3*second)"],
        "k_calcite_diss": [0.015, "1/year", "1/second"],
        "k_aragonite_precip": [0.8, "mol/(m^3*year)", "mol/(m^3*second)"],
        "k_aragonite_diss": [5.0, "1/year", "1/second"],
    }
    for k, v in constants_def.items():
        k_values[k] = Q_(f"{v[0]} {v[1]}").to(f"{v[2]}").magnitude
    return constants_def, k_values
