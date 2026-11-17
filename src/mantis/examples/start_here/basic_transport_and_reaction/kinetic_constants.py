"""Define reaction constants.

Note, Velde et al likely have a typo in the FeS equilibrium constant.
Rickard (2006) defines two reaction constants, one with Ksp = 3.5
and the other with Ksp = - 3.5. Here we use the latter since otherwise
we would never precipitate FeS unless H2S > 1 mmol/l
"""

from typing import Any


def get_reaction_constants(
    pH: float, phi: float | Any
) -> dict:
    """Return reaction constants with their initial units.

    The model operates in SI units (meter, second, mole), and unit conversion
    is automatically handled by diff_lib.py.
    """
    constants: dict = {
        "POC": "0.001 1/year",
        "K_SO4": "0.9 umol/cm^3",
    }

    return constants


if __name__ == "__main__":
    from mantis.diff_lib import convert_reaction_constants

    constants = get_reaction_constants(7.5, 0.8)
    k_values = convert_reaction_constants(constants)
    for k, v in k_values.items():
        print(f"{k} = {v:.2e}")
