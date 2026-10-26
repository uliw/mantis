"""Define reaction constants.

Note, Velde et al likely have a typo in the FeS equilibrium constant.
Rickard (2006) defines two reaction constants, one with Ksp = 3.5
and the other with Ksp = - 3.5. Here we use the latter since otherwise
we would never precipitate FeS unless H2S > 1 mmol/l
"""

from typing import Any


def get_reaction_constants(
    pH: float, phi: float | Any, k_values: dict = {}
) -> dict:
    """Convert reaction_constants into model units.

        It is assumed that the reaction constants are in phase specific units.
        This can be tricky, as Velde et al, state that their k-values are in
        bulk units, but then they convert it before using (see Table 3). E.g.
        they write:

       (1−ϕ)⋅kSIR⋅[FeOOH]⋅[HS−]
    ​
        which suggests that the (1−ϕ) is used to convert FeOOH from solid
        phase concentration to bulk concentration. Since this code
        tracks all concentrations in phase specific units, converting the
        k-value by deviding by (1−ϕ), would cancel out, and we can use the
        k-value as published.
    """
    import pint

    ureg = pint.UnitRegistry()
    Q_ = ureg.Quantity

    constants: dict = {
        # key, [value, unit, target unit]
        "POC": [0.001, "1/year", "1/second"],  #
        # Monod half-saturation & inhibition thresholds
        "K_SO4": [0.9, "umol/cm^3", "mol/m^3"],
    }

    for k, v in constants.items():
        k_values[k] = Q_(f"{v[0]} {v[1]}").to(f"{v[2]}").magnitude

    return constants, k_values

