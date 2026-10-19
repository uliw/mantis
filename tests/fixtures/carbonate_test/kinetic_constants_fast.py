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

    if hasattr(phi, "value"):
        phi = phi.value

    velde: dict = {
        "POC_fast": [10, "1/year", "1/second"],  # highly reactive OM
        "POC_slow": [0.1, "1/year", "1/second"],  # reactive OM, Velde et al
        "Fe2_p_eq": [
            696,
            "dimensionless",
            "dimensionless",
        ],  # sorbed vs Fe2+ Fe2+ liq.
        "FeS2_O2": [
            1e-10,
            "m^3/(mol*second)",
            "m^3/(mol*second)",
        ],  # FeS2 + O2 -> SO4, Halevy et al
        "FeS_S0": [  # FeS + S0 -> FeS2
            5e-8,
            "m^3/(mol*second)",
            "m^3/(mol*second)",
        ],  # FeS + HS -> FeS2
        "FeS_TS2": [
            5e-8,
            "m^3/(mol*second)",
            "m^3/(mol*second)",
        ],  # FeS + H2S -> FeS2, at 10C -> notes.org
        "Hplus": [10 ** (-pH), "mol/l", "mol/m^3"],
        # Oxidation reactions after Velde at eal 2016
        "TS2_O2": [1e7, "cm^3/(umol*year)", "m^3/(mol*second)"],
        "Fe3_hs": [
            494,
            "cm^3/(umol*year)",
            "m^3/(mol*second)",
        ],  # check with halevy
        "FeS_O2": [1e7, "cm^3/(umol*year)", "m^3/(mol*second)"],
        "Fe2_O2": [1e7, "cm^3/(umol*year)", "m^3/(mol*second)"],
        # "Fe2p_ox": [1e73, "cm^3/(umol*year)", "m^3/(mol*second)"],
        "S0_O2": [4e7, "cm^3/(umol*year)", "m^3/(mol*second)"],
        # Iron sulfide reactions after Velde at eal 2016
        "FeS_isp": [
            1e4,
            "umol/((cm^3*year))",
            "mol/(m^3 *second)",
        ],  # FeS precipitation
        "FeS_isd": [3, "1/year", "1/second"],  # FeS dissolution
        "FeS_sp": [  # FeS equilibrium constant
            10**-3.5,
            "mol/l",
            "mol/m^3 ",
        ],
        # S0 disproportionation rate (first-order in S0)
        "S0_dispro": [50, "1/year", "1/second"],
        # Monod half-saturation & inhibition thresholds
        "K_O2": [0.001, "umol/cm^3", "mol/m^3"],
        "K_O2_TS2": [0.001, "umol/cm^3", "mol/m^3"],
        "K_SO4": [0.9, "umol/cm^3", "mol/m^3"],
        "K_Fe3_diss_red": [10.4, "umol/cm^3", "mol/m^3"],
        "K_Fe3": [1e-3, "dimensionless", "dimensionless"],
        # Carbonate system constants
        # Calcite dissolution and precipitation rate constants:
        # Jourabchi, Van Cappellen, & Regnier (2005), Am. J. Sci. 305: 919-956, Table 2.
        # Wang & Van Cappellen (1996), Geochim. Cosmochim. Acta 60: 2993-3014.
        "pK1_carb": [5.86, "dimensionless", "dimensionless"],
        "pK2_carb": [8.92, "dimensionless", "dimensionless"],
        "Ksp_calcite": [10**-6.37, "(mol/l)^2", "(mol/m^3)^2"],
        "k_calcite_diss": [3.65, "1/year", "1/second"],  # 0.01 d^-1 = 3.65 yr^-1 (Jourabchi et al. 2005)
        "k_calcite_precip": [1.0, "mol/(m^3*year)", "mol/(m^3*second)"],  # 1.0 mol/(m^3*yr) (Jourabchi et al. 2005)
    }

    for k, v in velde.items():
        k_values[k] = Q_(f"{v[0]} {v[1]}").to(f"{v[2]}").magnitude

    return velde, k_values


if __name__ == "__main__":
    import pint

    ureg = pint.UnitRegistry()
    Q_ = ureg.Quantity

    velde, k_values = get_reaction_constants(7.5, 0.8)
    ureg = pint.UnitRegistry()
    for k, v in velde.items():
        n = Q_(f"{v[0]} {v[1]}").to(f"{v[2]}")
        print(f"{k} = {n.magnitude:.2e} [{n.units:~P}]")

    print(f"k.FeS_sp * k.Hplus: {k_values['FeS_sp'] * k_values['Hplus']:.2e}")
    print(f"k.Hplus: {k_values['Hplus']:.2e}")
