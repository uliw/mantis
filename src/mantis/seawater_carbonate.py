"""seawater_carbonate.py

Thermodynamic equilibrium calculations for the marine carbonate system as a function
of temperature, salinity, and hydrostatic pressure.

Formulations and Citations:
---------------------------
1. Carbonic Acid Dissociation Constants (K1, K2):
   - Lueker, T. J., Dickson, A. G., and Keeling, C. D. (2000). Ocean pCO2 calculated from
     dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based
     on laboratory measurements of CO2 in gas and seawater at equilibrium.
     Marine Chemistry, 70(1-3), 105-119.
     (Refit of Mehrbach et al. 1973 onto the Total pH scale).

2. Calcite Stoichiometric Solubility Product (K'sp_calcite):
   - Mucci, A. (1983). The solubility of calcite and aragonite in seawater at various
     salinities and temperatures, and at one atmosphere total pressure.
     American Journal of Science, 283(7), 780-799.

3. Hydrostatic Pressure Corrections (K1, K2, K'sp):
   - Millero, F. J. (1995). Thermodynamics of the carbon dioxide system in the oceans.
     Geochimica et Cosmochimica Acta, 59(4), 661-677.
   - Ingle, S. E. (1975). Solubility of calcite in the ocean.
     Marine Chemistry, 3(4), 301-319.
   - Millero, F. J. (1979). The thermodynamics of the carbonate system in seawater.
     Geochimica et Cosmochimica Acta, 43(10), 1651-1661.

4. Seawater Density (rho_sw):
   - Millero, F. J., and Poisson, A. (1981). International one-atmosphere equation of
     state of seawater. Deep Sea Research Part A, 28(6), 625-629.
     (UNESCO 1981 equation of state).
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple
import numpy as np


# Universal gas constant in cm^3 * bar / (mol * K)
# R = 8.314462618 J / (mol * K) = 83.14462618 cm^3 * bar / (mol * K)
R_GAS_BAR = 83.1446261815324


@dataclass
class CarbonateConstants:
    """Stores thermodynamic carbonate equilibrium constants and seawater density."""
    temp_c: float
    temp_k: float
    salinity: float
    pressure_bar: float
    rho_sw: float        # kg/m^3
    K1: float            # mol/m^3 (= mmol/L)
    K2: float            # mol/m^3 (= mmol/L)
    pK1: float
    pK2: float
    Ksp_calcite: float   # (mol/m^3)^2 (= (mmol/L)^2)
    K1_mol_kg: float     # mol/kg-SW
    K2_mol_kg: float     # mol/kg-SW
    Ksp_mol_kg: float    # (mol/kg-SW)^2
    ph_scale: str = "free"
    scale_factor: float = 1.0  # (1 + St / Ks) conversion factor between Total and Free scales


def calc_total_sulfate(salinity: float) -> float:
    """Calculates total sulfate concentration in mol/kg-SW from practical salinity.

    References:
    -----------
    - Morris, A. W., and Riley, J. P. (1966). The bromide/chlorinity and
      sulphate/chlorinity ratio in sea water. Deep Sea Research 13, 699-705.
    - Dickson, A. G., Sabine, C. L., and Christian, J. R. (2007). Guide to best
      practices for ocean CO2 measurements. PICES Special Publication 3, Chap. 5.
    """
    S = float(salinity)
    return (0.1400 / 96.062) * (S / 1.80655)


def calc_bisulfate_constant(
    temp_c: float,
    salinity: float,
    pressure_bar: float = 1.0,
) -> float:
    """Calculates stoichiometric dissociation constant of bisulfate (HSO4-) in seawater on the free scale.

    References:
    -----------
    - 1 atm constant (KS0):
      Dickson, A. G. (1990). Standard potential of the reaction:
      AgCl(s) + 1/2 H2(g) = Ag(s) + HCl(aq), and the standard acidity constant of the ion
      HSO4- in synthetic seawater from 273.15 to 318.15 K.
      The Journal of Chemical Thermodynamics, 22(2), 113-127.
    - Pressure correction:
      Millero, F. J. (1995). Thermodynamics of the carbon dioxide system in the oceans.
      Geochimica et Cosmochimica Acta, 59(4), 661-677 (Table 9).

    Parameters:
    -----------
    temp_c : float
        Temperature in degrees Celsius.
    salinity : float
        Practical salinity in PSU.
    pressure_bar : float, optional
        Hydrostatic pressure in bars (default: 1.0).

    Returns:
    --------
    float
        Stoichiometric KS in mol/kg-SW on the free scale.
    """
    T_c = float(temp_c)
    T_k = T_c + 273.15
    S = float(salinity)
    P_bar = float(pressure_bar)

    # Ionic strength on free scale (Dickson 1990)
    I = 19.924 * S / (1000.0 - 1.005 * S)

    ln_Ks_0 = (
        -4276.1 / T_k
        + 141.328
        - 23.093 * np.log(T_k)
        + (-13856.0 / T_k + 324.57 - 47.986 * np.log(T_k)) * np.sqrt(I)
        + (35474.0 / T_k - 771.54 + 114.723 * np.log(T_k)) * I
        - (2698.0 / T_k) * (I**1.5)
        + (1776.0 / T_k) * (I**2)
        + np.log(1.0 - 0.001005 * S)
    )
    Ks_0 = np.exp(ln_Ks_0)

    # Pressure correction (Millero 1995)
    deltaVS = -18.03 + 0.0466 * T_c + 0.000316 * (T_c**2)
    kappaS = (-4.53 + 0.09 * T_c) / 1000.0
    lnKsfac = (-deltaVS + 0.5 * kappaS * P_bar) * P_bar / (R_GAS_BAR * T_k)
    Ks_P = Ks_0 * np.exp(lnKsfac)

    return float(Ks_P)


def calc_seawater_density(temp_c: float, salinity: float = 34.0) -> float:
    """Calculates seawater density at 1 atm pressure using the UNESCO 1981 equation of state.

    Reference:
    ----------
    Millero, F. J., and Poisson, A. (1981). International one-atmosphere equation
    of state of seawater. Deep Sea Research Part A, 28(6), 625-629.

    Parameters:
    -----------
    temp_c : float
        Temperature in degrees Celsius.
    salinity : float, optional
        Practical salinity in PSU (default: 34.0).

    Returns:
    --------
    float
        Density in kg/m^3 (divide by 1000 for kg/L or g/cm^3).
    """
    T = float(temp_c)
    S = float(salinity)

    # Pure water density SMOW (kg/m^3)
    rho_w = (
        999.842594
        + 6.793952e-2 * T
        - 9.095290e-3 * T**2
        + 1.001685e-4 * T**3
        - 1.120083e-6 * T**4
        + 6.536332e-9 * T**5
    )

    # Salinity contribution terms
    A = (
        0.824493
        - 4.0899e-3 * T
        + 7.6438e-5 * T**2
        - 8.2467e-7 * T**3
        + 5.3875e-9 * T**4
    )
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
    C = 4.8314e-4

    rho = rho_w + A * S + B * (S**1.5) + C * (S**2)
    return float(rho)


def calc_carbonate_constants(
    temp_c: float,
    salinity: float = 34.0,
    pressure_bar: float = 1.0,
    ph_scale: str = "free",
) -> CarbonateConstants:
    """Calculates K1, K2, and K'sp(calcite) corrected for temperature, salinity, pressure, and pH scale.

    References:
    -----------
    - K1, K2 (Total scale): Lueker, Dickson, & Keeling (2000), Marine Chemistry 70: 105-119.
    - Ksp: Mucci (1983), American Journal of Science 283: 780-799.
    - Pressure corrections: Millero (1995), GCA 59: 661-677; Ingle (1975), Marine Chemistry 3: 301-319.
    - Density: Millero & Poisson (1981), Deep Sea Research 28A: 625-629.
    - Free vs Total pH scale conversion: Dickson (1990), J. Chem. Thermodyn. 22: 113-127;
      Dickson et al. (2007), Guide to best practices for ocean CO2 measurements, PICES Spec. Publ. 3.

    Parameters:
    -----------
    temp_c : float
        Temperature in degrees Celsius.
    salinity : float, optional
        Salinity in PSU (default: 34.0).
    pressure_bar : float, optional
        Hydrostatic pressure in bars (default: 1.0 bar = ~10 m depth / sea surface).
    ph_scale : str, optional
        pH scale convention: 'free' (default) or 'total'.
        - 'free': [H+] considers only free hydronium ions ([H+]_F).
        - 'total': [H+] includes bisulfate ion pairs ([H+]_T = [H+]_F * (1 + St / Ks)).

    Returns:
    --------
    CarbonateConstants
        Dataclass containing constants on both mol/m^3 (mmol/L) and mol/kg-SW scales.
    """
    T_c = float(temp_c)
    T_k = T_c + 273.15
    S = float(salinity)
    P_bar = float(pressure_bar)
    scale_norm = str(ph_scale).lower().strip()
    if scale_norm not in ("free", "total"):
        raise ValueError(f"Invalid ph_scale '{ph_scale}'. Supported scales are 'free' and 'total'.")

    # 1. First and Second Dissociation Constants of Carbonic Acid at 1 atm (LDK00)
    # Total pH scale, in mol/kg-SW
    pK1_0 = (
        3633.86 / T_k
        - 61.2172
        + 9.6777 * np.log(T_k)
        - 0.011555 * S
        + 0.0001152 * S**2
    )
    K1_0 = 10.0 ** (-pK1_0)

    pK2_0 = (
        471.78 / T_k
        + 25.929
        - 3.16967 * np.log(T_k)
        - 0.01781 * S
        + 0.0001122 * S**2
    )
    K2_0 = 10.0 ** (-pK2_0)

    # Pressure correction for K1 and K2 (Millero 1995)
    deltaV1 = -25.5 + 0.1271 * T_c
    Kappa1 = (-3.08 + 0.0877 * T_c) / 1000.0
    lnK1fac = (-deltaV1 + 0.5 * Kappa1 * P_bar) * P_bar / (R_GAS_BAR * T_k)
    K1_mol_kg_total = K1_0 * np.exp(lnK1fac)

    deltaV2 = -15.82 - 0.0219 * T_c
    Kappa2 = (1.13 - 0.1475 * T_c) / 1000.0
    lnK2fac = (-deltaV2 + 0.5 * Kappa2 * P_bar) * P_bar / (R_GAS_BAR * T_k)
    K2_mol_kg_total = K2_0 * np.exp(lnK2fac)

    # 2. Conversion between Total scale and Free scale:
    # [H+]_T = [H+]_F * (1 + St / Ks)
    # K_free = K_total / (1 + St / Ks)
    Ks_P = calc_bisulfate_constant(temp_c=T_c, salinity=S, pressure_bar=P_bar)
    St = calc_total_sulfate(salinity=S)
    scale_factor = 1.0 + (St / Ks_P)

    if scale_norm == "free":
        K1_mol_kg = K1_mol_kg_total / scale_factor
        K2_mol_kg = K2_mol_kg_total / scale_factor
    else:
        K1_mol_kg = K1_mol_kg_total
        K2_mol_kg = K2_mol_kg_total

    # 3. Calcite Stoichiometric Solubility Product at 1 atm (Mucci 1983)
    # in (mol/kg-SW)^2  (invariant with respect to pH scale)
    logKsp_0 = (
        -171.9065
        - 0.077993 * T_k
        + 2839.319 / T_k
        + 71.595 * np.log10(T_k)
        + (-0.77712 + 0.0028426 * T_k + 178.34 / T_k) * np.sqrt(S)
        - 0.07711 * S
        + 0.0041249 * (S**1.5)
    )
    Ksp_0 = 10.0**logKsp_0

    # Pressure correction for Calcite Ksp (Ingle 1975, Millero 1979/1995)
    deltaVKCa = -48.76 + 0.5304 * T_c
    KappaKCa = (-11.76 + 0.3692 * T_c) / 1000.0
    lnKspfac = (-deltaVKCa + 0.5 * KappaKCa * P_bar) * P_bar / (R_GAS_BAR * T_k)
    Ksp_mol_kg = Ksp_0 * np.exp(lnKspfac)

    # 4. Seawater Density and Conversion to MANTIS Units
    # For H+ speciation in [H+] = 10^-pH, K1 and K2 are in mol/kg-SW.
    # For Ksp = [Ca2+][CO32-], concentrations are in mol/m^3 (= mmol/L), so:
    # Ksp (mol/m^3)^2 = Ksp (mol/kg-SW)^2 * rho_sw^2
    rho_sw = calc_seawater_density(T_c, S)  # kg/m^3
    Ksp_m3 = Ksp_mol_kg * (rho_sw**2)

    return CarbonateConstants(
        temp_c=T_c,
        temp_k=T_k,
        salinity=S,
        pressure_bar=P_bar,
        rho_sw=rho_sw,
        K1=float(K1_mol_kg),
        K2=float(K2_mol_kg),
        pK1=float(-np.log10(K1_mol_kg)),
        pK2=float(-np.log10(K2_mol_kg)),
        Ksp_calcite=float(Ksp_m3),
        K1_mol_kg=float(K1_mol_kg),
        K2_mol_kg=float(K2_mol_kg),
        Ksp_mol_kg=float(Ksp_mol_kg),
        ph_scale=scale_norm,
        scale_factor=float(scale_factor),
    )


def solve_carbonate_speciation_from_dic_ph(
    dic: float, ph: float, const: CarbonateConstants
) -> Dict[str, float]:
    """Calculates full carbonate speciation given DIC and pH."""
    h_plus = 10.0 ** (-float(ph))
    K1 = const.K1
    K2 = const.K2
    dic_val = float(dic)

    denom = h_plus**2 + K1 * h_plus + K1 * K2
    co2 = dic_val * (h_plus**2) / denom
    hco3 = dic_val * (K1 * h_plus) / denom
    co3 = dic_val * (K1 * K2) / denom
    ta = hco3 + 2.0 * co3

    return {
        "DIC": dic_val,
        "pH": float(ph),
        "Hplus": float(h_plus),
        "TA": float(ta),
        "CO2": float(co2),
        "HCO3": float(hco3),
        "CO3": float(co3),
    }


def solve_carbonate_speciation_from_ta_dic(
    ta: float, dic: float, const: CarbonateConstants
) -> Dict[str, float]:
    """Calculates full carbonate speciation given TA and DIC by solving the exact quadratic."""
    ta_val = float(ta)
    dic_val = float(dic)
    K1 = const.K1
    K2 = const.K2

    a = max(ta_val, 1e-12)
    b = K1 * (ta_val - dic_val)
    c_term = K1 * K2 * (ta_val - 2.0 * dic_val)

    disc = max(b**2 - 4.0 * a * c_term, 0.0)
    h_plus = (-b + np.sqrt(disc)) / (2.0 * a)
    h_plus = max(h_plus, 1e-14)
    ph = float(-np.log10(h_plus))

    denom = h_plus**2 + K1 * h_plus + K1 * K2
    co2 = dic_val * (h_plus**2) / denom
    hco3 = dic_val * (K1 * h_plus) / denom
    co3 = dic_val * (K1 * K2) / denom

    return {
        "DIC": dic_val,
        "pH": ph,
        "Hplus": float(h_plus),
        "TA": ta_val,
        "CO2": float(co2),
        "HCO3": float(hco3),
        "CO3": float(co3),
    }


def solve_carbonate_speciation_from_ta_ph(
    ta: float, ph: float, const: CarbonateConstants
) -> Dict[str, float]:
    """Calculates full carbonate speciation given TA and pH."""
    ta_val = float(ta)
    h_plus = 10.0 ** (-float(ph))
    K1 = const.K1
    K2 = const.K2

    denom = h_plus**2 + K1 * h_plus + K1 * K2
    alpha1 = K1 * h_plus / denom
    alpha2 = K1 * K2 / denom
    denom_ta = alpha1 + 2.0 * alpha2
    dic_val = ta_val / max(denom_ta, 1e-14)

    co2 = dic_val * (h_plus**2) / denom
    hco3 = dic_val * (K1 * h_plus) / denom
    co3 = dic_val * (K1 * K2) / denom

    return {
        "DIC": float(dic_val),
        "pH": float(ph),
        "Hplus": float(h_plus),
        "TA": ta_val,
        "CO2": float(co2),
        "HCO3": float(hco3),
        "CO3": float(co3),
    }


def adjust_for_calcite_equilibrium(
    ca2_top: float,
    dic: Optional[float] = None,
    ta: Optional[float] = None,
    ph: Optional[float] = None,
    const: Optional[CarbonateConstants] = None,
) -> Tuple[Dict[str, float], str]:
    """Adjusts carbonate boundary conditions to achieve exact calcite equilibrium (Omega = 1.0).

    Thermodynamic Principle:
    -------------------------
    At calcite equilibrium:
        Omega = [Ca2+] * [CO32-] / K'sp = 1.0
        ==> [CO32-]_eq = K'sp / [Ca2+]

    Supported Combinations (exit if Ca2 is not defined or <= 0):
    1. pH & DIC given  -> adjust TA
    2. pH & TA given   -> adjust DIC
    3. DIC & TA given  -> adjust pH (and align TA/DIC to exact equilibrium)

    Returns:
    --------
    Tuple[Dict[str, float], str]:
        (speciation_dict, human_readable_adjustment_description)
    """
    if ca2_top is None or float(ca2_top) <= 0.0:
        raise ValueError(
            "[MANTIS ERROR] Carbonate equilibrium chemistry requires species 'Ca2' "
            "with a positive concentration at the upper boundary to compute calcite saturation state Omega."
        )

    if const is None:
        raise ValueError("[MANTIS ERROR] CarbonateConstants must be provided.")

    ca2_val = float(ca2_top)
    co3_eq = const.Ksp_calcite / ca2_val
    K1 = const.K1
    K2 = const.K2

    has_ph = ph is not None
    has_dic = dic is not None
    has_ta = ta is not None

    # Case 1: pH & DIC given -> adjust TA
    if has_ph and has_dic and not has_ta:
        dic_val = float(dic)
        # alpha2 required to produce co3_eq
        alpha2_target = co3_eq / max(dic_val, 1e-12)
        if alpha2_target >= 1.0:
            raise ValueError(
                f"[MANTIS ERROR] DIC={dic_val} mol/m^3 is too low to satisfy [CO32-]_eq={co3_eq:.4e} mol/m^3."
            )
        # Solve for H+ from alpha2: K1*K2 / (h+^2 + K1*h+ + K1*K2) = alpha2
        a = alpha2_target
        b = alpha2_target * K1
        c_term = (alpha2_target - 1.0) * K1 * K2
        disc = max(b**2 - 4.0 * a * c_term, 0.0)
        h_plus_eq = (-b + np.sqrt(disc)) / (2.0 * a)
        ph_eq = float(-np.log10(h_plus_eq))

        denom = h_plus_eq**2 + K1 * h_plus_eq + K1 * K2
        co2_eq = dic_val * (h_plus_eq**2) / denom
        hco3_eq = dic_val * (K1 * h_plus_eq) / denom
        ta_eq = float(hco3_eq + 2.0 * co3_eq)

        spec = {
            "DIC": dic_val,
            "pH": ph_eq,
            "Hplus": float(h_plus_eq),
            "TA": ta_eq,
            "CO2": float(co2_eq),
            "HCO3": float(hco3_eq),
            "CO3": float(co3_eq),
        }
        msg = (
            f"Given pH={ph:.2f} & DIC={dic_val:.2f} mol/m^3: "
            f"Adjusted TA to {ta_eq:.3f} mol/m^3 (and equilibrium pH to {ph_eq:.2f}) "
            f"to achieve Omega_calcite = 1.00 ([CO32-] = {co3_eq:.4f} mol/m^3)."
        )
        return spec, msg

    # Case 2: pH & TA given -> adjust DIC
    elif has_ph and has_ta and not has_dic:
        ta_val = float(ta)
        # From TA = [HCO3-] + 2[CO32-]eq  ==> [HCO3-]eq = TA - 2[CO32-]eq
        hco3_eq = ta_val - 2.0 * co3_eq
        if hco3_eq <= 0:
            raise ValueError(
                f"[MANTIS ERROR] TA={ta_val} mol/m^3 is too low for [CO32-]_eq={co3_eq:.4e} mol/m^3."
            )
        # From K2 = h+ * [CO32-] / [HCO3-]  ==> h+ = K2 * [HCO3-] / [CO32-]
        h_plus_eq = K2 * hco3_eq / co3_eq
        ph_eq = float(-np.log10(h_plus_eq))

        # co2 = h+ * [HCO3-] / K1
        co2_eq = h_plus_eq * hco3_eq / K1
        dic_eq = float(co2_eq + hco3_eq + co3_eq)

        spec = {
            "DIC": dic_eq,
            "pH": ph_eq,
            "Hplus": float(h_plus_eq),
            "TA": ta_val,
            "CO2": float(co2_eq),
            "HCO3": float(hco3_eq),
            "CO3": float(co3_eq),
        }
        msg = (
            f"Given pH={ph:.2f} & TA={ta_val:.2f} mol/m^3: "
            f"Adjusted DIC to {dic_eq:.3f} mol/m^3 (and equilibrium pH to {ph_eq:.2f}) "
            f"to achieve Omega_calcite = 1.00 ([CO32-] = {co3_eq:.4f} mol/m^3)."
        )
        return spec, msg

    # Case 3: DIC & TA given -> adjust pH
    elif has_dic and has_ta and not has_ph:
        dic_val = float(dic)
        ta_orig = float(ta)
        # alpha2 required to produce co3_eq
        alpha2_target = co3_eq / max(dic_val, 1e-12)
        if alpha2_target >= 1.0:
            raise ValueError(
                f"[MANTIS ERROR] DIC={dic_val} mol/m^3 is too low to satisfy [CO32-]_eq={co3_eq:.4e} mol/m^3."
            )
        a = alpha2_target
        b = alpha2_target * K1
        c_term = (alpha2_target - 1.0) * K1 * K2
        disc = max(b**2 - 4.0 * a * c_term, 0.0)
        h_plus_eq = (-b + np.sqrt(disc)) / (2.0 * a)
        ph_eq = float(-np.log10(h_plus_eq))

        denom = h_plus_eq**2 + K1 * h_plus_eq + K1 * K2
        co2_eq = dic_val * (h_plus_eq**2) / denom
        hco3_eq = dic_val * (K1 * h_plus_eq) / denom
        ta_eq = float(hco3_eq + 2.0 * co3_eq)

        spec = {
            "DIC": dic_val,
            "pH": ph_eq,
            "Hplus": float(h_plus_eq),
            "TA": ta_eq,
            "CO2": float(co2_eq),
            "HCO3": float(hco3_eq),
            "CO3": float(co3_eq),
        }
        msg = (
            f"Given DIC={dic_val:.2f} & TA={ta_orig:.2f} mol/m^3: "
            f"Adjusted pH to {ph_eq:.2f} (and TA from {ta_orig:.2f} to {ta_eq:.3f} mol/m^3) "
            f"to achieve Omega_calcite = 1.00 ([CO32-] = {co3_eq:.4f} mol/m^3)."
        )
        return spec, msg

    else:
        raise ValueError(
            "[MANTIS ERROR] Invalid combination in carbonate_equilibrium_chemistry. "
            "Please specify exactly two of (pH, DIC, TA)."
        )
