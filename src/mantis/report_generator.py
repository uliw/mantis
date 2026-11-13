"""MANTIS simulation report and metadata generator.

Generates human-readable Org-mode documents and machine-readable CSV files
documenting simulation parameters, biogeochemical reactions, rate laws,
boundary conditions, fluxes, kinetic constants, and elemental mass balances.
"""

from __future__ import annotations

import csv
from datetime import datetime
import os
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd

from mantis.diff_lib import get_flux


def stoichiometry_to_mhchem(rxn_str: str) -> str:
    """Convert reaction string (e.g. 'POC_fast + 4 Fe3 -> 4 Fe2_total') to mhchem format."""
    parts = rxn_str.split("->")
    if len(parts) != 2:
        parts = rxn_str.split("=")
    if len(parts) != 2:
        return f"\\ce{{{rxn_str}}}"

    def format_half(half: str) -> str:
        terms = [t.strip() for t in half.split("+") if t.strip()]
        formatted = []
        for term in terms:
            tokens = term.split()
            if len(tokens) == 2:
                coeff, name = tokens
            elif len(tokens) == 1:
                coeff, name = "", tokens[0]
            else:
                formatted.append(term)
                continue

            chem_map = {
                "POC_fast": "(CH2O)_f",
                "POC_slow": "(CH2O)_s",
                "POC": "CH2O",
                "Fe3": "FeOOH",
                "Fe2_total": "Fe^{2+}_{tot}",
                "Fe2_liq": "Fe^{2+}",
                "Fe2_sorbed": "X:Fe^{2+}",
                "Fe2": "Fe^{2+}",
                "SO4": "SO4^{2-}",
                "TS2": "H2S",
                "HS": "HS^-",
                "FeS": "FeS",
                "FeS2": "FeS2",
                "O2": "O2",
                "CO2": "CO2",
                "HCO3": "HCO3^-",
                "CO3": "CO3^{2-}",
                "H2O": "H2O",
                "Hplus": "H^+",
                "H": "H^+",
                "CH4": "CH4",
                "NH4": "NH4^+",
                "NO3": "NO3^-",
            }
            ch_name = chem_map.get(name, name)
            if coeff and coeff != "1":
                formatted.append(f"{coeff} {ch_name}")
            else:
                formatted.append(ch_name)
        return " + ".join(formatted)

    return f"$\\ce{{{format_half(parts[0])} -> {format_half(parts[1])}}}$"


def formula_to_elements(formula: str, name: str = "") -> dict[str, float]:
    """Parse elemental composition from chemical formula or species name."""
    overrides = {
        "POC_fast": {"C": 1.0, "H": 2.0, "O": 1.0},
        "POC_slow": {"C": 1.0, "H": 2.0, "O": 1.0},
        "POC": {"C": 1.0, "H": 2.0, "O": 1.0},
        "Fe3": {"Fe": 1.0, "O": 2.0, "H": 1.0},
        "FeOOH": {"Fe": 1.0, "O": 2.0, "H": 1.0},
        "Fe2_total": {"Fe": 1.0},
        "Fe2_liq": {"Fe": 1.0},
        "Fe2_sorbed": {"Fe": 1.0},
        "Fe2": {"Fe": 1.0},
        "SO4": {"S": 1.0, "O": 4.0},
        "TS2": {"S": 1.0, "H": 2.0},
        "H2S": {"S": 1.0, "H": 2.0},
        "HS": {"S": 1.0, "H": 1.0},
        "FeS": {"Fe": 1.0, "S": 1.0},
        "FeS2": {"Fe": 1.0, "S": 2.0},
        "O2": {"O": 2.0},
        "CO2": {"C": 1.0, "O": 2.0},
        "HCO3": {"H": 1.0, "C": 1.0, "O": 3.0},
        "CO3": {"C": 1.0, "O": 3.0},
        "H2O": {"H": 2.0, "O": 1.0},
        "H": {"H": 1.0},
        "Hplus": {"H": 1.0},
    }
    if name in overrides:
        return dict(overrides[name])
    clean = formula.replace("{", "").replace("}", "").strip()
    if clean in overrides:
        return dict(overrides[clean])

    tokens = re.findall(r"([A-Z][a-z]?)(?:_?(\d*\.?\d*))?", clean)
    counts: dict[str, float] = {}
    for el, num in tokens:
        if not el:
            continue
        n = float(num) if num else 1.0
        counts[el] = counts.get(el, 0.0) + n
    return counts


def generate_simulation_report(
    result: Any,
    output_dir: str | Path | None = None,
    save_csv_tables: bool = True,
    save_org_file: bool = True,
    save_txt_summary: bool = True,
) -> dict[str, Any]:
    """Generate simulation documentation reports and tables.

    Parameters
    ----------
    result : DiageneticModelResult
        MANTIS simulation result object.
    output_dir : str, Path, or None, optional
        Target directory to store logs and tables. If None, defaults to
        '<model_dir>/logs/'.
    save_csv_tables : bool, optional
        Whether to export individual CSV tables for machine learning/LLM access (default True).
    save_org_file : bool, optional
        Whether to export human-readable Org-mode report with LaTeX equations (default True).
    save_txt_summary : bool, optional
        Whether to export formatted plaintext summary similar to reactran 05_info.txt (default True).

    Returns
    -------
    dict
        Paths of created files and table dataframes.
    """
    mp = result.mp
    c = result.c
    k = result.k
    df = getattr(result, "df", None)
    if df is None:
        if callable(getattr(result, "diagenetic_reactions", None)):
            try:
                df, _ = result.save_csv()
            except Exception:
                df = None
        if df is None:
            # Construct a minimal DataFrame from c and z
            data = {"z": result.z}
            for s in getattr(result, "species_list", []):
                if hasattr(c, s):
                    v = getattr(c, s)
                    data[f"c_{s}"] = v.value if hasattr(v, "value") else v
            df = pd.DataFrame(data)

    base_dir = getattr(mp, "base_dir", None) or Path.cwd()
    if output_dir is None:
        logs_path = Path(base_dir) / "logs"
    else:
        logs_path = Path(output_dir)
        if not logs_path.is_absolute():
            logs_path = Path(base_dir) / logs_path

    logs_path.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # 1. Run Metadata
    # -------------------------------------------------------------------------
    now = datetime.now()
    run_date_str = getattr(mp, "run_start_time", now.strftime("%Y-%m-%d %H:%M:%S"))
    cpu_time = getattr(mp, "cpu_time", getattr(result, "cpu_time", None))
    wall_time = getattr(mp, "wall_time", getattr(result, "wall_time", None))
    if cpu_time is None:
        cpu_time = 0.0
    if wall_time is None:
        wall_time = 0.0

    def _clean_fname(p: str | None) -> str:
        if not p:
            return ""
        name = Path(p).name
        return f"~{name}~" if name else ""

    input_files = {
        "species_file": _clean_fname(str(getattr(mp, "species_file", getattr(mp, "_orig_species", "")))),
        "reactions_file": _clean_fname(str(getattr(mp, "reactions_file", getattr(mp, "_orig_reactions", "")))),
        "constants_file": _clean_fname(str(getattr(mp, "constants_file", getattr(mp, "_orig_constants", "")))),
        "equations_file": _clean_fname(str(getattr(mp, "equations_file", getattr(mp, "_orig_equations", "")))),
        "layout_file": _clean_fname(str(getattr(mp, "layout_file", getattr(mp, "_orig_layout", "")))),
    }
    output_files = {
        "csv_path": _clean_fname(str(getattr(result, "fqfn", ""))),
        "state_archive": _clean_fname(f"{mp.plot_name}_state.npz"),
        "pdf_plot": _clean_fname(f"{mp.plot_name}.pdf"),
    }

    meta_rows = [
        {"Parameter": "Simulation Name", "Value": str(getattr(mp, "plot_name", "mantis_sim"))},
        {"Parameter": "Run Date", "Value": str(run_date_str)},
        {"Parameter": "CPU Time", "Value": f"{cpu_time:.3f} s"},
        {"Parameter": "Wall Clock Time", "Value": f"{wall_time:.3f} s"},
        {"Parameter": "Total Simulation Time", "Value": f"{getattr(result, 'total_time', 0.0):.2e} s"},
        {"Parameter": "Total Timesteps", "Value": str(getattr(result, "step", 0))},
        {"Parameter": "Convergence Status", "Value": str(getattr(result, "converged", "Yes"))},
        {"Parameter": "Grid Points (N)", "Value": str(len(result.z))},
        {"Parameter": "Domain Length", "Value": f"{float(result.z[-1]):.4f} m"},
    ]
    for ik, iv in input_files.items():
        if iv:
            meta_rows.append({"Parameter": f"Input: {ik}", "Value": iv})
    for ok, ov in output_files.items():
        if ov:
            meta_rows.append({"Parameter": f"Output: {ok}", "Value": ov})

    df_meta = pd.DataFrame(meta_rows)

    # -------------------------------------------------------------------------
    # 2. Table 1: Model & Environmental Parameters
    # -------------------------------------------------------------------------
    def _extract_val(v: Any, default: Any = None) -> Any:
        if v is None:
            return default
        if hasattr(v, "value"):
            v = v.value
        if isinstance(v, (np.ndarray, list)):
            arr = np.asarray(v)
            if arr.ndim > 0 and len(arr) > 1:
                return f"{float(arr[0]):.4g} to {float(arr[-1]):.4g}"
            return float(arr.flat[0])
        return v

    phi_val = _extract_val(getattr(mp, "phi", 0.8))
    w_val_raw = getattr(mp, "w", 0.0)
    w_yr = float(_extract_val(w_val_raw, 0.0)) * 365.25 * 86400 * 100  # m/s -> cm/yr
    bt0_raw = getattr(mp, "BT0", 0.0)
    bt0_yr = float(_extract_val(bt0_raw, 0.0)) * 365.25 * 86400 * 10000  # m^2/s -> cm^2/yr
    bi0_raw = getattr(mp, "BI0", 0.0)
    bi0_yr = float(_extract_val(bi0_raw, 0.0)) * 365.25 * 86400  # 1/s -> 1/yr

    temp_val = _extract_val(getattr(mp, "temp", getattr(mp, "TC", 10.0)))
    if isinstance(temp_val, list):
        temp_str = f"{temp_val[0]} to {temp_val[-1]}"
    else:
        temp_str = str(temp_val)

    table1_rows = [
        {"Parameter": "Salinity", "Symbol": "$S$", "Depth Profile": "Constant", "Value": str(getattr(mp, "S", "35")), "Unit": "–"},
        {"Parameter": "Temperature", "Symbol": "$T_C$", "Depth Profile": "Linear / Constant", "Value": temp_str, "Unit": "°C"},
        {"Parameter": "Pressure", "Symbol": "$P$", "Depth Profile": "Constant", "Value": str(getattr(mp, "P", "1.013")), "Unit": "bar"},
        {"Parameter": "Porewater pH", "Symbol": r"$\mathrm{pH}$", "Depth Profile": "Constant", "Value": str(getattr(mp, "pH", "7.5")), "Unit": "–"},
        {"Parameter": "Porosity", "Symbol": r"$\phi$", "Depth Profile": "Profile / Constant", "Value": str(phi_val), "Unit": "–"},
        {"Parameter": "Sediment Dry Density", "Symbol": r"$\rho_{\mathrm{sed}}$", "Depth Profile": "Constant", "Value": str(getattr(mp, "rho_sed", "2.6")), "Unit": r"$\mathrm{g\,cm^{-3}}$"},
        {"Parameter": "Sedimentation Velocity", "Symbol": "$w$", "Depth Profile": "Constant / Compaction", "Value": f"{w_yr:.4g}", "Unit": r"$\mathrm{cm\,yr^{-1}}$"},
        {"Parameter": "Bioturbation Coeff at SWI", "Symbol": "$D_{b,0}$", "Depth Profile": str(getattr(mp, "BT_profile", "sigmoidal")), "Value": f"{bt0_yr:.4g}", "Unit": r"$\mathrm{cm^2\,yr^{-1}}$"},
        {"Parameter": "Bioturbation Mixing Depth", "Symbol": "$x_L$", "Depth Profile": "–", "Value": f"{float(getattr(mp, 'BT_depth', 0.06))*100:.2f}", "Unit": "cm"},
        {"Parameter": "Bioturbation Attenuation", "Symbol": "$x_{\\mathrm{bm}}$", "Depth Profile": "–", "Value": f"{float(getattr(mp, 'BT_attenuation', 0.02))*100:.2f}", "Unit": "cm"},
        {"Parameter": "Bio-irrigation Coeff at SWI", "Symbol": r"$\alpha_0$", "Depth Profile": "exponential", "Value": f"{bi0_yr:.4g}", "Unit": r"$\mathrm{yr^{-1}}$"},
        {"Parameter": "Bio-irrigation Attenuation", "Symbol": "$x_{\\mathrm{irr}}$", "Depth Profile": "–", "Value": f"{float(getattr(mp, 'BI_depth', 0.03))*100:.2f}", "Unit": "cm"},
    ]
    df_table1 = pd.DataFrame(table1_rows)

    # -------------------------------------------------------------------------
    # 3. Table 2: Biogeochemical Reactions
    # -------------------------------------------------------------------------
    raw_rxns = getattr(mp, "reactions", [])
    if not raw_rxns and hasattr(mp, "reactions_file") and mp.reactions_file and Path(mp.reactions_file).exists():
        import importlib.util

        spec = importlib.util.spec_from_file_location("rxn_temp", mp.reactions_file)
        if spec and spec.loader:
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            raw_rxns = getattr(mod, "reactions", [])

    def _format_k_val(k_val: Any) -> str:
        if isinstance(k_val, dict):
            return ", ".join(f"~{v}~" for v in k_val.values())
        elif isinstance(k_val, (list, tuple, set)):
            return ", ".join(f"~{v}~" for v in k_val)
        elif k_val:
            return f"~{k_val}~"
        return "–"

    table2_rows = []
    if raw_rxns:
        for idx, r_entry in enumerate(raw_rxns, 1):
            r_name = r_entry.get("reaction_name", f"reaction_{idx}")
            r_eq = r_entry.get("reaction", "")
            table2_rows.append({
                "Reaction Name": f"~{r_name}~",
                "Chemical Equation": stoichiometry_to_mhchem(r_eq),
                "Rate Constant": _format_k_val(r_entry.get("k_value_name")),
            })
    else:
        table2_rows.append({
            "Reaction Name": "–",
            "Chemical Equation": "–",
            "Rate Constant": "–",
        })
    df_table2 = pd.DataFrame(table2_rows)

    # -------------------------------------------------------------------------
    # 4. Table 3: Kinetic Limiters & Rate Laws
    # -------------------------------------------------------------------------
    raw_limiters = getattr(mp, "limiters", {})
    if not raw_limiters and hasattr(mp, "reactions_module") and hasattr(mp.reactions_module, "limiters"):
        raw_limiters = mp.reactions_module.limiters

    chem_vars = {
        "POC_fast": r"$\ce{[(CH2O)_f]}$",
        "POC_slow": r"$\ce{[(CH2O)_s]}$",
        "POC": r"$\ce{[CH2O]}$",
        "Fe3_diss_red": r"$\ce{[FeOOH]_{re}}$",
        "Fe3": r"$\ce{[FeOOH]}$",
        "Fe2_total": r"$\ce{[Fe^{2+}_{tot}]}$",
        "Fe2_liq": r"$\ce{[Fe^{2+}]}$",
        "Fe2_sorbed": r"$\ce{[X:Fe^{2+}]}$",
        "Fe2": r"$\ce{[Fe^{2+}]}$",
        "SO4": r"$\ce{[SO4^{2-}]}$",
        "TS2": r"$\ce{[H2S]}$",
        "HS": r"$\ce{[HS^-]}$",
        "FeS": r"$\ce{[FeS]}$",
        "FeS2": r"$\ce{[FeS2]}$",
        "O2": r"$\ce{[O2]}$",
        "CO2": r"$\ce{[CO2]}$",
        "HCO3": r"$\ce{[HCO3^-]}$",
        "Hplus": r"$\ce{[H^+]}$",
        "H": r"$\ce{[H^+]}$",
    }

    def _clean_rate_expr(expr_str: str) -> str:
        s = re.sub(r"\bnp\.maximum\b", "max", str(expr_str))
        for var, ch in chem_vars.items():
            s = re.sub(r"\bc\." + var + r"\b", lambda m, c=ch: c, s)
            s = re.sub(r"\bmp\." + var + r"\b", lambda m, c=ch: c, s)
        s = re.sub(r"\bc\.", "", s)
        s = re.sub(r"\bmp\.", "", s)
        return s

    table3_rows = []
    # Add limiters
    for lim_name, lim_expr in raw_limiters.items():
        table3_rows.append({
            "Type": "Limiter",
            "Target / Reaction": f"~{lim_name}~",
            "Expression": _clean_rate_expr(lim_expr),
            "Notes": "Kinetic substrate limiter or inhibition factor",
        })

    # Add reaction rate links
    for idx, r_entry in enumerate(raw_rxns, 1):
        r_name = r_entry.get("reaction_name", f"R{idx}")
        lims = r_entry.get("limiters", [])
        inhibs = r_entry.get("inhibitors", [])
        lim_str = ", ".join(f"~{l}~" for l in lims) if lims else "None"
        inhib_str = ", ".join(f"~{i}~" for i in inhibs) if inhibs else "None"
        table3_rows.append({
            "Type": "Rate Law",
            "Target / Reaction": f"~{r_name}~",
            "Expression": f"k * [Substrates] * L({lim_str}) * I({inhib_str})",
            "Notes": f"Limiters: {lim_str}; Inhibitors: {inhib_str}",
        })
    df_table3 = pd.DataFrame(table3_rows)

    # -------------------------------------------------------------------------
    # 5. Table 4: Boundary Conditions & Fluxes
    # -------------------------------------------------------------------------
    species_defs = getattr(mp, "species_defs", None) or getattr(mp, "species", None) or {}
    bc_map = getattr(mp, "bc_map", {})
    all_species = getattr(result, "species_list", list(species_defs.keys()))

    table4_rows = []
    fluxes_up: dict[str, float] = {}
    fluxes_down: dict[str, float] = {}

    for sp in all_species:
        sp_meta = species_defs.get(sp, {})
        sp_phase = sp_meta.get("phase", "solid" if getattr(result, "is_solid", False) else "liquid")
        bc_info = bc_map.get(sp, {})

        top_type = bc_info.get("type", "Dirichlet")
        top_val = bc_info.get("top", 0.0)
        bot_type = bc_info.get("bot_type", "Neumann")
        bot_val = bc_info.get("bot", 0.0)

        # Determine units for top/bot BC
        top_unit = "umol/(cm^2*yr)" if top_type in ("Robin", "particulate") else ("mol/m^3" if top_type in ("Dirichlet", "dissolved") else "mol/m^4")
        bot_unit = "umol/(cm^2*yr)" if bot_type in ("Robin", "particulate") else ("mol/m^3" if bot_type in ("Dirichlet", "dissolved") else "mol/m^4")

        # Calculate fluxes at top (depth=0) and bottom (depth=None)
        try:
            f_up = get_flux(sp, depth=0.0, df=df, c=c, mp=mp, z=result.z, D_mol=result.D_mol, units="umol/(cm^2*yr)")
        except Exception:
            f_up = 0.0
        try:
            f_down = get_flux(sp, depth=None, df=df, c=c, mp=mp, z=result.z, D_mol=result.D_mol, units="umol/(cm^2*yr)")
        except Exception:
            f_down = 0.0

        fluxes_up[sp] = f_up
        fluxes_down[sp] = f_down

        # Concentrations at boundaries
        c_col = f"c_{sp}" if f"c_{sp}" in df else (sp if sp in df else None)
        c_up_val = float(df[c_col].iloc[0]) if c_col else 0.0
        c_down_val = float(df[c_col].iloc[-1]) if c_col else 0.0

        table4_rows.append({
            "Species": sp,
            "Phase": sp_phase,
            "Top BC Type": top_type,
            "Top BC Value": f"{top_val:.4g}",
            "Top BC Unit": top_unit,
            "C_up (mol/m^3)": f"{c_up_val:.4g}",
            "Flux_up (umol/cm^2/yr)": f"{f_up:.4g}",
            "Bot BC Type": bot_type,
            "Bot BC Value": f"{bot_val:.4g}",
            "Bot BC Unit": bot_unit,
            "C_down (mol/m^3)": f"{c_down_val:.4g}",
            "Flux_down (umol/cm^2/yr)": f"{f_down:.4e}",
        })
    df_table4 = pd.DataFrame(table4_rows)

    # -------------------------------------------------------------------------
    # 6. Table 5: Kinetic Constants & Rate Coefficients
    # -------------------------------------------------------------------------
    table5_rows = []
    # If kinetic_constants file provided raw dict
    raw_k_dict = getattr(mp, "raw_reaction_constants", {})
    if not raw_k_dict and hasattr(mp, "constants_module") and hasattr(mp.constants_module, "get_reaction_constants"):
        try:
            raw_k_dict = mp.constants_module.get_reaction_constants(getattr(mp, "pH", 7.5), getattr(mp, "phi", 0.8))
        except Exception:
            raw_k_dict = {}

    if hasattr(k, "items"):
        k_items = list(k.items())
    elif isinstance(k, dict):
        k_items = list(k.items())
    else:
        k_items = [(k_name, getattr(k, k_name)) for k_name in dir(k) if not k_name.startswith("_") and not callable(getattr(k, k_name))]

    for k_name, si_val in k_items:
        orig_val = raw_k_dict.get(k_name, "–")
        table5_rows.append({
            "Constant": k_name,
            "Original Value / Unit": str(orig_val),
            "SI Value (m-s-mol)": f"{float(si_val):.4e}" if isinstance(si_val, (int, float, np.floating)) else str(si_val),
        })
    df_table5 = pd.DataFrame(table5_rows)

    # -------------------------------------------------------------------------
    # 7. Table 6: Elemental Mass Balance
    # -------------------------------------------------------------------------
    # Elements to balance: C, S, Fe
    # Note: O and H are components of bulk solvent H2O (55.5 mol/L porewater), which is not an
    # actively solved transport variable. Thus O and H are excluded from transport flux balance.
    elements = ["C", "S", "Fe"]
    elem_up: dict[str, float] = {e: 0.0 for e in elements}
    elem_down: dict[str, float] = {e: 0.0 for e in elements}

    for sp, f_up in fluxes_up.items():
        sp_meta = species_defs.get(sp, {})
        # If species is an isotope companion (e.g. SO4_S32), skip to avoid double counting base element
        if "isotope_of" in sp_meta or sp.endswith("_S32"):
            continue
        formula = sp_meta.get("formula", sp)
        el_counts = formula_to_elements(formula, sp)
        for el, count in el_counts.items():
            if el in elem_up:
                elem_up[el] += count * f_up

    for sp, f_down in fluxes_down.items():
        sp_meta = species_defs.get(sp, {})
        if "isotope_of" in sp_meta or sp.endswith("_S32"):
            continue
        formula = sp_meta.get("formula", sp)
        el_counts = formula_to_elements(formula, sp)
        for el, count in el_counts.items():
            if el in elem_down:
                elem_down[el] += count * f_down

    table6_rows = []
    for el in elements:
        up_val = elem_up[el]
        down_val = elem_down[el]
        imbalance = up_val - down_val
        rel_imb = (imbalance / (abs(up_val) + 1e-30)) * 100.0
        table6_rows.append({
            "Element": el,
            "Up Flux (umol/cm^2/yr)": f"{up_val:.4g}",
            "Down Flux (umol/cm^2/yr)": f"{down_val:.4g}",
            "Imbalance (umol/cm^2/yr)": f"{imbalance:.4e}",
            "Unit": "umol/(cm^2*yr)",
            "Imbalance (%)": f"{rel_imb:.2f}%",
        })
    df_table6 = pd.DataFrame(table6_rows)

    # -------------------------------------------------------------------------
    # 8. Table 7: Isotope Fractionation and Equilibrium Effects
    # -------------------------------------------------------------------------
    table7_rows: list[dict[str, Any]] = []
    has_isotopes = bool(getattr(mp, "isotopes", False))
    if has_isotopes:
        rxn_params = getattr(mp, "reaction_parameters", {})
        if not isinstance(rxn_params, dict):
            rxn_params = {}

        R_std = getattr(mp, "VCDT", getattr(mp, "d34S_ratio", 0.044162589))

        # 1. Kinetic fractionation from reactions
        reactions_list = getattr(mp, "reactions", [])
        seen_processes = set()

        for r in reactions_list:
            iso_frac = r.get("isotope_fractionation")
            if not iso_frac or not isinstance(iso_frac, dict):
                continue
            rxn_name = r.get("reaction_name", "")
            proc_name = rxn_name.replace("_fast", "").replace("_slow", "").replace("_", " ").strip()
            if proc_name in seen_processes:
                continue
            seen_processes.add(proc_name)

            for f_key, f_val in iso_frac.items():
                alpha_val = None
                if isinstance(f_val, (list, tuple)) and len(f_val) > 0:
                    alpha_ref = f_val[0]
                    if isinstance(alpha_ref, (int, float)):
                        alpha_val = float(alpha_ref)
                    elif isinstance(alpha_ref, str):
                        if hasattr(mp, alpha_ref):
                            alpha_val = float(getattr(mp, alpha_ref))
                        else:
                            alpha_val = float(rxn_params.get(alpha_ref, 1.0))
                elif isinstance(f_val, (int, float)):
                    alpha_val = float(f_val)

                if alpha_val is not None:
                    eps = (1.0 - alpha_val) * 1000.0
                    eps_str = f"{eps:.0f} mUr" if abs(eps - round(eps)) < 1e-4 else f"{eps:.1f} mUr"
                    table7_rows.append({
                        "Process": proc_name,
                        "Enrichment factor": eps_str,
                        "Type": "kinetic",
                    })

        # 2. Equilibrium fractionation from speciation
        speciation_dict = getattr(mp, "speciation", {})
        if isinstance(speciation_dict, dict):
            for sys_name, sys_cfg in speciation_dict.items():
                eq_frac = sys_cfg.get("equilibrium_fractionation")
                if not eq_frac or not isinstance(eq_frac, dict):
                    continue
                proc_name = "sulfide speciation" if "sulfide" in sys_name.lower() else sys_name.replace("_", " ").strip()
                alpha_entry = eq_frac.get("alpha", {})
                if isinstance(alpha_entry, dict):
                    for a_k, a_v in alpha_entry.items():
                        a_val = float(a_v)
                        # Check whether parameter is 32S formulation (a_val ~ 0.99991542) or heavy isotope alpha (e.g. 1.002)
                        clean_ak = str(a_k)
                        if ("32" in clean_ak or abs(a_val - 1.0) < 0.0005) and a_val != 1.0:
                            R_frac = R_std / (1.0 + R_std)
                            eps_34 = (1.0 - a_val) / R_frac * 1000.0
                            eps_str = f"{round(eps_34):+d} mUr" if abs(eps_34 - round(eps_34)) < 0.05 else f"{eps_34:+.1f} mUr"
                        else:
                            eps = (a_val - 1.0) * 1000.0
                            eps_str = f"{eps:+.0f} mUr" if abs(eps - round(eps)) < 1e-4 else f"{eps:+.1f} mUr"
                        table7_rows.append({
                            "Process": proc_name,
                            "Enrichment factor": eps_str,
                            "Type": "equilibrium",
                        })

    df_table7 = pd.DataFrame(table7_rows if table7_rows else [{"Process": "–", "Enrichment factor": "–", "Type": "–"}])

    # -------------------------------------------------------------------------
    # Save CSV tables
    # -------------------------------------------------------------------------
    csv_paths: dict[str, Path] = {}
    if save_csv_tables:
        csv_paths["01_simulation_metadata"] = logs_path / "01_simulation_metadata.csv"
        csv_paths["02_parameters"] = logs_path / "02_parameters.csv"
        csv_paths["03_reactions"] = logs_path / "03_reactions.csv"
        csv_paths["04_kinetic_rate_laws"] = logs_path / "04_kinetic_rate_laws.csv"
        csv_paths["05_boundary_conditions_and_fluxes"] = logs_path / "05_boundary_conditions_and_fluxes.csv"
        csv_paths["06_kinetic_constants"] = logs_path / "06_kinetic_constants.csv"
        csv_paths["07_elemental_mass_balance"] = logs_path / "07_elemental_mass_balance.csv"
        if has_isotopes and table7_rows:
            csv_paths["08_isotope_fractionation"] = logs_path / "08_isotope_fractionation.csv"

        df_meta.to_csv(csv_paths["01_simulation_metadata"], index=False)
        df_table1.to_csv(csv_paths["02_parameters"], index=False)
        df_table2.to_csv(csv_paths["03_reactions"], index=False)
        df_table3.to_csv(csv_paths["04_kinetic_rate_laws"], index=False)
        df_table4.to_csv(csv_paths["05_boundary_conditions_and_fluxes"], index=False)
        df_table5.to_csv(csv_paths["06_kinetic_constants"], index=False)
        df_table6.to_csv(csv_paths["07_elemental_mass_balance"], index=False)
        if has_isotopes and table7_rows:
            df_table7.to_csv(csv_paths["08_isotope_fractionation"], index=False)

    # -------------------------------------------------------------------------
    # Save Plaintext Summary (05_info.txt style)
    # -------------------------------------------------------------------------
    txt_path = logs_path / f"{mp.plot_name}_info.txt"
    if save_txt_summary:
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write("======================================================================\n")
            f.write(f"              MANTIS SIMULATION REPORT: {mp.plot_name}\n")
            f.write("======================================================================\n")
            f.write(f"  Run Date        : {run_date_str}\n")
            f.write(f"  CPU Time        : {cpu_time:.3f} s\n")
            f.write(f"  Wall Time       : {wall_time:.3f} s\n")
            f.write(f"  Grid Points (N) : {len(result.z)}\n")
            f.write(f"  Domain Depth    : {float(result.z[-1]):.4f} m\n")
            f.write(f"  Convergence     : {getattr(result, 'converged', 'Yes')}\n\n")

            if has_isotopes and table7_rows:
                f.write("======================================================================\n")
                f.write("             ISOTOPE FRACTIONATION & EQUILIBRIUM EFFECTS\n")
                f.write("======================================================================\n")
                f.write(f"  {'Process':<30} | {'Enrichment factor':<20} | {'Type':<12}\n")
                f.write("  " + "-" * 66 + "\n")
                for r in table7_rows:
                    f.write(f"  {r['Process']:<30} | {r['Enrichment factor']:<20} | {r['Type']:<12}\n")
                f.write("\n")

            f.write("======================================================================\n")
            f.write("                       ELEMENTAL MASS BALANCE\n")
            f.write("======================================================================\n")
            for r in table6_rows:
                f.write(
                    f"  Element {r['Element']:<4} | Up Flux: {r['Up Flux (umol/cm^2/yr)']:<10} | "
                    f"Down Flux: {r['Down Flux (umol/cm^2/yr)']:<10} | "
                    f"Imbalance: {r['Imbalance (umol/cm^2/yr)']} {r['Unit']} ({r['Imbalance (%)']})\n"
                )
            f.write("\n======================================================================\n")
            f.write("                SPECIES BOUNDARY CONDITIONS & FLUXES\n")
            f.write("======================================================================\n")
            for r in table4_rows:
                f.write(
                    f"  {r['Species']:<12} ({r['Phase']:<6}) | C.up: {r['C_up (mol/m^3)']:<9} mol/m^3 | "
                    f"C.down: {r['C_down (mol/m^3)']:<9} mol/m^3 | Flux.up: {r['Flux_up (umol/cm^2/yr)']:<10} umol/(cm^2*yr) | "
                    f"Flux.down: {r['Flux_down (umol/cm^2/yr)']} umol/(cm^2*yr)\n"
                )
            f.write("\n======================================================================\n")
            f.write("                          REACTION NETWORK\n")
            f.write("======================================================================\n")
            f.write(f"  Reactions modeled ({len(table2_rows)}):\n")
            for r in table2_rows:
                f.write(f"  {r['Reaction Name']:<35} : {r['Chemical Equation']}\n")
            f.write("\n")

    # -------------------------------------------------------------------------
    # Save Org-Mode Document with Dynamic Python CSV Import
    # -------------------------------------------------------------------------
    org_path = logs_path / f"{mp.plot_name}_report.org"
    if save_org_file:
        def make_dynamic_org_table(block_name: str, csv_filename: str, caption: str) -> str:
            return (
                f"#+NAME: {block_name}\n"
                f"#+BEGIN_SRC python :results value :colnames no :exports results\n"
                f"import csv\n"
                f"with open('{csv_filename}', 'r', newline='') as f:\n"
                f"    rows = list(csv.reader(f))\n"
                f"return [rows[0], None] + rows[1:] if len(rows) > 1 else rows\n"
                f"#+END_SRC\n"
                f"#+attr_latex: :font \\footnotesize :placement [t]\n"
                f"#+caption: {caption}\n"
                f"#+RESULTS: {block_name}\n"
            )

        with open(org_path, "w", encoding="utf-8") as f:
            f.write(f"#+TITLE: MANTIS Diagenetic Simulation Report: {mp.plot_name}\n")
            f.write(f"#+DATE: {run_date_str}\n")
            f.write("#+OPTIONS: toc:nil num:nil ^:nil\n")
            f.write("#+LATEX_HEADER: \\usepackage{amsmath,amssymb,booktabs}\n")
            f.write("#+LATEX_HEADER: \\usepackage[version=4]{mhchem}\n")
            f.write("#+LATEX_HEADER: \\usepackage{array,ragged2e}\n")
            f.write("#+LATEX_HEADER: \\newcolumntype{L}[1]{>{\\RaggedRight\\arraybackslash}p{#1}}\n")
            f.write("#+LATEX_HEADER: \\newcolumntype{C}[1]{>{\\Centering\\arraybackslash}p{#1}}\n")
            f.write("#+LATEX_HEADER: \\newcolumntype{R}[1]{>{\\RaggedLeft\\arraybackslash}p{#1}}\n\n")

            f.write(make_dynamic_org_table(
                "meta",
                "01_simulation_metadata.csv",
                "Simulation metadata and runtime performance.",
            ))
            f.write("\n")

            f.write(make_dynamic_org_table(
                "parameters",
                "02_parameters.csv",
                "Model and environmental parameters (cf. Van de Velde et al. 2016 Table 1).",
            ))
            f.write("\n")

            f.write(
                "#+NAME: reactions\n"
                "#+BEGIN_SRC python :results value :colnames no :exports results\n"
                "import csv\n"
                "with open('03_reactions.csv', 'r', newline='') as f:\n"
                "    rows = list(csv.reader(f))\n"
                "return [rows[0], None] + rows[1:] if len(rows) > 1 else rows\n"
                "#+END_SRC\n"
                "#+attr_latex: :font \\footnotesize :placement [t] :align L{40mm}L{60mm}L{40mm}\n"
                "#+caption: Biogeochemical reaction network and stoichiometry (cf. Van de Velde et al. 2016 Table 2).\n"
                "#+RESULTS: reactions\n"
            )
            f.write("\n")

            f.write(
                "#+NAME: rate_laws\n"
                "#+BEGIN_SRC python :results value :colnames no :exports results\n"
                "import csv\n"
                "with open('04_kinetic_rate_laws.csv', 'r', newline='') as f:\n"
                "    rows = list(csv.reader(f))\n"
                "# Drop Notes column (column index 3)\n"
                "cleaned = [[r[0], r[1], r[2]] for r in rows]\n"
                "return [cleaned[0], None] + cleaned[1:] if len(cleaned) > 1 else cleaned\n"
                "#+END_SRC\n"
                "#+attr_latex: :font \\footnotesize :placement [t] :align L{25mm}L{35mm}L{80mm}\n"
                "#+caption: Kinetic rate formulations and substrate limiters (cf. Van de Velde et al. 2016 Table 3).\n"
                "#+RESULTS: rate_laws\n"
            )
            f.write("\n")

            f.write(
                "#+NAME: boundary_fluxes\n"
                "#+BEGIN_SRC python :results value :colnames no :exports results\n"
                "import csv\n"
                "with open('05_boundary_conditions_and_fluxes.csv', 'r', newline='') as f:\n"
                "    rows = list(csv.reader(f))\n"
                "# Keep only the first 5 columns and drop 'Top BC ' from headers\n"
                "h_clean = [h.replace('Top BC ', '') for h in rows[0][:5]]\n"
                "cleaned = [h_clean, None] + [r[:5] for r in rows[1:]]\n"
                "return cleaned if len(rows) > 1 else rows\n"
                "#+END_SRC\n"
                "#+attr_latex: :font \\footnotesize :placement [t]\n"
                "#+caption: Boundary conditions and concentrations (cf. Van de Velde et al. 2016 Table 4).\n"
                "#+RESULTS: boundary_fluxes\n"
            )
            f.write("\n")

            f.write(make_dynamic_org_table(
                "constants",
                "06_kinetic_constants.csv",
                "Kinetic constants and parameter values (cf. Van de Velde et al. 2016 Table 5).",
            ))
            f.write("\n")

            if has_isotopes and table7_rows:
                f.write(make_dynamic_org_table(
                    "isotopes",
                    "08_isotope_fractionation.csv",
                    "Isotope fractionation and equilibrium effects.",
                ))
                f.write("\n")

    return {
        "logs_dir": logs_path,
        "csv_paths": csv_paths,
        "org_path": org_path if save_org_file else None,
        "txt_path": txt_path if save_txt_summary else None,
        "tables": {
            "metadata": df_meta,
            "parameters": df_table1,
            "reactions": df_table2,
            "rate_laws": df_table3,
            "boundary_fluxes": df_table4,
            "kinetic_constants": df_table5,
            "mass_balance": df_table6,
            "isotope_fractionation": df_table7 if (has_isotopes and table7_rows) else None,
        },
    }
