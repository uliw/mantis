"""species_parser.py

Strict parser and validator for Mantis species definitions and isotope standards.
Enforces the modern declarative schema without legacy backwards compatibility:
- isotope_standards must be keyed by delta notation (e.g. 'd34S', 'd18O')
- Boundary conditions must use 'bc_top' and canonical types ('Dirichlet', 'Neumann', 'Robin')
- Companion isotope species are declared via 'isotopes' and automatically synthesized.
"""

import sys
import copy
import inspect
import ast
import re
from pathlib import Path
from typing import Dict, Any, Tuple, Optional
import numpy as np
import pint

_default_ureg = pint.UnitRegistry()


CANONICAL_BC_TYPES = {"Dirichlet", "Neumann", "Robin"}


def _find_definition_location(var_name: str = "isotope_standards", key: Optional[str] = None) -> Optional[str]:
    """Finds the source file name and line number where var_name (and optionally key) was defined."""
    frame = inspect.currentframe()
    candidate_path: Optional[Path] = None
    fallback_loc: Optional[str] = None

    while frame:
        if "self" in frame.f_locals:
            s = frame.f_locals["self"]
            sp = getattr(s, "species_file", None) or getattr(s, "species_path", None)
            if sp:
                candidate_path = Path(sp)
                break
        for mod_name in ("mod_species", "_species_mod", "species_mod"):
            if mod_name in frame.f_locals:
                m = frame.f_locals[mod_name]
                if hasattr(m, "__file__") and m.__file__:
                    candidate_path = Path(m.__file__)
                    break
        if candidate_path:
            break
        fn = frame.f_code.co_filename
        if "species_parser.py" not in fn and fallback_loc is None:
            fallback_loc = f"{fn}:{frame.f_lineno}"
        frame = frame.f_back

    if candidate_path and candidate_path.exists():
        try:
            tree = ast.parse(candidate_path.read_text(encoding="utf-8", errors="ignore"), filename=str(candidate_path))
            for node in ast.walk(tree):
                target_name = None
                if isinstance(node, ast.Assign):
                    for t in node.targets:
                        if isinstance(t, ast.Name):
                            target_name = t.id
                elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                    target_name = node.target.id

                if target_name == var_name and isinstance(node.value, ast.Dict):
                    if key is not None:
                        for k in node.value.keys:
                            if isinstance(k, ast.Constant) and k.value == key:
                                return f"{candidate_path}:{k.lineno}"
                    return f"{candidate_path}:{node.lineno}"
        except Exception:
            pass
        return str(candidate_path)

    return fallback_loc


def validate_isotope_standards(isotope_standards: Any) -> Dict[str, Dict[str, Any]]:
    """Strictly validates the isotope_standards dictionary."""
    if not isinstance(isotope_standards, dict):
        loc = _find_definition_location("isotope_standards")
        loc_str = f" in {loc}" if loc else ""
        print(
            f"\n[MANTIS ERROR]{loc_str} 'isotope_standards' must be a dictionary. Got: {type(isotope_standards)}.\n",
            file=sys.stderr,
        )
        sys.exit(1)

    validated = {}
    for std_key, std_meta in isotope_standards.items():
        if not isinstance(std_meta, dict):
            loc = _find_definition_location("isotope_standards", std_key)
            loc_str = f" at {loc}" if loc else ""
            print(
                f"\n[MANTIS ERROR]{loc_str} Isotope standard definition for '{std_key}' must be a dictionary. Got: {type(std_meta)}.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        # Reject legacy single-element key without element attribute
        if std_key in ("S", "C", "N", "O") and "element" not in std_meta:
            loc = _find_definition_location("isotope_standards", std_key)
            loc_str = f" at {loc}" if loc else ""
            print(
                f"\n[MANTIS ERROR]{loc_str} Legacy isotope standard key '{std_key}' detected in isotope_standards.\n"
                f"Mantis requires isotope standards to be keyed by delta notation (e.g. 'd34S', 'd18O')\n"
                f"and include 'standard', 'ratio', 'unit', and 'element'.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        required_fields = ["standard", "ratio", "unit", "element"]
        missing = [f for f in required_fields if f not in std_meta]
        if missing:
            loc = _find_definition_location("isotope_standards", std_key)
            loc_str = f" at {loc}" if loc else ""
            print(
                f"\n[MANTIS ERROR]{loc_str} Isotope standard '{std_key}' is missing required fields: {missing}.\n"
                f"Expected 'standard', 'ratio', 'unit', and 'element'. Got: {std_meta}\n",
                file=sys.stderr,
            )
            sys.exit(1)

        try:
            ratio_val = float(std_meta["ratio"])
        except (ValueError, TypeError):
            loc = _find_definition_location("isotope_standards", std_key)
            loc_str = f" at {loc}" if loc else ""
            print(
                f"\n[MANTIS ERROR]{loc_str} Isotope standard '{std_key}' has invalid 'ratio': {std_meta.get('ratio')}.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        validated[std_key] = {
            "standard": str(std_meta["standard"]),
            "ratio": ratio_val,
            "unit": str(std_meta["unit"]),
            "element": str(std_meta["element"]),
        }

    return validated


def resolve_isotope_standard(
    ref: str, isotope_standards: Dict[str, Dict[str, Any]]
) -> Tuple[str, Dict[str, Any]]:
    """Resolves an isotope standard by its delta key (e.g. 'd34S') or element (e.g. 'S')."""
    if ref in isotope_standards:
        return ref, isotope_standards[ref]

    # Match by element
    matched = [
        (k, v) for k, v in isotope_standards.items() if v.get("element") == ref
    ]
    if len(matched) == 1:
        return matched[0]
    elif len(matched) > 1:
        keys = [m[0] for m in matched]
        print(
            f"\n[MANTIS ERROR] Ambiguous isotope reference '{ref}'. Multiple isotope standards match element '{ref}': {keys}.\n"
            f"Please specify the explicit delta standard key (e.g. {keys}) instead of just the element symbol.\n",
            file=sys.stderr,
        )
        sys.exit(1)

    print(
        f"\n[MANTIS ERROR] Isotope reference '{ref}' not found in isotope_standards.\n"
        f"Available standards: {list(isotope_standards.keys())}.\n",
        file=sys.stderr,
    )
    sys.exit(1)


def parse_boundary_condition(
    bc_dict: Any,
    sp_name: str,
    location: str,
    ureg: Any = None,
) -> Tuple[str, float, Dict[str, Any]]:
    """Parses and strictly validates a boundary condition dictionary ('bc_top' or 'bc_bot').
    
    Returns:
        (bc_type, parsed_value, extra_properties)
    """
    if not isinstance(bc_dict, dict):
        print(
            f"\n[MANTIS ERROR] Boundary condition '{location}' for species '{sp_name}' must be a dictionary.\n"
            f"Got: {type(bc_dict)} ({bc_dict}).\n",
            file=sys.stderr,
        )
        sys.exit(1)

    # Reject legacy type specifiers
    if "type" in bc_dict or any(
        isinstance(v, str) and v.lower() in ("dissolved", "particulate")
        for v in bc_dict.values()
    ):
        print(
            f"\n[MANTIS ERROR] Species '{sp_name}' uses legacy boundary type syntax in '{location}'.\n"
            f"Legacy 'type': 'dissolved' / 'particulate' is deprecated and unsupported.\n"
            f"Please use canonical types as keys: 'Dirichlet', 'Neumann', or 'Robin'.\n"
            f"Example: '{location}': {{'Dirichlet': 28.2}}\n",
            file=sys.stderr,
        )
        sys.exit(1)

    canonical_keys = [k for k in bc_dict if k in CANONICAL_BC_TYPES]
    if len(canonical_keys) != 1:
        print(
            f"\n[MANTIS ERROR] Species '{sp_name}' '{location}' must specify exactly one canonical boundary condition type "
            f"('Dirichlet', 'Neumann', or 'Robin'). Got: {list(bc_dict.keys())}.\n",
            file=sys.stderr,
        )
        sys.exit(1)

    bc_type = canonical_keys[0]
    raw_val = bc_dict[bc_type]

    # Convert units if string with Pint
    if raw_val is None:
        val = None
    elif isinstance(raw_val, str):
        if ureg is None:
            ureg = _default_ureg
        try:
            q_val = ureg(raw_val)
            if any(time_u in raw_val for time_u in ("yr", "year", "second", "s", "day", "hour")):
                val = float(q_val.to("mol/(m^2 * second)").magnitude)
            else:
                val = float(q_val.to("mol/m^3").magnitude)
        except Exception as e:
            print(
                f"\n[MANTIS ERROR] Failed to parse unit string '{raw_val}' in '{location}' for '{sp_name}': {e}\n",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        try:
            val = float(raw_val)
        except (ValueError, TypeError):
            print(
                f"\n[MANTIS ERROR] Value for {bc_type} in '{location}' for species '{sp_name}' must be numeric or a unit string. Got: {raw_val}\n",
                file=sys.stderr,
            )
            sys.exit(1)

    extras = {k: v for k, v in bc_dict.items() if k != bc_type}
    return bc_type, val, extras


def get_l_mass(bulk: float, delta: float, standard_ratio: float) -> float:
    """Computes light isotope mass from bulk mass, delta (mUr/permil), and standard ratio."""
    r_val = (delta / 1000.0 + 1.0) * standard_ratio
    return bulk / (1.0 + r_val)


def parse_and_validate_species_dict(
    raw_species: Dict[str, Any],
    isotope_standards: Dict[str, Any],
    ureg: Any = None,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Dict[str, Any]]]:
    """Strictly parses species definitions, validates boundary conditions, and synthesizes companion species.
    
    Returns:
        (species_defs, bc_map, validated_isotope_standards)
    """
    if not isinstance(raw_species, dict):
        print(
            f"\n[MANTIS ERROR] 'species' must be a dictionary. Got: {type(raw_species)}.\n",
            file=sys.stderr,
        )
        sys.exit(1)

    val_stds = validate_isotope_standards(isotope_standards)

    # 1. First pass: check for deprecated legacy keys
    for sp_name, meta in raw_species.items():
        if not isinstance(meta, dict):
            print(
                f"\n[MANTIS ERROR] Species definition for '{sp_name}' must be a dictionary. Got: {type(meta)}.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        if "bc" in meta:
            print(
                f"\n[MANTIS ERROR] Species '{sp_name}' uses legacy 'bc' syntax.\n"
                f"Please update to 'bc_top': {{'Dirichlet'|'Neumann'|'Robin': value, ...}}.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        legacy_iso_keys = [k for k in ("isotope", "delta_name", "isotope_of") if k in meta]
        if legacy_iso_keys:
            print(
                f"\n[MANTIS ERROR] Species '{sp_name}' defines deprecated legacy isotope keys: {legacy_iso_keys}.\n"
                f"Isotope companions must now be declared inside the parent species via 'isotopes': {{'<companion_name>': '<element_or_delta_key>'}}.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        legacy_spec_keys = [k for k in ("partition", "equilibrium_of", "alias_of") if k in meta]
        if legacy_spec_keys:
            print(
                f"\n[MANTIS ERROR] Species '{sp_name}' declares legacy speciation keys: {legacy_spec_keys}.\n"
                f"Speciation must now be declared in the top-level 'speciation' dictionary in species.py.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        inc_ta = meta.get("include_in_TA", meta.get("include_in_ta"))
        if inc_ta is not None and not isinstance(inc_ta, (int, float, np.number)):
            print(
                f"\n[MANTIS ERROR] Species '{sp_name}' has invalid 'include_in_TA': {inc_ta}.\n"
                f"Expected a numerical factor (e.g., 1.0, 2.0, -1.0).\n",
                file=sys.stderr,
            )
            sys.exit(1)

    # 2. Build full species dictionary with synthesized companions
    species_defs = copy.deepcopy(raw_species)
    synthesized_companions: Dict[str, Dict[str, Any]] = {}

    for sp_name, meta in raw_species.items():
        iso_dict = meta.get("isotopes")
        if not iso_dict:
            continue

        if not isinstance(iso_dict, dict):
            print(
                f"\n[MANTIS ERROR] 'isotopes' in species '{sp_name}' must be a dictionary mapping companion name to element or standard key.\n"
                f"Got: {type(iso_dict)} ({iso_dict}).\n",
                file=sys.stderr,
            )
            sys.exit(1)

        for comp_name, iso_ref in iso_dict.items():
            std_key, std_meta = resolve_isotope_standard(iso_ref, val_stds)

            # Synthesize full companion species definition
            comp_entry: Dict[str, Any] = {
                "phase": meta.get("phase", "liquid"),
                "include": meta.get("include", True),
                "formula": meta.get("formula", sp_name),
                "diffusion": meta.get("diffusion"),
                "diff_coeff_function": meta.get("diff_coeff_function"),
                "isotope_of": sp_name,
                "isotope_standard": std_key,
                "element": std_meta["element"],
                "prefill": meta.get("prefill", False),
            }

            # Inherit equilibrium relation if parent is an equilibrium species
            if "equilibrium_of" in meta:
                parent_base = meta["equilibrium_of"]
                parent_meta = raw_species.get(parent_base, {})
                parent_iso_dict = parent_meta.get("isotopes", {})
                # Find parent companion tracking the same standard
                parent_comp = None
                for p_comp, p_ref in parent_iso_dict.items():
                    p_std, _ = resolve_isotope_standard(p_ref, val_stds)
                    if p_std == std_key:
                        parent_comp = p_comp
                        break
                if parent_comp:
                    comp_entry["equilibrium_of"] = parent_comp
                else:
                    comp_entry["equilibrium_of"] = parent_base
                comp_entry["fraction"] = meta.get("fraction")
                if "constants" in meta:
                    comp_entry["constants"] = meta["constants"]

            synthesized_companions[comp_name] = comp_entry

    species_defs.update(synthesized_companions)

    # 3. Parse and construct boundary condition map (bc_map)
    bc_map: Dict[str, Dict[str, Any]] = {}

    for sp_name, meta in species_defs.items():
        is_transported = meta.get("include", True)
        is_companion = "isotope_of" in meta

        # Process top boundary condition
        if "bc_top" in meta:
            bc_type, top_val, extras = parse_boundary_condition(
                meta["bc_top"], sp_name, "bc_top", ureg=ureg
            )
            if is_transported and top_val is None:
                print(
                    f"\n[MANTIS ERROR] Transported species '{sp_name}' must have a numerical boundary value in 'bc_top'. Got: None.\n",
                    file=sys.stderr,
                )
                sys.exit(1)
            bc_entry: Dict[str, Any] = {
                "type": bc_type,
                "top": top_val,
                "phase": meta.get("phase", "liquid"),
            }
            # Check for diffusive boundary layer (DBL)
            for k_dbl in ("dbl", "diffusive_boundary_layer", "dbl_thickness"):
                if k_dbl in extras:
                    bc_entry["dbl"] = float(extras[k_dbl])
            # Preserve delta values for isotope companions
            for k_extra, v_extra in extras.items():
                if k_extra in val_stds:
                    bc_entry[k_extra] = float(v_extra)

            bc_map[sp_name] = bc_entry
        elif is_transported and not is_companion:
            print(
                f"\n[MANTIS ERROR] Species '{sp_name}' has include=True but is missing 'bc_top'.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        # Process bottom boundary condition (defaults to Neumann 0.0 if omitted)
        if "bc_bot" in meta or "bc_bottom" in meta:
            bot_dict = meta.get("bc_bot", meta.get("bc_bottom"))
            bot_type, bot_val, _ = parse_boundary_condition(
                bot_dict, sp_name, "bc_bot", ureg=ureg
            )
            if sp_name in bc_map:
                bc_map[sp_name]["bot_type"] = bot_type
                bc_map[sp_name]["bot"] = bot_val
        elif is_transported and not is_companion:
            if sp_name in bc_map:
                bc_map[sp_name]["bot_type"] = "Neumann"
                bc_map[sp_name]["bot"] = 0.0

    # 4. Synthesize boundary conditions for isotope companion species
    for comp_name, comp_entry in synthesized_companions.items():
        parent_name = comp_entry["isotope_of"]
        parent_bc = bc_map.get(parent_name)
        if not parent_bc:
            continue

        std_key = comp_entry["isotope_standard"]
        std_ratio = val_stds[std_key]["ratio"]
        parent_meta = raw_species.get(parent_name, {})
        stoich = float(parent_meta.get("stoichiometry", 1.0))

        # Extract delta value from parent bc_top or default 0.0
        delta_val = float(parent_bc.get(std_key, 0.0))
        top_val_bulk = parent_bc.get("top", 0.0)
        if top_val_bulk is None:
            top_val_bulk = 0.0

        top_val_comp = get_l_mass(top_val_bulk * stoich, delta_val, std_ratio)

        comp_bc_entry: Dict[str, Any] = {
            "type": parent_bc["type"],
            "top": top_val_comp,
            "phase": comp_entry.get("phase", "liquid"),
            "bot_type": parent_bc.get("bot_type", "Neumann"),
            "bot": parent_bc.get("bot", 0.0),
        }
        if "dbl" in parent_bc:
            comp_bc_entry["dbl"] = parent_bc["dbl"]

        bc_map[comp_name] = comp_bc_entry

    return species_defs, bc_map, val_stds


VALID_SPECIATION_TYPES = {
    "carbonate_system",
    "carbonate",
    "acid_base",
    "partition",
    "complexation",
    "custom",
}


def validate_and_normalize_speciation(
    raw_speciation: Optional[Dict[str, Any]],
    species_defs: Dict[str, Any],
    ureg: Any = None,
) -> Dict[str, Dict[str, Any]]:
    """Validates explicit speciation blocks and normalizes legacy species speciation.

    Supports unified declarative speciation blocks:
      - 'carbonate_system': thermodynamic seawater carbonate speciation (DIC + TA -> CO3, HCO3, CO2, Hplus, pH)
      - 'acid_base': acid-base conjugate systems (e.g. TS2 -> HS, H2S)
      - 'partition': solid-liquid phase partitioning of lumped pools (e.g. Fe2_total -> dissolved, sorbed)
      - 'complexation': aqueous ligand complexation and free ion fractions (e.g. Fe2_liq, FeHS+)

    Ensures 100% backwards compatibility:
      Any legacy speciation constructs declared directly on species (e.g. 'partition',
      'equilibrium_of', 'fraction', 'include_in_TA') are automatically synthesized and
      normalized into the unified speciation dictionary.
    """
    speciation_defs: Dict[str, Dict[str, Any]] = {}

    # 1. Process explicit speciation dictionary if provided
    if raw_speciation is not None:
        if not isinstance(raw_speciation, dict):
            print(
                f"\n[MANTIS ERROR] 'speciation' must be a dictionary. Got: {type(raw_speciation)}.\n",
                file=sys.stderr,
            )
            sys.exit(1)

        for sys_name, sys_cfg in raw_speciation.items():
            if not isinstance(sys_cfg, dict):
                print(
                    f"\n[MANTIS ERROR] Speciation system '{sys_name}' must be a dictionary. Got: {type(sys_cfg)}.\n",
                    file=sys.stderr,
                )
                sys.exit(1)

            stype = sys_cfg.get("type", "acid_base")
            if stype not in VALID_SPECIATION_TYPES:
                print(
                    f"\n[MANTIS ERROR] Speciation system '{sys_name}' has invalid type '{stype}'.\n"
                    f"Valid types: {sorted(VALID_SPECIATION_TYPES)}.\n",
                    file=sys.stderr,
                )
                sys.exit(1)

            # Normalize canonical type
            if stype == "carbonate":
                stype = "carbonate_system"

            norm_sys = copy.deepcopy(sys_cfg)
            norm_sys["type"] = stype

            # If partition fractions are declared, propagate to primary species if not already set
            fractions = norm_sys.get("fractions", {})
            if isinstance(fractions, dict):
                primary_sp = norm_sys.get("primary") or norm_sys.get("parent")
                if isinstance(primary_sp, str) and primary_sp in species_defs:
                    if "partition" not in species_defs[primary_sp]:
                        species_defs[primary_sp]["partition"] = copy.deepcopy(fractions)

            # Register any sub-species declared in this speciation block into species_defs
            sub_sp = norm_sys.get("species", {})
            primary_val = norm_sys.get("primary", norm_sys.get("parent"))
            eq_parent = primary_val if isinstance(primary_val, str) else None

            # Support equilibrium_fractionation (preferred) or isotope_fractionation (alias)
            eq_frac_cfg = norm_sys.get("equilibrium_fractionation", norm_sys.get("isotope_fractionation"))

            if isinstance(sub_sp, dict):
                for sp_name, sp_meta in sub_sp.items():
                    if isinstance(sp_meta, dict):
                        if sp_name not in species_defs:
                            species_defs[sp_name] = {
                                "phase": sp_meta.get("phase", "liquid"),
                                "include": False,
                                "formula": sp_meta.get("formula", sp_name),
                                "equilibrium_of": eq_parent,
                            }
                        if "fraction" in sp_meta and sp_meta["fraction"] is not None:
                            species_defs[sp_name]["fraction"] = sp_meta["fraction"]
                        if "phase" in sp_meta:
                            species_defs[sp_name]["phase"] = sp_meta["phase"]
                        if "alias_of" in sp_meta:
                            species_defs[sp_name]["alias_of"] = sp_meta["alias_of"]
                        elif norm_sys.get("type") == "partition" and eq_parent:
                            if "diss" in sp_name.lower() or "liq" in sp_name.lower():
                                species_defs[sp_name].setdefault("alias_of", f"{eq_parent}.dissolved")
                            elif "sorb" in sp_name.lower():
                                species_defs[sp_name].setdefault("alias_of", f"{eq_parent}.sorbed")
                        if "include_in_TA" in sp_meta:
                            species_defs[sp_name]["include_in_TA"] = float(sp_meta["include_in_TA"])
                        if "fraction" in sp_meta:
                            species_defs[sp_name]["fraction"] = sp_meta["fraction"]
                        if "equilibrium_of" not in species_defs[sp_name] and eq_parent:
                            species_defs[sp_name]["equilibrium_of"] = eq_parent
                        if "constants" in norm_sys and "constants" not in species_defs[sp_name]:
                            species_defs[sp_name]["constants"] = copy.deepcopy(norm_sys["constants"])
                        if eq_frac_cfg and "equilibrium_fractionation" not in species_defs[sp_name]:
                            species_defs[sp_name]["equilibrium_fractionation"] = copy.deepcopy(eq_frac_cfg)
                        if "isotopes" in sp_meta and isinstance(sp_meta["isotopes"], dict):
                            species_defs[sp_name]["isotopes"] = copy.deepcopy(sp_meta["isotopes"])
                            for iso_sp_name, iso_ref in sp_meta["isotopes"].items():
                                if iso_sp_name not in species_defs:
                                    # Find parent companion tracking the same element/standard
                                    parent_base = eq_parent
                                    parent_meta = species_defs.get(parent_base, {}) if parent_base else {}
                                    parent_iso_dict = parent_meta.get("isotopes", {})
                                    parent_comp = None
                                    for p_comp, p_ref in parent_iso_dict.items():
                                        if p_ref == iso_ref:
                                            parent_comp = p_comp
                                            break
                                    species_defs[iso_sp_name] = {
                                        "phase": sp_meta.get("phase", "liquid"),
                                        "include": False,
                                        "formula": sp_meta.get("formula", sp_name),
                                        "isotope_of": sp_name,
                                        "equilibrium_of": parent_comp or parent_base,
                                        "fraction": sp_meta.get("fraction"),
                                    }
                                    if "constants" in norm_sys:
                                        species_defs[iso_sp_name]["constants"] = copy.deepcopy(norm_sys["constants"])
                                    if eq_frac_cfg:
                                        species_defs[iso_sp_name]["equilibrium_fractionation"] = copy.deepcopy(eq_frac_cfg)

            # Validate equilibrium_fractionation species keys if specified as a dict
            if eq_frac_cfg is not None:
                alpha_val = eq_frac_cfg.get("alpha") if isinstance(eq_frac_cfg, dict) else eq_frac_cfg
                if isinstance(alpha_val, dict):
                    for a_key in alpha_val.keys():
                        clean_k = str(a_key).strip()
                        # Allow format "A_to_B", "A / B", "A/B"
                        tokens = None
                        if "_to_" in clean_k:
                            tokens = clean_k.split("_to_")
                        elif " / " in clean_k:
                            tokens = clean_k.split(" / ")
                        elif "/" in clean_k:
                            tokens = clean_k.split("/")

                        if tokens and len(tokens) == 2:
                            s1, s2 = tokens[0].strip(), tokens[1].strip()
                            # Check that both s1 and s2 exist in species_defs
                            for s in (s1, s2):
                                if s not in species_defs:
                                    loc_str = f" in speciation system '{sys_name}'"
                                    print(
                                        f"\n[MANTIS ERROR]{loc_str} 'equilibrium_fractionation' references unknown species or isotope companion '{s}'.\n"
                                        f"Key '{clean_k}' must refer to recognized species declared in the model (e.g. 'H2S_S32_to_HS_S32').\n"
                                        f"Available species: {sorted(list(species_defs.keys()))}.\n",
                                        file=sys.stderr,
                                    )
                                    sys.exit(1)
            elif isinstance(sub_sp, (list, tuple)):
                for sp_name in sub_sp:
                    if sp_name not in species_defs:
                        species_defs[sp_name] = {
                            "phase": "liquid",
                            "include": False,
                            "formula": sp_name,
                            "equilibrium_of": eq_parent,
                        }

            speciation_defs[sys_name] = norm_sys

    return speciation_defs

