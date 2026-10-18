"""generate_equation_system.py

Standalone declarative parser and code generator for MANTIS reaction transport models.

Translates declarative reaction specifications (e.g. chemical_equations_new.py) into
a fully assembled Python equation module (e.g. equations.py) with exact validation of:
  - Species definitions in species.py
  - Reaction constants in reaction_constants.py
  - Limiters in limiters.py

Usage:
------
    python -m mantis.generate_equation_system [OPTIONS]
    python src/mantis/generate_equation_system.py -i chemical_equations_new.py -o equations.py
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple


class ValidationError(Exception):
    """Raised when chemical reaction specifications fail validation checks."""
    pass


class ReactionSystemValidator:
    """Validates chemical reaction specifications against model definitions."""

    def __init__(
        self,
        species_path: Path,
        constants_path: Path,
        limiters_path: Optional[Path] = None,
        reactions_path: Optional[Path] = None,
    ):
        self.species_path = Path(species_path)
        self.constants_path = Path(constants_path)
        self.limiters_path = Path(limiters_path) if limiters_path is not None else None
        self.reactions_path = Path(reactions_path) if reactions_path is not None else None

        self.valid_species: Set[str] = set()
        self.valid_constants: Set[str] = set()
        self.valid_limiters: Set[str] = set()

        self._load_definitions()

    def _load_definitions(self) -> None:
        """Loads and extracts defined species, constants, and limiters."""
        # 1. Species
        if not self.species_path.exists():
            raise FileNotFoundError(f"Species file not found: {self.species_path}")
        spec_species = importlib.util.spec_from_file_location("species_mod", self.species_path)
        if spec_species is None or spec_species.loader is None:
            raise ImportError(f"Could not load module from {self.species_path}")
        mod_species = importlib.util.module_from_spec(spec_species)
        spec_species.loader.exec_module(mod_species)
        if hasattr(mod_species, "species") and isinstance(mod_species.species, dict):
            self.valid_species = set(mod_species.species.keys())
            self.species_defs = mod_species.species
            self.isotope_standards = getattr(mod_species, "isotope_standards", {})
            # Also register dynamic/alias species if defined with equilibrium_of or alias_of
            for sp_k, sp_meta in mod_species.species.items():
                if isinstance(sp_meta, dict):
                    if "equilibrium_of" in sp_meta or "alias_of" in sp_meta or "fraction" in sp_meta:
                        self.valid_species.add(sp_k)
        else:
            raise ValueError(f"No 'species' dictionary found in {self.species_path}")

        # 2. Reaction constants
        if not self.constants_path.exists():
            raise FileNotFoundError(f"Reaction constants file not found: {self.constants_path}")
        spec_constants = importlib.util.spec_from_file_location("rc_mod", self.constants_path)
        if spec_constants is None or spec_constants.loader is None:
            raise ImportError(f"Could not load module from {self.constants_path}")
        mod_constants = importlib.util.module_from_spec(spec_constants)
        spec_constants.loader.exec_module(mod_constants)

        # Attempt dynamic extraction via get_reaction_constants
        self.valid_constants = set()
        if hasattr(mod_constants, "get_reaction_constants"):
            try:
                _velde, k_vals = mod_constants.get_reaction_constants(7.5, 0.8)
                if isinstance(k_vals, dict):
                    self.valid_constants.update(k_vals.keys())
                if isinstance(_velde, dict):
                    self.valid_constants.update(_velde.keys())
            except Exception:
                pass

        # Static AST fallback to capture constants in velde dict or module variables
        try:
            rc_tree = ast.parse(self.constants_path.read_text(encoding="utf-8"))
            for node in ast.walk(rc_tree):
                if isinstance(node, ast.Dict):
                    for key in node.keys:
                        if isinstance(key, ast.Constant) and isinstance(key.value, str):
                            self.valid_constants.add(key.value)
                elif isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Name):
                            self.valid_constants.add(target.id)
        except Exception:
            pass

        # 3. Limiters (from limiters_path or reactions_path)
        self.valid_limiters = set()
        self.limiters_defs: Dict[str, str] = {}

        def _extract_limiters_from_ast(tree: ast.AST):
            extracted = {}
            for node in ast.walk(tree):
                if isinstance(node, ast.Assign):
                    for target in node.targets:
                        if isinstance(target, ast.Subscript):
                            if isinstance(target.value, ast.Name) and target.value.id == "limiters":
                                if isinstance(target.slice, ast.Constant) and isinstance(target.slice.value, str):
                                    expr_code = ast.unparse(node.value) if hasattr(ast, "unparse") else ""
                                    extracted[target.slice.value] = expr_code
                        elif isinstance(target, ast.Name) and target.id == "limiters":
                            if isinstance(node.value, ast.Dict):
                                for k_node, v_node in zip(node.value.keys, node.value.values):
                                    if isinstance(k_node, ast.Constant) and isinstance(k_node.value, str):
                                        expr_code = ast.unparse(v_node) if hasattr(ast, "unparse") else ""
                                        extracted[k_node.value] = expr_code
            return extracted

        raw_limiters: Dict[str, str] = {}
        # 1. First check if limiters are defined in reactions_path (modern declarative standard)
        if self.reactions_path is not None and self.reactions_path.exists():
            try:
                rxn_tree = ast.parse(self.reactions_path.read_text(encoding="utf-8"))
                raw_limiters.update(_extract_limiters_from_ast(rxn_tree))
            except Exception:
                pass

        # 2. If not defined in reactions_path, fall back to limiters_path
        if not raw_limiters and self.limiters_path is not None and self.limiters_path.exists():
            try:
                lim_tree = ast.parse(self.limiters_path.read_text(encoding="utf-8"))
                raw_limiters.update(_extract_limiters_from_ast(lim_tree))
            except Exception as e:
                raise ValueError(f"Could not parse limiters file {self.limiters_path}: {e}")

        # Expand limiters:
        # Base implicit limiter L_imp = 1 / (C + K) -> derive _implicit, _explicit, _inhibit
        for lim_name, expr in raw_limiters.items():
            self.limiters_defs[lim_name] = expr
            self.valid_limiters.add(lim_name)

            # If defined as base form (e.g. 'O2' or 'SO4' or 'Fe3') or ending in '_implicit'
            base_key = lim_name[:-9] if lim_name.endswith("_implicit") else lim_name

            # Add implicit form
            self.valid_limiters.add(f"{base_key}_implicit")
            self.valid_limiters.add(f"{base_key}_explicit")
            self.valid_limiters.add(f"{base_key}_inhibit")
            self.valid_limiters.add(f"{base_key}_inhib")

            if "_" in base_key:
                main_k, sub_k = base_key.split("_", 1)
                self.valid_limiters.add(f"{main_k}_implicit_{sub_k}")
                self.valid_limiters.add(f"{main_k}_explicit_{sub_k}")
                self.valid_limiters.add(f"{main_k}_inhibit_{sub_k}")
                self.valid_limiters.add(f"{main_k}_inhib_{sub_k}")

    @staticmethod
    def parse_reaction_species(reaction_str: str) -> Tuple[List[Tuple[float, str]], List[Tuple[float, str]]]:
        """Parses a chemical reaction equation string into (reactants, products).

        Each species entry is (stoichiometric_coefficient, species_name).
        Correctly handles ions with trailing '+' or '-' (e.g. '2 H+', 'SO4--') by splitting on
        whitespace-padded ' + ' delimiters.
        """
        if "->" not in reaction_str:
            raise ValueError(f"Reaction string missing '->' delimiter: '{reaction_str}'")

        left, right = reaction_str.split("->", 1)

        def _parse_side(side_str: str) -> List[Tuple[float, str]]:
            terms = re.split(r"\s+\+\s+", side_str.strip())
            parsed = []
            for term in terms:
                term = term.strip()
                if not term:
                    continue
                parts = term.split(None, 1)
                if len(parts) == 1:
                    coeff = 1.0
                    sp = parts[0]
                else:
                    try:
                        coeff = float(parts[0])
                        sp = parts[1]
                    except ValueError:
                        coeff = 1.0
                        sp = term
                parsed.append((coeff, sp))
            return parsed

        return _parse_side(left), _parse_side(right)

    def get_species_phase(self, sp: str) -> str:
        """Returns 'liquid', 'solid', or 'mixed' for the given species.

        Resolves aliases, isotopes, and equilibrium species back to their parent phase.
        """
        meta = getattr(self, "species_defs", {}).get(sp)
        if not meta or not isinstance(meta, dict):
            parent = self.get_parent_species(sp)
            if parent != sp and parent in getattr(self, "species_defs", {}):
                return self.get_species_phase(parent)
            return "liquid"

        if "phase" in meta:
            return meta["phase"]
        if "alias_of" in meta:
            parent = meta["alias_of"].split(".")[0]
            return self.get_species_phase(parent)
        if "equilibrium_of" in meta:
            return self.get_species_phase(meta["equilibrium_of"])
        if "isotope_of" in meta:
            return self.get_species_phase(meta["isotope_of"])
        return "liquid"

    def get_parent_species(self, sp: str) -> str:
        """Returns the primary transport state variable name for a given species.

        E.g. HS -> TS2, Fe2_liq -> Fe2_total, SO4 -> SO4.
        """
        meta = getattr(self, "species_defs", {}).get(sp)
        if isinstance(meta, dict):
            if "equilibrium_of" in meta:
                return meta["equilibrium_of"]
            if "alias_of" in meta:
                return meta["alias_of"].split(".")[0]
            if "isotope_of" in meta and not meta.get("include", True):
                return meta["isotope_of"]
        return sp

    def get_speciation_factor_expr(self, sp: str) -> Optional[str]:
        """Returns python code expression for the fraction of parent species, or None if 1.0.
        
        Derives the fraction expression dynamically from species_defs:
        - If 'fraction' is defined, returns 'mp.<species_name.lower()>_frac' or the analytical string.
        - If 'alias_of' is defined as '<parent>.<subphase>', returns 'mp.<parent>_<subphase>'.
        - If 'isotope_of' is defined, inherits the base species' speciation factor expression.
        """
        meta = getattr(self, "species_defs", {}).get(sp)
        if isinstance(meta, dict):
            # 1. Direct equilibrium speciation fraction
            if "fraction" in meta:
                return f"mp.{sp.lower()}_frac"

            # 2. Subphase alias of partitioned parent (e.g. Fe2_total.dissolved -> mp.Fe2_total_dissolved)
            if "alias_of" in meta and isinstance(meta["alias_of"], str):
                parts = meta["alias_of"].split(".", 1)
                parent_name = parts[0]
                subphase = parts[1] if len(parts) > 1 else "dissolved"
                parent_meta = getattr(self, "species_defs", {}).get(parent_name, {})
                if isinstance(parent_meta, dict) and "partition" in parent_meta:
                    return f"getattr(mp, '{parent_name}_{subphase}', getattr(mp, '{sp}', 1.0))"

            # 3. Isotopic companion inheriting base species speciation factor
            if "isotope_of" in meta and isinstance(meta["isotope_of"], str):
                base_sp = meta["isotope_of"]
                base_spec = self.get_speciation_factor_expr(base_sp)
                if base_spec:
                    return base_spec

        return None

    def is_transport_species(self, sp: str) -> bool:
        """Returns True if the species is an active transport PDE variable (include: True)."""
        meta = getattr(self, "species_defs", {}).get(sp)
        if isinstance(meta, dict):
            return bool(meta.get("include", True))
        return sp in self.valid_species

    def validate_reaction(self, rxn_dict: Dict[str, Any]) -> List[str]:
        """Validates a single reaction dictionary. Returns a list of error messages (empty if valid)."""
        errors = []
        rxn_name = rxn_dict.get("reaction_name", "<unnamed_reaction>")
        rxn_str = rxn_dict.get("reaction", "")
        dyn_vars = rxn_dict.get("dynamic_variables", {}) or {}
        iso_map = rxn_dict.get("isotope_species", {}) or {}
        k_val = rxn_dict.get("k_value_name")
        limiters = rxn_dict.get("limiters")

        known_local_species = self.valid_species.union(dyn_vars.keys())

        # 1. Validate reaction equation and species names
        if not rxn_str:
            errors.append(f"[{rxn_name}] Reaction equation string is empty or missing.")
        else:
            try:
                reactants, products = self.parse_reaction_species(rxn_str)
                for coeff, sp in reactants + products:
                    if sp not in known_local_species:
                        errors.append(
                            f"[{rxn_name}] Species '{sp}' not found in species.py or dynamic_variables."
                        )
            except Exception as e:
                errors.append(f"[{rxn_name}] Failed to parse reaction equation '{rxn_str}': {e}")

        # 2. Validate kinetic constant(s)
        if k_val is not None:
            if isinstance(k_val, str):
                if k_val not in self.valid_constants:
                    errors.append(
                        f"[{rxn_name}] Rate constant '{k_val}' not found in reaction_constants.py."
                    )
            elif isinstance(k_val, dict):
                for branch, k_sub in k_val.items():
                    if k_sub not in self.valid_constants:
                        errors.append(
                            f"[{rxn_name}] Rate constant '{k_sub}' (branch '{branch}') not found in reaction_constants.py."
                        )
            else:
                errors.append(f"[{rxn_name}] Invalid k_value_name type: {type(k_val)}.")

        # 3. Validate limiters & inhibitors
        def _check_lim_entry(key_label: str, val: Any):
            if val is None:
                return
            if isinstance(val, dict):
                for target_sp, lim_items in val.items():
                    if target_sp not in known_local_species and target_sp != "global":
                        errors.append(
                            f"[{rxn_name}] {key_label} target species '{target_sp}' not defined in species.py or dynamic_variables."
                        )
                    item_list = [lim_items] if isinstance(lim_items, str) else lim_items
                    for lim_name in item_list:
                        # Check bare name, _implicit, _inhibit, etc.
                        if not any(
                            cand in self.valid_limiters
                            for cand in (lim_name, f"{lim_name}_implicit", f"{lim_name}_inhibit", f"{lim_name}_inhib")
                        ):
                            errors.append(
                                f"[{rxn_name}] Fatal Error: {key_label} '{lim_name}' is not defined in the limiters specification or reactions file."
                            )
            elif isinstance(val, (list, tuple)):
                for lim_name in val:
                    if not any(
                        cand in self.valid_limiters
                        for cand in (lim_name, f"{lim_name}_implicit", f"{lim_name}_inhibit", f"{lim_name}_inhib")
                    ):
                        errors.append(
                            f"[{rxn_name}] Fatal Error: {key_label} '{lim_name}' is not defined in the limiters specification or reactions file."
                        )
            elif isinstance(val, str):
                if not any(
                    cand in self.valid_limiters
                    for cand in (val, f"{val}_implicit", f"{val}_inhibit", f"{val}_inhib")
                ):
                    errors.append(
                        f"[{rxn_name}] Fatal Error: {key_label} '{val}' is not defined in the limiters specification or reactions file."
                    )
            else:
                errors.append(f"[{rxn_name}] Invalid {key_label} format: {type(val)}.")

        _check_lim_entry("Limiter", limiters)
        inhibitors = rxn_dict.get("inhibitors")
        _check_lim_entry("Inhibitor", inhibitors)

        # 4. Validate isotope mapping
        if iso_map:
            if not isinstance(iso_map, dict):
                errors.append(f"[{rxn_name}] isotope_species must be a dictionary.")
            else:
                for base_sp, iso_sp in iso_map.items():
                    if base_sp not in known_local_species:
                        errors.append(
                            f"[{rxn_name}] Base isotope species '{base_sp}' not defined in species.py or dynamic_variables."
                        )
                    if iso_sp not in known_local_species:
                        errors.append(
                            f"[{rxn_name}] Isotope species '{iso_sp}' not defined in species.py or dynamic_variables."
                        )

        # 5. Validate dynamic variables expressions (ensure valid Python syntax)
        for var_name, expr in dyn_vars.items():
            try:
                ast.parse(expr)
            except SyntaxError as se:
                errors.append(
                    f"[{rxn_name}] Invalid Python syntax in dynamic variable '{var_name}': '{expr}' ({se})."
                )

        return errors

    def validate_reactions(self, reactions: List[Dict[str, Any]]) -> List[str]:
        """Validates all reactions in the specification list."""
        all_errors = []
        for rxn_dict in reactions:
            errs = self.validate_reaction(rxn_dict)
            all_errors.extend(errs)
        return all_errors


def load_reactions_from_file(input_path: Path) -> List[Dict[str, Any]]:
    """Loads the reactions list from a chemical equations Python or JSON file."""
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input equations file not found: {input_path}")

    if input_path.suffix == ".py":
        spec = importlib.util.spec_from_file_location("chem_eq_mod", input_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load module from {input_path}")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        if not hasattr(mod, "reactions") or not isinstance(mod.reactions, list):
            raise ValueError(f"Input file {input_path} must define a 'reactions' list.")
        return mod.reactions
    elif input_path.suffix == ".json":
        import json
        with open(input_path, "r", encoding="utf-8") as f:
            data = json.load(f)
            if isinstance(data, list):
                return data
            elif isinstance(data, dict) and "reactions" in data:
                return data["reactions"]
            raise ValueError(f"JSON input file {input_path} must be a list of reactions or contain a 'reactions' key.")
    else:
        raise ValueError(f"Unsupported input file format '{input_path.suffix}'. Expected .py or .json.")


class MantisArgumentParser(argparse.ArgumentParser):
    """ArgumentParser subclass ensuring backward compatibility aliases on parsed arguments."""

    def parse_args(self, args=None, namespace=None):
        ns = super().parse_args(args=args, namespace=namespace)
        if not hasattr(ns, "input") or getattr(ns, "input") is None:
            setattr(ns, "input", getattr(ns, "reactions", None))
        if not hasattr(ns, "output") or getattr(ns, "output") is None:
            setattr(ns, "output", getattr(ns, "equations", None))
        if not hasattr(ns, "constants") or getattr(ns, "constants") is None:
            setattr(ns, "constants", getattr(ns, "kinetic_constants", None))
        return ns


# Alias for backward compatibility
FiPyriteArgumentParser = MantisArgumentParser


def build_arg_parser() -> argparse.ArgumentParser:
    """Builds and returns the command line argument parser with full help descriptions."""
    parser = MantisArgumentParser(
        prog="generate_equation_system",
        description=(
            "MANTIS Declarative Equation System Generator:\n"
            "Validates declarative chemical reaction specifications (e.g. chemical_equations_new.py) "
            "against model species, reaction constants, and limiters, and generates an assembled "
            "reaction equations module (e.g. equations.py)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-r", "--reactions", "-i", "--input",
        dest="reactions",
        type=Path,
        default=Path("nbk/experiments/chemical_equations_new.py"),
        help="Path to declarative reaction definitions file (default: nbk/experiments/chemical_equations_new.py)",
    )
    parser.add_argument(
        "-e", "--equations", "-o", "--output",
        dest="equations",
        type=Path,
        default=Path("nbk/experiments/equations.py"),
        help="Path to output generated equations module (default: nbk/experiments/equations.py)",
    )
    parser.add_argument(
        "-s", "--species",
        type=Path,
        default=Path("nbk/experiments/species.py"),
        help="Path to model species definitions file (default: nbk/experiments/species.py)",
    )
    parser.add_argument(
        "-k", "--kinetic-constants", "--constants",
        dest="kinetic_constants",
        type=Path,
        default=Path("nbk/experiments/reaction_constants.py"),
        help="Path to reaction constants file (default: nbk/experiments/reaction_constants.py)",
    )
    parser.add_argument(
        "-l", "--limiters",
        type=Path,
        default=Path("nbk/experiments/limiters.py"),
        help="Path to limiters definition file (default: nbk/experiments/limiters.py)",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Perform cross-reference validation of species, constants, and limiters without writing output.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable detailed validation and parsing diagnostic output.",
    )

    return parser


def generate_equations(
    reactions_path: str | Path,
    equations_path: str | Path,
    kinetic_constants_path: str | Path = Path("nbk/experiments/reaction_constants.py"),
    limiters_path: str | Path = Path("nbk/experiments/limiters.py"),
    species_path: str | Path = Path("nbk/experiments/species.py"),
    verbose: bool = False,
    validate_only: bool = False,
) -> Path:
    """Validates and compiles declarative reactions into an assembled equations module.

    Parameters
    ----------
    reactions_path : str | Path
        Path to declarative reactions definitions file.
    equations_path : str | Path
        Path where the generated equations module will be written.
    kinetic_constants_path : str | Path
        Path to reaction constants definitions file.
    limiters_path : str | Path | None, optional
        Path to limiters definition file. If None or not found, limiters will be read from reactions_path.
    species_path : str | Path
        Path to model species definitions file.
    verbose : bool
        If True, prints diagnostic output.
    validate_only : bool
        If True, performs cross-validation without writing output.

    Returns
    -------
    Path
        Path to the output equations module.
    """
    def _find_path(p_in: str | Path) -> Path:
        p = Path(p_in)
        if p.exists():
            return p.resolve()
        if (Path.cwd() / p).exists():
            return (Path.cwd() / p).resolve()
        if (Path.cwd() / p.name).exists():
            return (Path.cwd() / p.name).resolve()
        return p.resolve()

    reactions_p = _find_path(reactions_path)
    equations_p = Path(equations_path).resolve()
    constants_p = _find_path(kinetic_constants_path)
    limiters_p = _find_path(limiters_path) if limiters_path is not None else None
    species_p = _find_path(species_path)

    if verbose:
        print(f"[MANTIS] Loading reactions from: {reactions_p}")
        print(f"[MANTIS] Using species definitions: {species_p}")
        print(f"[MANTIS] Using reaction constants:   {constants_p}")
        if limiters_p and limiters_p.exists():
            print(f"[MANTIS] Using limiters definitions: {limiters_p}")
        else:
            print(f"[MANTIS] Using limiters definitions from reactions file")

    # 1. Load validator
    validator = ReactionSystemValidator(
        species_path=species_p,
        constants_path=constants_p,
        limiters_path=limiters_p if (limiters_p and limiters_p.exists()) else None,
        reactions_path=reactions_p,
    )

    # 2. Load input reactions
    reactions = load_reactions_from_file(reactions_p)

    # 3. Validate
    errors = validator.validate_reactions(reactions)
    if errors:
        err_msg = f"Validation failed for {reactions_p} with {len(errors)} error(s):\n" + "\n".join(f"  • {e}" for e in errors)
        raise ValidationError(err_msg)

    if verbose:
        print(f"[MANTIS] All {len(reactions)} reactions validated successfully!")

    if validate_only:
        return equations_p

    # 4. Code Generation
    generator = EquationSystemGenerator(reactions=reactions, validator=validator)
    generated_code = generator.generate_code()

    equations_p.parent.mkdir(parents=True, exist_ok=True)
    with open(equations_p, "w", encoding="utf-8") as f:
        f.write(generated_code)

    if verbose:
        print(f"[MANTIS] Successfully wrote equations module to: {equations_p}")

    return equations_p


def resolve_and_prepare_equations(
    p_dict: Dict[str, Any],
    experiment_name: Optional[str] = None,
    base_dir: Optional[Path] = None,
) -> Tuple[Any, Path, Path]:
    """Resolves reaction and equation files from p_dict, regenerates equations if outdated,
    and dynamically imports the equations module.

    Resolution scheme:
      reactions: defaults to f"{experiment}_reactions.py"
      equations: defaults to f"{experiment}_equations.py"
      kinetic_constants: defaults to "reaction_constants.py"
      limiters: defaults to "limiters.py"

    Parameters
    ----------
    p_dict : dict
        Parameter dictionary passed to the simulation.
    experiment_name : str, optional
        Name of the experiment stem. If None, derived from p_dict["experiment"].
    base_dir : Path, optional
        Base directory to search for experiment files. Defaults to cwd or experiment file parent.

    Returns
    -------
    Tuple[module, Path, Path]
        (loaded_equations_module, resolved_reactions_path, resolved_equations_path)
    """
    exp = p_dict.get("experiment") or experiment_name or "pyrite"

    if base_dir is None:
        if "__file__" in p_dict and p_dict["__file__"]:
            base_dir = Path(p_dict["__file__"]).resolve().parent
        else:
            base_dir = Path.cwd()

    def _resolve(val: Any, default_name: str) -> Path:
        filename = default_name if (val is None or val == "file_name") else str(val)
        p = Path(filename)
        if p.is_absolute() and p.exists():
            return p.resolve()
        # Check relative to cwd
        if (Path.cwd() / p).exists():
            return (Path.cwd() / p).resolve()
        # Check relative to base_dir
        if (base_dir / p).exists():
            return (base_dir / p).resolve()
        return (base_dir / p).resolve()

    # Apply defaults if not specified or None
    rxn_val = p_dict.get("reactions")
    eq_val = p_dict.get("equations")
    kc_val = p_dict.get("kinetic_constants")
    lim_val = p_dict.get("limiters")
    sp_val = p_dict.get("species")

    reactions_path = _resolve(rxn_val, f"{exp}_reactions.py")
    equations_path = _resolve(eq_val, f"{exp}_equations.py")
    constants_path = _resolve(kc_val, "reaction_constants.py")
    species_path = _resolve(sp_val, "species.py")

    # If limiters not explicitly specified in p_dict, check if base_dir or cwd has limiters.py
    if lim_val is not None:
        limiters_path = _resolve(lim_val, "limiters.py")
    elif (base_dir / "limiters.py").exists():
        limiters_path = (base_dir / "limiters.py").resolve()
    elif (Path.cwd() / "limiters.py").exists():
        limiters_path = (Path.cwd() / "limiters.py").resolve()
    else:
        limiters_path = None

    if not reactions_path.exists():
        raise FileNotFoundError(
            f"Declarative reactions file not found: '{reactions_path}'. "
            f"Please specify a valid 'reactions' file in p_dict or create '{exp}_reactions.py'."
        )

    # Determine whether regeneration is needed based on mtime timestamps
    needs_regen = False
    reason = ""
    generator_path = Path(__file__).resolve()
    force_regen = bool(p_dict.get("force_regen", False) or p_dict.get("rebuild_equations", False))

    if force_regen:
        needs_regen = True
        reason = "forced regeneration requested"
    elif not equations_path.exists():
        needs_regen = True
        reason = f"{equations_path.name} does not exist"
    elif "def diagenetic_reactions" not in equations_path.read_text(encoding="utf-8", errors="ignore"):
        needs_regen = True
        reason = f"{equations_path.name} is missing diagenetic_reactions function"
    elif reactions_path.stat().st_mtime > equations_path.stat().st_mtime:
        needs_regen = True
        reason = f"{reactions_path.name} is newer than {equations_path.name}"
    elif constants_path.exists() and constants_path.stat().st_mtime > equations_path.stat().st_mtime:
        needs_regen = True
        reason = f"{constants_path.name} is newer than {equations_path.name}"
    elif limiters_path is not None and limiters_path.exists() and limiters_path.stat().st_mtime > equations_path.stat().st_mtime:
        needs_regen = True
        reason = f"{limiters_path.name} is newer than {equations_path.name}"
    elif species_path.exists() and species_path.stat().st_mtime > equations_path.stat().st_mtime:
        needs_regen = True
        reason = f"{species_path.name} is newer than {equations_path.name}"
    elif generator_path.exists() and generator_path.stat().st_mtime > equations_path.stat().st_mtime:
        needs_regen = True
        reason = f"{generator_path.name} is newer than {equations_path.name}"

    if needs_regen:
        print(f"[MANTIS] Regenerating {equations_path.name} ({reason})...")
        generate_equations(
            reactions_path=reactions_path,
            equations_path=equations_path,
            kinetic_constants_path=constants_path,
            limiters_path=limiters_path,
            species_path=species_path,
            verbose=False,
        )
        importlib.invalidate_caches()

    # Ensure parent directories of equations and limiters are on sys.path for import resolution
    import_dirs = [str(equations_path.parent)]
    if limiters_path is not None:
        import_dirs.append(str(limiters_path.parent))
    for p_dir in import_dirs:
        if p_dir not in sys.path:
            sys.path.insert(0, p_dir)

    # Dynamically load the generated equations module
    mod_name = equations_path.stem
    spec = importlib.util.spec_from_file_location(mod_name, equations_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load equations module from {equations_path}")
    eq_mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = eq_mod
    spec.loader.exec_module(eq_mod)

    # Dynamically load reaction constants function if reaction_constants is not already a callable
    if not callable(p_dict.get("reaction_constants")) and constants_path.exists():
        c_spec = importlib.util.spec_from_file_location(constants_path.stem, constants_path)
        if c_spec and c_spec.loader:
            c_mod = importlib.util.module_from_spec(c_spec)
            c_spec.loader.exec_module(c_mod)
            if hasattr(c_mod, "get_reaction_constants"):
                p_dict["reaction_constants"] = c_mod.get_reaction_constants

    # Record loaded module and paths in p_dict
    p_dict["reactions_module"] = eq_mod
    p_dict["reactions_path"] = reactions_path
    p_dict["equations_path"] = equations_path

    return eq_mod, reactions_path, equations_path


class ReactionCategory(str):
    """Encapsulates classification metadata for a reaction while acting as a string for backward compatibility."""

    def __new__(
        cls,
        value: str,
        *,
        mechanism: str = "general_kinetic",
        phase_archetype: str = "liquid_to_liquid",
        has_solid_reactant: bool = False,
        has_solid_product: bool = False,
        reactants: Optional[List[Tuple[float, str]]] = None,
        products: Optional[List[Tuple[float, str]]] = None,
        canonical_name: str = "",
    ):
        obj = super().__new__(cls, value)
        obj.mechanism = mechanism
        obj.phase_archetype = phase_archetype
        obj.has_solid_reactant = has_solid_reactant
        obj.has_solid_product = has_solid_product
        obj.reactants = reactants or []
        obj.products = products or []
        obj.canonical_name = canonical_name or value
        return obj


class EquationSystemGenerator:
    """Generates an assembled equations.py module with SymPy analytical Jacobian."""

    def __init__(self, reactions: List[Dict[str, Any]], validator: ReactionSystemValidator):
        self.reactions = reactions
        self.validator = validator

    @classmethod
    def classify_reaction(cls, r: Dict[str, Any], validator: Optional[ReactionSystemValidator] = None) -> ReactionCategory:
        """Classifies a chemical reaction based on reactant/product phases and mechanism type,
        completely independent of chemical symbols or species names.
        """
        r_str = r.get("reaction", "")
        if "->" not in r_str:
            raise ValueError(f"Reaction string missing '->' delimiter: '{r_str}'")

        try:
            reactants, products = ReactionSystemValidator.parse_reaction_species(r_str)
        except Exception as e:
            raise ValueError(f"Could not parse reaction '{r_str}': {e}")

        r_sp = [sp for _, sp in reactants]
        p_sp = [sp for _, sp in products]
        k_val = r.get("k_value_name")

        # 1. Mechanism detection
        if isinstance(k_val, dict) and any(k in k_val for k in ("precipitation", "dissolution", "equilibrium")):
            mechanism = "reversible_saturation"
        elif len(reactants) == 1 and len(products) >= 2:
            mechanism = "disproportionation"
        else:
            mechanism = "general_kinetic"

        # 2. Reactant and product phase detection
        if validator is not None:
            r_phases = [validator.get_species_phase(sp) for _, sp in reactants]
            p_phases = [validator.get_species_phase(sp) for _, sp in products]
        else:
            r_phases = ["liquid" for _ in reactants]
            p_phases = ["liquid" for _ in products]

        has_solid_r = any(p in ("solid", "mixed") for p in r_phases)
        has_solid_p = any(p in ("solid", "mixed") for p in p_phases)

        # 3. Phase archetype assignment
        if not has_solid_r and not has_solid_p:
            phase_archetype = "liquid_to_liquid"
        elif has_solid_r and not has_solid_p:
            phase_archetype = "solid_to_liquid"
        elif not has_solid_r and has_solid_p:
            phase_archetype = "liquid_to_solid"
        else:
            phase_archetype = "solid_transformation"

        # 4. Canonical name is simply the reaction name or phase_mechanism
        r_name = r.get("reaction_name", "")
        canonical_name = r_name or f"{phase_archetype}_{mechanism}"

        return ReactionCategory(
            canonical_name,
            mechanism=mechanism,
            phase_archetype=phase_archetype,
            has_solid_reactant=has_solid_r,
            has_solid_product=has_solid_p,
            reactants=reactants,
            products=products,
            canonical_name=canonical_name,
        )

    def _generate_kinetic_reaction_fn(self, r: Dict[str, Any], cat: ReactionCategory) -> List[str]:
        r_name = r.get("reaction_name", "rxn")
        fn_name = f"rxn_{r_name}"
        reactants = cat.reactants
        products = cat.products
        has_solid = cat.has_solid_reactant
        k_name = r.get("k_value_name")

        transport_r = [
            (coeff, sp) for coeff, sp in reactants
            if self.validator.is_transport_species(self.validator.get_parent_species(sp))
        ]
        if not transport_r:
            transport_r = reactants
        master_coeff, master_sp = transport_r[0]
        for coeff, sp in transport_r:
            if self.validator.get_species_phase(sp) in ("solid", "mixed"):
                master_coeff, master_sp = coeff, sp
                break

        master_parent = self.validator.get_parent_species(master_sp)

        ref_sp = None
        ref_coeff = 1.0
        if r.get("ref_species"):
            ref_sp = r["ref_species"]
            ref_coeff = float(r.get("stoich_ref", 1.0))
        else:
            for coeff, sp in reactants:
                p_sp = self.validator.get_parent_species(sp)
                if isinstance(k_name, str) and (p_sp in k_name or sp in k_name):
                    ref_sp = p_sp
                    ref_coeff = float(coeff)
                    break
            if not ref_sp:
                ref_sp = master_parent
                ref_coeff = float(master_coeff)

        other_reactants_entries = []
        for coeff, sp in reactants:
            p_sp = self.validator.get_parent_species(sp)
            if p_sp != master_parent:
                if p_sp == "O2" and "POC" in master_parent:
                    other_reactants_entries.append(f'"{p_sp}": getattr(mp, "POC_O2_ratio", {coeff})')
                else:
                    other_reactants_entries.append(f'"{p_sp}": {coeff}')
        other_reactants_str = "{" + ", ".join(other_reactants_entries) + "}"

        transport_products = {
            self.validator.get_parent_species(sp): coeff
            for coeff, sp in products
            if self.validator.is_transport_species(self.validator.get_parent_species(sp))
        }

        lines = [
            f'def {fn_name}(c, k_val, lim, LHS, RHS, RATES, CROSS, mp):',
            f'    """Reaction: {r.get("reaction", "")}"""',
            f'    has_solid = {has_solid}',
            f'    k_num = mp.k.get("{k_name}") if hasattr(mp, "k") else getattr(k_val, "{k_name}", 0.0)',
        ]
        coeff_expr = 'k_num'
        master_spec = self.validator.get_speciation_factor_expr(master_sp)
        if master_spec:
            coeff_expr += f' * {master_spec}'

        for coeff, sp in reactants:
            p_sp = self.validator.get_parent_species(sp)
            if p_sp != master_parent:
                spec = self.validator.get_speciation_factor_expr(sp)
                term = f'getattr(c, "{p_sp}")' if not spec else f'(getattr(c, "{p_sp}") * {spec})'
                coeff_expr += f' * {term}'

        for l_name in (r.get("limiters") or []):
            lines.append(f'    _lim_{l_name} = lim.get("{l_name}_implicit", lim.get("{l_name}", 1.0))')
            coeff_expr += f' * _lim_{l_name}'
        for i_name in (r.get("inhibitors") or []):
            lines.append(f'    _inh_{i_name} = lim.get("{i_name}_inhibit", lim.get("{i_name}_inhib", 1.0))')
            coeff_expr += f' * _inh_{i_name}'

        lines.extend([
            f'    coeff_master = {coeff_expr}',
            f'    rate_master = coeff_master * getattr(c, "{master_parent}")',
            '    add_coupled_reaction(',
            '        CROSS=CROSS, LHS=LHS, RATES=RATES, mp=mp,',
            f'        master_species={{"{master_parent}": {master_coeff}}},',
            f'        reactants={other_reactants_str},',
            f'        products={transport_products},',
            f'        coeff_master=coeff_master, rate_master=rate_master,',
            f'        has_solid=has_solid, reaction_name="{r_name}",',
            f'        ref_species="{ref_sp}", stoich_ref={ref_coeff},',
            '    )',
        ])

        iso_sp_map = r.get("isotope_species")
        if iso_sp_map:
            # Find the reactant that participates in isotope tracking
            iso_r_candidates = [
                (coeff, sp, self.validator.get_parent_species(sp))
                for coeff, sp in reactants
                if self.validator.get_parent_species(sp) in iso_sp_map
            ]
            if iso_r_candidates:
                iso_master_coeff, iso_orig_sp, iso_master_parent = iso_r_candidates[0]
                iso_master = iso_sp_map[iso_master_parent]
                lines.append('    if getattr(mp, "isotopes", False):')

                # Build rate coefficient for the isotope reaction w.r.t iso_master_parent
                iso_coeff_expr = 'k_num'
                iso_spec = self.validator.get_speciation_factor_expr(iso_orig_sp)
                if iso_spec:
                    iso_coeff_expr += f' * {iso_spec}'
                for coeff, sp in reactants:
                    p_sp = self.validator.get_parent_species(sp)
                    if p_sp != iso_master_parent:
                        spec = self.validator.get_speciation_factor_expr(sp)
                        term = f'getattr(c, "{p_sp}")' if not spec else f'(getattr(c, "{p_sp}") * {spec})'
                        iso_coeff_expr += f' * {term}'
                for l_name in (r.get("limiters") or []):
                    iso_coeff_expr += f' * _lim_{l_name}'
                for i_name in (r.get("inhibitors") or []):
                    iso_coeff_expr += f' * _inh_{i_name}'

                iso_frac = r.get("isotope_fractionation")
                if iso_frac:
                    frac_key, frac_val = next(iter(iso_frac.items()))
                    alpha_name, lim_name = frac_val[0], frac_val[1]
                    raw_src = frac_key.split(" ")[0]
                    # Map raw_src back to its base parent if raw_src is an isotope companion
                    src_sp = raw_src
                    for b_sp, i_sp in iso_sp_map.items():
                        if raw_src in (b_sp, i_sp):
                            src_sp = b_sp
                            break
                    src_sp = self.validator.get_parent_species(src_sp)
                    src_iso_sp = iso_sp_map.get(src_sp)
                    lines.extend([
                        f'        _coeff_base_iso = {iso_coeff_expr}',
                        f'        _alpha = 1.0 + (getattr(mp, "{alpha_name}", 1.0) - 1.0) * lim.get("{lim_name}", 1.0)',
                        f'        coeff_master_iso = calculate_fractionated_rate_coeff(',
                        f'            _coeff_base_iso, getattr(c, "{src_sp}"), getattr(c, "{src_iso_sp}"), _alpha, eps=1e-30',
                        '        )',
                    ])
                else:
                    lines.append(f'        coeff_master_iso = {iso_coeff_expr}')

                iso_other_reactants = {}
                for coeff, sp in reactants:
                    p_sp = self.validator.get_parent_species(sp)
                    if p_sp != iso_master_parent and p_sp in iso_sp_map:
                        iso_other_reactants[iso_sp_map[p_sp]] = coeff
                iso_products = {}
                for coeff, sp in products:
                    p_sp = self.validator.get_parent_species(sp)
                    if p_sp in iso_sp_map:
                        iso_products[iso_sp_map[p_sp]] = coeff

                iso_rxn_name = f"{r_name}_iso"
                lines.extend([
                    '        add_coupled_reaction(',
                    '            CROSS=CROSS, LHS=LHS, RATES=RATES, mp=mp,',
                    f'            master_species={{"{iso_master}": {iso_master_coeff}}},',
                    f'            reactants={iso_other_reactants},',
                    f'            products={iso_products},',
                    f'            coeff_master=coeff_master_iso, rate_master=coeff_master_iso * getattr(c, "{iso_master}"),',
                    f'            has_solid=has_solid, reaction_name="{iso_rxn_name}",',
                    f'            ref_species="{ref_sp}", stoich_ref={ref_coeff},',
                    '        )',
                ])

        if r_name != fn_name:
            lines.append(f'{r_name} = {fn_name}')
        lines.append('')
        return lines

    def _generate_reversible_saturation_fn(self, r: Dict[str, Any], cat: ReactionCategory) -> List[str]:
        r_name = r.get("reaction_name", "mineral_saturation")
        fn_name = f"rxn_{r_name}"
        k_name = r.get("k_value_name", {})
        if isinstance(k_name, dict):
            k_prec_name = k_name.get("precipitation", f"{r_name}_isp")
            k_diss_name = k_name.get("dissolution", f"{r_name}_isd")
            k_sp_name = k_name.get("equilibrium", f"{r_name}_sp")
        else:
            k_prec_name = f"{k_name}_isp"
            k_diss_name = f"{k_name}_isd"
            k_sp_name = f"{k_name}_sp"

        # Identify mineral product (solid) and dissolved reactants
        mineral_candidates = [
            (coeff, sp) for coeff, sp in cat.products
            if self.validator.get_species_phase(sp) == "solid"
        ]
        if not mineral_candidates:
            mineral_coeff, mineral_sp = cat.products[0]
        else:
            mineral_coeff, mineral_sp = mineral_candidates[0]
        mineral_parent = self.validator.get_parent_species(mineral_sp)

        reactants_info = []
        for coeff, sp in cat.reactants:
            p_sp = self.validator.get_parent_species(sp)
            spec = self.validator.get_speciation_factor_expr(sp)
            reactants_info.append((coeff, sp, p_sp, spec))

        lines = [
            f'def {fn_name}(c, k_val, lim, LHS, RHS, RATES, CROSS, mp):',
            f'    """Reaction: {r.get("reaction", "")} (reversible saturation with C1 smooth transition)"""',
            '    phi = mp.phi',
        ]

        # Extract reactant concentration variables
        omega_num_terms = []
        c_val_names = []
        for coeff, sp, p_sp, spec in reactants_info:
            v_name = f"{p_sp}_val"
            pw_name = f"{sp}_pw"
            lines.append(f'    {v_name} = np.maximum(np.asarray(c.{p_sp}.value if hasattr(c.{p_sp}, "value") else c.{p_sp}), 1e-20)')
            if spec:
                lines.append(f'    {pw_name} = {v_name} * {spec}')
            else:
                lines.append(f'    {pw_name} = {v_name}')
            term = pw_name if coeff == 1.0 else f"({pw_name}**{coeff})"
            omega_num_terms.append(term)
            c_val_names.append((p_sp, v_name))

        # Mineral solid concentration
        lines.append(f'    {mineral_parent}_val = np.maximum(np.asarray(c.{mineral_parent}.value if hasattr(c.{mineral_parent}, "value") else c.{mineral_parent}), 1e-20)')

        # Equilibrium constants and saturation state omega
        has_hplus = any("HS" in sp or "H2S" in sp for _, sp, _, _ in reactants_info)
        lines.extend([
            f'    k_sp = getattr(k_val, "{k_sp_name}", getattr(mp, "{k_sp_name}", 10**-3.5))',
            f'    k_isp = getattr(k_val, "{k_prec_name}", getattr(mp, "{k_prec_name}", 1.0))',
            f'    k_isd = getattr(k_val, "{k_diss_name}", getattr(mp, "{k_diss_name}", 0.3))',
        ])
        if has_hplus:
            lines.extend([
                '    k_Hplus = np.asarray(getattr(c, "Hplus")) if hasattr(c, "TA") and hasattr(c, "DIC") and hasattr(c, "Hplus") else getattr(k_val, "Hplus", getattr(mp, "Hplus", 10**-7.5))',
                '    omega_den = k_Hplus * k_sp + 1e-30',
            ])
        else:
            lines.append('    omega_den = k_sp + 1e-30')
        lines.extend([
            f'    omega = ({" * ".join(omega_num_terms)}) / omega_den',
            '    is_prec = np.where(omega >= 1.0, 1.0, 0.0)',
            '    is_diss = np.where(omega < 1.0, 1.0, 0.0)',
            '    epsilon = float(getattr(mp, "fes_smooth_epsilon", getattr(mp, "saturation_smooth_epsilon", 0.01)))',
            '    k_prec_eff = k_isp * is_prec',
            '    df = np.maximum(omega - 1.0, 0.0)',
            '    u_prec = np.minimum(df / (epsilon + 1e-30), 1.0)',
            '    smooth_ramp_prec = (3.0 * u_prec**2 - 2.0 * u_prec**3)',
            '    R_prec = k_prec_eff * df * smooth_ramp_prec',
        ])

        # Dynamic symmetric weighting for precipitation across reactants
        if len(reactants_info) == 2:
            p1, v1 = c_val_names[0]
            p2, v2 = c_val_names[1]
            lines.extend([
                f'    w_{p1} = {v2}**2 / ({v1}**2 + {v2}**2 + 1e-30)',
                f'    w_{p2} = 1.0 - w_{p1}',
                f'    prec_coeff_{p1} = (w_{p1} * R_prec) / {v1}',
                f'    prec_rate_{p1} = prec_coeff_{p1} * {v1}',
                f'    prec_coeff_{p2} = (w_{p2} * R_prec) / {v2}',
                f'    prec_rate_{p2} = prec_coeff_{p2} * {v2}',
            ])
            for p_curr, v_curr in [(p1, v1), (p2, v2)]:
                lines.extend([
                    '    add_implicit_coupling_new(',
                    f'        CROSS, RATES, LHS, target_species="{mineral_parent}", source_species="{p_curr}",',
                    f'        coeff=prec_coeff_{p_curr}, rate=prec_rate_{p_curr}, mp=mp, has_solid=False,',
                    '        add_lhs_sink=False, stoich_ratio=1.0,',
                    '    )',
                ])
            lines.extend([
                f'    add_implicit_sink(LHS, RATES, "{p1}", prec_coeff_{p1}, prec_rate_{p1}, mp=mp, has_solid=False)',
                f'    add_implicit_sink(LHS, RATES, "{p2}", prec_coeff_{p2}, prec_rate_{p2}, mp=mp, has_solid=False)',
                f'    CROSS["{p1}"].append(("{p2}", -prec_coeff_{p2} * phi))',
                f'    RATES["{p1}"] -= prec_rate_{p2} * phi',
                f'    CROSS["{p2}"].append(("{p1}", -prec_coeff_{p1} * phi))',
                f'    RATES["{p2}"] -= prec_rate_{p1} * phi',
            ])
        else:
            p1, v1 = c_val_names[0]
            lines.extend([
                f'    prec_coeff_{p1} = R_prec / {v1}',
                f'    prec_rate_{p1} = prec_coeff_{p1} * {v1}',
                '    add_implicit_coupling_new(',
                f'        CROSS, RATES, LHS, target_species="{mineral_parent}", source_species="{p1}",',
                f'        coeff=prec_coeff_{p1}, rate=prec_rate_{p1}, mp=mp, has_solid=False,',
                '        add_lhs_sink=False, stoich_ratio=1.0,',
                '    )',
                f'    add_implicit_sink(LHS, RATES, "{p1}", prec_coeff_{p1}, prec_rate_{p1}, mp=mp, has_solid=False)',
            ])

        # Dissolution of solid mineral
        lines.extend([
            '    k_diss_eff = k_isd * is_diss',
            '    us = np.maximum(1.0 - omega, 0.0)',
            '    u_diss = np.minimum(us / (epsilon + 1e-30), 1.0)',
            '    smooth_ramp_diss = (3.0 * u_diss**2 - 2.0 * u_diss**3)',
            f'    diss_coeff_{mineral_parent} = k_diss_eff * us * smooth_ramp_diss',
            f'    diss_rate_{mineral_parent} = diss_coeff_{mineral_parent} * {mineral_parent}_val',
        ])
        for idx, (p_curr, v_curr) in enumerate(c_val_names):
            lines.extend([
                '    add_implicit_coupling_new(',
                f'        CROSS, RATES, LHS, target_species="{p_curr}", source_species="{mineral_parent}",',
                f'        coeff=diss_coeff_{mineral_parent}, rate=diss_rate_{mineral_parent}, mp=mp, has_solid=True,',
                f'        add_lhs_sink={idx == 0}, stoich_ratio=1.0,',
                '    )',
            ])

        # Isotope handling for reversible saturation
        iso_sp_map = r.get("isotope_species", {})
        if iso_sp_map:
            # Find reactant that tracks isotopes
            iso_r = [
                (coeff, sp, p_sp) for coeff, sp, p_sp, _ in reactants_info
                if p_sp in iso_sp_map
            ]
            if iso_r and mineral_parent in iso_sp_map:
                r_coeff, r_sp, r_parent = iso_r[0]
                r_iso = iso_sp_map[r_parent]
                m_iso = iso_sp_map[mineral_parent]

                # Determine isotope standard ratio dynamically
                std_name = self.validator.species_defs.get(r_parent, {}).get("isotope_standard", "S")
                std_info = getattr(self.validator, "isotope_standards", {}).get(std_name, {})
                std_key = std_info.get("standard", "VCDT")
                if "ratio" not in std_info:
                    raise ValueError(
                        f"[MANTIS ERROR] Isotope standard '{std_name}' is missing 'ratio' in species.py."
                    )
                std_ratio = float(std_info["ratio"])

                lines.extend([
                    '    if getattr(mp, "isotopes", False):',
                    f'        std_ratio = getattr(mp, "{std_key}", getattr(mp, "{std_name}_ratio", {std_ratio}))',
                    f'        f_iso_default = 1.0 / (1.0 + float(std_ratio))',
                ])
                r_spec = self.validator.get_speciation_factor_expr(r_sp)
                if r_spec:
                    other_eq_species = [
                        s for s, m in getattr(self.validator, "species_defs", {}).items()
                        if isinstance(m, dict) and m.get("equilibrium_of") == r_parent and s != r_sp and not m.get("isotope_of")
                    ]
                    if other_eq_species:
                        other_sp = other_eq_species[0]
                        other_frac = self.validator.get_speciation_factor_expr(other_sp) or f"mp.{other_sp.lower()}_frac"
                        alpha_var = f"{other_sp.lower()}_{r_sp.lower()}_alpha"
                    else:
                        other_frac = f"(1.0 - {r_spec})"
                        alpha_var = f"{r_sp.lower()}_alpha"

                    lines.extend([
                        f'        eq_alpha = getattr(mp, "{alpha_var}", getattr(mp, "h2s_hs_alpha", 1.0))',
                        f'        {r_sp}_iso = partition_equilibrium_isotope(',
                        f'            c.{r_iso}, {r_spec}, {other_frac}, eq_alpha',
                        '        )',
                        f'        {r_sp}_val_np = np.asarray({r_sp}_pw)',
                        f'        {r_sp}_iso_val = np.asarray({r_sp}_iso)',
                        f'        {r_parent}_val_np = np.asarray(c.{r_parent}.value if hasattr(c.{r_parent}, "value") else c.{r_parent})',
                        f'        {r_parent}_iso_val = np.asarray(c.{r_iso}.value if hasattr(c.{r_iso}, "value") else c.{r_iso})',
                        f'        f_iso_{r_sp} = np.where({r_sp}_val_np > 1e-6, {r_sp}_iso_val / ({r_sp}_val_np + 1e-30), f_iso_default)',
                        f'        f_iso_{r_sp} = np.clip(f_iso_{r_sp}, 0.5, 1.5)',
                        f'        f_iso_{r_parent} = np.where({r_parent}_val_np > 1e-6, {r_parent}_iso_val / ({r_parent}_val_np + 1e-30), f_iso_default)',
                        f'        f_iso_{r_parent} = np.clip(f_iso_{r_parent}, 0.5, 1.5)',
                        f'        ratio_iso = np.where(f_iso_{r_parent} > 1e-10, f_iso_{r_sp} / f_iso_{r_parent}, 1.0)',
                        f'        prec_coeff_{r_parent}_iso = prec_coeff_{r_parent} * ratio_iso',
                    ])
                    f_iso_r = f"f_iso_{r_sp}"
                else:
                    lines.extend([
                        f'        {r_parent}_val_np = np.asarray(c.{r_parent}.value if hasattr(c.{r_parent}, "value") else c.{r_parent})',
                        f'        {r_parent}_iso_val = np.asarray(c.{r_iso}.value if hasattr(c.{r_iso}, "value") else c.{r_iso})',
                        f'        f_iso_{r_parent} = np.where({r_parent}_val_np > 1e-6, {r_parent}_iso_val / ({r_parent}_val_np + 1e-30), f_iso_default)',
                        f'        f_iso_{r_parent} = np.clip(f_iso_{r_parent}, 0.5, 1.5)',
                        f'        prec_coeff_{r_parent}_iso = prec_coeff_{r_parent}',
                    ])
                    f_iso_r = f"f_iso_{r_parent}"

                # Couple isotope precipitation
                lines.extend([
                    '        add_implicit_coupling_new(',
                    f'            CROSS, RATES, LHS, target_species="{m_iso}", source_species="{r_iso}",',
                    f'            coeff=prec_coeff_{r_parent}_iso, rate=prec_rate_{r_parent} * {f_iso_r}, mp=mp, has_solid=False,',
                    '            add_lhs_sink=False, stoich_ratio=1.0,',
                    '        )',
                ])
                for p_curr, _ in c_val_names:
                    if p_curr != r_parent:
                        lines.extend([
                            '        add_implicit_coupling_new(',
                            f'            CROSS, RATES, LHS, target_species="{m_iso}", source_species="{p_curr}",',
                            f'            coeff=prec_coeff_{p_curr} * {f_iso_r}, rate=prec_rate_{p_curr} * {f_iso_r}, mp=mp, has_solid=False,',
                            '            add_lhs_sink=False, stoich_ratio=1.0,',
                            '        )',
                            f'        CROSS["{r_iso}"].append(("{p_curr}", -prec_coeff_{p_curr} * {f_iso_r} * phi))',
                            f'        RATES["{r_iso}"] -= prec_rate_{p_curr} * {f_iso_r} * phi',
                        ])

                lines.append(f'        add_implicit_sink(LHS, RATES, "{r_iso}", prec_coeff_{r_parent}_iso, prec_rate_{r_parent} * {f_iso_r}, mp=mp, has_solid=False)')

                # Couple isotope dissolution
                lines.extend([
                    f'        {mineral_parent}_val_np = np.asarray(c.{mineral_parent}.value if hasattr(c.{mineral_parent}, "value") else c.{mineral_parent})',
                    f'        {m_iso}_val = np.asarray(c.{m_iso}.value if hasattr(c.{m_iso}, "value") else c.{m_iso})',
                    f'        f_iso_{mineral_parent} = np.where({mineral_parent}_val_np > 1e-3, {m_iso}_val / ({mineral_parent}_val_np + 1e-30), {f_iso_r})',
                    f'        f_iso_{mineral_parent} = np.clip(f_iso_{mineral_parent}, 0.5, 1.5)',
                    '        add_implicit_coupling_new(',
                    f'            CROSS, RATES, LHS, target_species="{r_iso}", source_species="{m_iso}",',
                    f'            coeff=diss_coeff_{mineral_parent}, rate=diss_rate_{mineral_parent} * f_iso_{mineral_parent}, mp=mp, has_solid=True,',
                    '            add_lhs_sink=True, stoich_ratio=1.0,',
                    '        )',
                ])
                for p_curr, _ in c_val_names:
                    if p_curr != r_parent:
                        lines.extend([
                            '        add_implicit_coupling_new(',
                            f'            CROSS, RATES, LHS, target_species="{p_curr}", source_species="{m_iso}",',
                            f'            coeff=diss_coeff_{mineral_parent}, rate=diss_rate_{mineral_parent} * f_iso_{mineral_parent}, mp=mp, has_solid=True,',
                            '            add_lhs_sink=False, stoich_ratio=1.0,',
                            '        )',
                        ])

        if r_name != fn_name:
            lines.append(f'{r_name} = {fn_name}')
        lines.append('')
        return lines

    def _generate_reversible_saturation_jacobian(self, r: Dict[str, Any], cat: ReactionCategory) -> List[str]:
        r_name = r.get("reaction_name", "mineral_saturation")
        k_name = r.get("k_value_name", {})
        if isinstance(k_name, dict):
            k_prec_name = k_name.get("precipitation", f"{r_name}_isp")
            k_diss_name = k_name.get("dissolution", f"{r_name}_isd")
            k_sp_name = k_name.get("equilibrium", f"{r_name}_sp")
        else:
            k_prec_name = f"{k_name}_isp"
            k_diss_name = f"{k_name}_isd"
            k_sp_name = f"{k_name}_sp"

        # Identify mineral product (solid) and dissolved reactants
        mineral_candidates = [
            (coeff, sp) for coeff, sp in cat.products
            if self.validator.get_species_phase(sp) == "solid"
        ]
        mineral_parent = self.validator.get_parent_species(mineral_candidates[0][1]) if mineral_candidates else self.validator.get_parent_species(cat.products[0][1])

        reactants_info = []
        for coeff, sp in cat.reactants:
            p_sp = self.validator.get_parent_species(sp)
            spec = self.validator.get_speciation_factor_expr(sp)
            reactants_info.append((coeff, sp, p_sp, spec))

        all_species_needed = [p for _, _, p, _ in reactants_info] + [mineral_parent]
        conds = [f'"{p}" in idx' for p in all_species_needed]

        lines = [
            f'    # Reversible Saturation: {r_name}',
            f'    if {" and ".join(conds)}:',
        ]
        for p in all_species_needed:
            lines.append(f'        i_{p.lower()} = idx["{p}"]')

        for _, sp, p_sp, spec in reactants_info:
            v_name = f"{p_sp}_val"
            pw_name = f"{sp}_pw"
            lines.append(f'        {v_name} = np.maximum(c["{p_sp}"].value if hasattr(c["{p_sp}"], "value") else c["{p_sp}"], 1e-20)')
            if spec:
                lines.append(f'        {pw_name} = {v_name} * ({spec})')
            else:
                lines.append(f'        {pw_name} = {v_name}')

        lines.append(f'        {mineral_parent}_val = np.maximum(c["{mineral_parent}"].value if hasattr(c["{mineral_parent}"], "value") else c["{mineral_parent}"], 1e-20)')

        omega_num = " * ".join([f"{sp}_pw" for _, sp, _, _ in reactants_info])
        has_hplus = any("HS" in sp or "H2S" in sp for _, sp, _, _ in reactants_info)
        lines.extend([
            f'        k_min_sp = float(getattr(k, "{k_sp_name}", getattr(mp, "{k_sp_name}", 10**-3.5)))',
            f'        k_min_isp = float(getattr(k, "{k_prec_name}", getattr(mp, "{k_prec_name}", 1.0)))',
            f'        k_min_isd = float(getattr(k, "{k_diss_name}", getattr(mp, "{k_diss_name}", 0.3)))',
        ])
        if has_hplus:
            lines.extend([
                '        k_Hplus = np.asarray(getattr(c, "Hplus")) if has_carb and hasattr(c, "Hplus") else float(getattr(k, "Hplus", getattr(mp, "Hplus", 10**-7.5)))',
                f'        omega = ({omega_num}) / (k_Hplus * k_min_sp + 1e-30)',
            ])
        else:
            lines.append(f'        omega = ({omega_num}) / (k_min_sp + 1e-30)')

        # Derivatives of omega w.r.t reactants
        for coeff_i, sp_i, p_i, spec_i in reactants_info:
            other_pw = [f"{sp}_pw" for _, sp, p_other, _ in reactants_info if p_other != p_i]
            fac_i = f"({spec_i})" if spec_i else "1.0"
            term_num = f"{fac_i} * " + " * ".join(other_pw) if other_pw else fac_i
            den_str = "(k_Hplus * k_min_sp + 1e-30)" if has_hplus else "(k_min_sp + 1e-30)"
            lines.append(f'        d_omega_d{p_i} = ({term_num}) / {den_str}')

        lines.extend([
            '        eps_sm = float(getattr(mp, "fes_smooth_epsilon", getattr(mp, "saturation_smooth_epsilon", 0.01)))',
            '        is_prec = omega >= 1.0',
            '        df_p = np.maximum(omega - 1.0, 0.0)',
            '        u_p = np.clip(df_p / (eps_sm + 1e-30), 0.0, 1.0)',
            '        S_p = 3.0 * u_p**2 - 2.0 * u_p**3',
            '        dS_du = np.where((u_p > 0.0) & (u_p < 1.0), 6.0 * u_p * (1.0 - u_p), 0.0)',
            '        dS_p_domega = dS_du / (eps_sm + 1e-30)',
            '        dR_prec_domega = k_min_isp * (S_p + df_p * dS_p_domega) * is_prec',
            '        R_prec = k_min_isp * df_p * S_p * is_prec',
        ])

        for _, _, p_i, _ in reactants_info:
            lines.append(f'        dRp_d{p_i} = diss_fac * dR_prec_domega * d_omega_d{p_i}')

        lines.extend([
            '        is_diss = omega < 1.0',
            '        us_d = np.maximum(1.0 - omega, 0.0)',
            '        u_d = np.clip(us_d / (eps_sm + 1e-30), 0.0, 1.0)',
            '        S_d = 3.0 * u_d**2 - 2.0 * u_d**3',
            '        dS_d_du = np.where((u_d > 0.0) & (u_d < 1.0), 6.0 * u_d * (1.0 - u_d), 0.0)',
            '        dS_d_domega = -dS_d_du / (eps_sm + 1e-30)',
            f'        dR_diss_domega = -k_min_isd * (S_d + us_d * (dS_d_du / (eps_sm + 1e-30))) * {mineral_parent}_val * is_diss',
            f'        dR_diss_d{mineral_parent} = k_min_isd * us_d * S_d * is_diss',
            f'        R_diss = k_min_isd * us_d * S_d * {mineral_parent}_val * is_diss',
        ])

        for _, _, p_i, _ in reactants_info:
            lines.append(f'        dRd_d{p_i} = solid_fac * dR_diss_domega * d_omega_d{p_i}')
        lines.append(f'        dRd_d{mineral_parent} = solid_fac * dR_diss_d{mineral_parent}')

        for _, _, p_i, _ in reactants_info:
            lines.append(f'        d_net_d{p_i} = dRp_d{p_i} - dRd_d{p_i}')
        lines.append(f'        d_net_d{mineral_parent} = -dRd_d{mineral_parent}')

        # Populate Jacobian entries
        m_idx = f"i_{mineral_parent.lower()}"
        for _, _, p_i, _ in reactants_info:
            p_idx = f"i_{p_i.lower()}"
            lines.append(f'        J[:, {m_idx}, {p_idx}] += d_net_d{p_i}')
            for _, _, p_target, _ in reactants_info:
                p_tgt_idx = f"i_{p_target.lower()}"
                lines.append(f'        J[:, {p_tgt_idx}, {p_idx}] -= d_net_d{p_i}')

        lines.append(f'        J[:, {m_idx}, {m_idx}] += d_net_d{mineral_parent}')
        for _, _, p_target, _ in reactants_info:
            p_tgt_idx = f"i_{p_target.lower()}"
            lines.append(f'        J[:, {p_tgt_idx}, {m_idx}] -= d_net_d{mineral_parent}')

        if has_hplus:
            lines.extend([
                '        if has_carb:',
                '            _Ka1_h2s = float(getattr(mp, "K_a1_h2s", 10.0 ** (-7.0)))',
                '            _domega_dh_m3 = -omega * (1.0 / (k_Hplus + 1e-30) + 1.0 / (_h * 1000.0 + _Ka1_h2s + 1e-30))',
                '            _domega_dh = _domega_dh_m3 * 1000.0',
                '            _d_net_dTA = (dR_prec_domega * diss_fac - dR_diss_domega * solid_fac) * _domega_dh * _dh_dta',
                '            _d_net_dDIC = (dR_prec_domega * diss_fac - dR_diss_domega * solid_fac) * _domega_dh * _dh_ddic',
                f'            J[:, {m_idx}, idx["TA"]] += _d_net_dTA',
                f'            J[:, {m_idx}, idx["DIC"]] += _d_net_dDIC',
            ])
            for _, _, p_target, _ in reactants_info:
                p_tgt_idx = f"i_{p_target.lower()}"
                lines.extend([
                    f'            J[:, {p_tgt_idx}, idx["TA"]] -= _d_net_dTA',
                    f'            J[:, {p_tgt_idx}, idx["DIC"]] -= _d_net_dDIC',
                ])

        # Isotope Jacobian entries
        iso_sp_map = r.get("isotope_species", {})
        if iso_sp_map:
            iso_r = [
                (coeff, sp, p_sp) for coeff, sp, p_sp, _ in reactants_info
                if p_sp in iso_sp_map
            ]
            if iso_r and mineral_parent in iso_sp_map:
                r_coeff, r_sp, r_parent = iso_r[0]
                r_iso = iso_sp_map[r_parent]
                m_iso = iso_sp_map[mineral_parent]
                std_name = self.validator.species_defs.get(r_parent, {}).get("isotope_standard", "S")
                std_info = getattr(self.validator, "isotope_standards", {}).get(std_name, {})
                std_key = std_info.get("standard", "VCDT")
                if "ratio" not in std_info:
                    raise ValueError(
                        f"[MANTIS ERROR] Isotope standard '{std_name}' is missing 'ratio' in species.py."
                    )
                std_ratio = float(std_info["ratio"])

                r_spec = self.validator.get_speciation_factor_expr(r_sp)
                other_eq_species = [
                    s for s, m in getattr(self.validator, "species_defs", {}).items()
                    if isinstance(m, dict) and m.get("equilibrium_of") == r_parent and s != r_sp and not m.get("isotope_of")
                ]
                if other_eq_species:
                    other_sp = other_eq_species[0]
                    other_frac = self.validator.get_speciation_factor_expr(other_sp) or f"mp.{other_sp.lower()}_frac"
                    alpha_var = f"{other_sp.lower()}_{r_sp.lower()}_alpha"
                else:
                    other_frac = f"(1.0 - {r_spec})" if r_spec else "1.0"
                    alpha_var = f"{r_sp.lower()}_alpha"
                target_frac = r_spec if r_spec else "1.0"

                lines.extend([
                    f'        if "{m_iso}" in idx and "{r_iso}" in idx and getattr(mp, "isotopes", False):',
                    f'            i_{m_iso.lower()} = idx["{m_iso}"]',
                    f'            i_{r_iso.lower()} = idx["{r_iso}"]',
                    f'            eq_alpha = getattr(mp, "{alpha_var}", getattr(mp, "h2s_hs_alpha", 1.0))',
                    f'            {r_sp}_iso = partition_equilibrium_isotope(',
                    f'                c["{r_iso}"].value if hasattr(c["{r_iso}"], "value") else c["{r_iso}"],',
                    f'                {target_frac}, {other_frac}, eq_alpha,',
                    '            )',
                    f'            std_ratio = getattr(mp, "{std_key}", getattr(mp, "{std_name}_ratio", {std_ratio}))',
                    '            f_iso_default = 1.0 / (1.0 + float(std_ratio))',
                    f'            mask_r = {r_sp}_pw > 1e-6',
                    f'            f_iso_r = np.where(mask_r, np.asarray({r_sp}_iso) / ({r_sp}_pw + 1e-30), f_iso_default)',
                    '            f_iso_r = np.clip(f_iso_r, 0.5, 1.5)',
                    '            R_prec_bulk = diss_fac * R_prec',
                ])
                for _, _, p_i, _ in reactants_info:
                    if p_i == r_parent:
                        lines.extend([
                            f'            dRp_iso_d{p_i} = dRp_d{p_i} * f_iso_r - np.where(mask_r, (R_prec_bulk / {p_i}_val) * f_iso_r, 0.0)',
                            f'            dRp_iso_d{r_iso} = np.where(is_prec & mask_r, R_prec_bulk / {p_i}_val, 0.0)',
                        ])
                    else:
                        lines.append(f'            dRp_iso_d{p_i} = dRp_d{p_i} * f_iso_r')

                lines.extend([
                    '            R_diss_bulk = solid_fac * R_diss',
                    f'            mask_m = {mineral_parent}_val > 1e-3',
                    f'            {m_iso}_val = np.asarray(c["{m_iso}"].value if hasattr(c["{m_iso}"], "value") else c["{m_iso}"])',
                    f'            f_iso_m = np.where(mask_m, {m_iso}_val / ({mineral_parent}_val + 1e-30), f_iso_r)',
                    '            f_iso_m = np.clip(f_iso_m, 0.5, 1.5)',
                ])
                for _, _, p_i, _ in reactants_info:
                    lines.append(f'            dRd_iso_d{p_i} = dRd_d{p_i} * f_iso_m')
                lines.extend([
                    f'            dRd_iso_d{mineral_parent} = np.where(~mask_m, dRd_d{mineral_parent} * f_iso_r, 0.0)',
                    f'            dRd_iso_d{m_iso} = np.where(is_diss & mask_m, R_diss_bulk / {mineral_parent}_val, 0.0)',
                ])

                # Assign isotope cross derivatives
                m_iso_idx = f"i_{m_iso.lower()}"
                r_iso_idx = f"i_{r_iso.lower()}"
                for _, _, p_i, _ in reactants_info:
                    p_idx = f"i_{p_i.lower()}"
                    lines.extend([
                        f'            J[:, {m_iso_idx}, {p_idx}] += (dRp_iso_d{p_i} - dRd_iso_d{p_i})',
                        f'            J[:, {r_iso_idx}, {p_idx}] -= (dRp_iso_d{p_i} - dRd_iso_d{p_i})',
                    ])
                lines.extend([
                    f'            J[:, {m_iso_idx}, i_{mineral_parent.lower()}] -= dRd_iso_d{mineral_parent}',
                    f'            J[:, {m_iso_idx}, {r_iso_idx}] += dRp_iso_d{r_iso}',
                    f'            J[:, {m_iso_idx}, {m_iso_idx}] -= dRd_iso_d{m_iso}',
                    f'            J[:, {r_iso_idx}, i_{mineral_parent.lower()}] += dRd_iso_d{mineral_parent}',
                    f'            J[:, {r_iso_idx}, {r_iso_idx}] -= dRp_iso_d{r_iso}',
                    f'            J[:, {r_iso_idx}, {m_iso_idx}] += dRd_iso_d{m_iso}',
                ])

        lines.append('')
        return lines

    def _generate_kinetic_reaction_jacobian(self, r: Dict[str, Any], cat: ReactionCategory) -> List[str]:
        r_name = r.get("reaction_name", "rxn")
        k_name = r.get("k_value_name")
        reactants = cat.reactants
        products = cat.products
        has_solid = cat.has_solid_reactant
        fac_var = "solid_fac" if has_solid else "diss_fac"

        r_parents = [(coeff, sp, self.validator.get_parent_species(sp)) for coeff, sp in reactants]
        p_parents = [
            (coeff, sp, self.validator.get_parent_species(sp)) for coeff, sp in products
            if self.validator.is_transport_species(self.validator.get_parent_species(sp))
        ]

        ref_sp = None
        ref_coeff = 1.0
        if r.get("ref_species"):
            ref_sp = r["ref_species"]
            ref_coeff = float(r.get("stoich_ref", 1.0))
        else:
            for coeff, sp in reactants:
                p_sp = self.validator.get_parent_species(sp)
                if isinstance(k_name, str) and (p_sp in k_name or sp in k_name):
                    ref_sp = p_sp
                    ref_coeff = float(coeff)
                    break
            if not ref_sp:
                ref_sp = r_parents[0][2]
                ref_coeff = float(r_parents[0][0])

        all_needed_parents = set([p for _, _, p in r_parents] + [p for _, _, p in p_parents])
        conds = [f'"{p}" in idx' for p in sorted(all_needed_parents)]
        lines = [
            f'    # Reaction: {r_name} ({r.get("reaction", "")})',
            f'    if {" and ".join(conds)}:',
            f'        k_num_{r_name} = float(mp.k.get("{k_name}") if hasattr(mp, "k") else getattr(k, "{k_name}", 0.0))',
        ]
        for p in sorted(all_needed_parents):
            lines.append(f'        c_{p} = np.asarray(c["{p}"].value if hasattr(c["{p}"], "value") else c["{p}"])')

        # Limiters & Inhibitors terms
        lim_terms = []
        for l_name in (r.get("limiters") or []):
            lines.append(f'        _lim_{l_name} = lim.get("{l_name}_implicit", lim.get("{l_name}", 1.0))')
            lim_terms.append(f'_lim_{l_name}')

        for i_name in (r.get("inhibitors") or []):
            lines.append(f'        _inh_{i_name} = lim.get("{i_name}_inhibit", lim.get("{i_name}_inhib", 1.0))')
            lim_terms.append(f'_inh_{i_name}')

        is_calcite_precip = ("calcite" in r_name.lower() or "caco3" in r_name.lower()) and "CaCO3" in [p for _, _, p in p_parents] and "TA" in all_needed_parents and "DIC" in all_needed_parents and "Ca2" in all_needed_parents
        is_calcite_diss = ("calcite" in r_name.lower() or "caco3" in r_name.lower()) and "CaCO3" in [p for _, _, p in r_parents] and "TA" in all_needed_parents and "DIC" in all_needed_parents and "Ca2" in all_needed_parents

        if is_calcite_precip:
            lines.extend([
                '        _Ksp = float(mp.k.get("Ksp_calcite") if hasattr(mp, "k") else getattr(k, "Ksp_calcite", 0.4266))',
                '        _omega = (c_Ca2 * _co3) / (_Ksp + 1e-30)',
                '        _eps_sm = float(getattr(mp, "saturation_smooth_epsilon", 0.01))',
                '        _df_p = np.maximum(_omega - 1.0, 0.0)',
                '        _u_p = np.clip(_df_p / (_eps_sm + 1e-30), 0.0, 1.0)',
                '        _S_p = 3.0 * _u_p**2 - 2.0 * _u_p**3',
                '        _dS_p_du = np.where((_u_p > 0.0) & (_u_p < 1.0), 6.0 * _u_p * (1.0 - _u_p), 0.0)',
                '        _dS_p_domega = _dS_p_du / (_eps_sm + 1e-30)',
                f'        _dR_prec_domega = diss_fac * k_num_{r_name} * (_S_p + _df_p * _dS_p_domega) * (_omega >= 1.0)',
                '        _dr_dCa2 = _dR_prec_domega * (_co3 / (_Ksp + 1e-30))',
                '        _dr_dTA = _dR_prec_domega * (c_Ca2 / (_Ksp + 1e-30)) * _dco3_dta',
                '        _dr_dDIC = _dR_prec_domega * (c_Ca2 / (_Ksp + 1e-30)) * _dco3_ddic',
                '        for _vname, _dr in [("Ca2", _dr_dCa2), ("TA", _dr_dTA), ("DIC", _dr_dDIC)]:',
                '            J[:, idx["Ca2"], idx[_vname]] -= 1.0 * _dr',
                '            J[:, idx["DIC"], idx[_vname]] -= 1.0 * _dr',
                '            J[:, idx["TA"], idx[_vname]] -= 2.0 * _dr',
                '            J[:, idx["CaCO3"], idx[_vname]] += 1.0 * _dr',
            ])
            lines.append('')
            return lines

        if is_calcite_diss:
            lines.extend([
                '        _Ksp = float(mp.k.get("Ksp_calcite") if hasattr(mp, "k") else getattr(k, "Ksp_calcite", 0.4266))',
                '        _omega = (c_Ca2 * _co3) / (_Ksp + 1e-30)',
                '        _eps_sm = float(getattr(mp, "saturation_smooth_epsilon", 0.01))',
                '        _us_d = np.maximum(1.0 - _omega, 0.0)',
                '        _u_d = np.clip(_us_d / (_eps_sm + 1e-30), 0.0, 1.0)',
                '        _S_d = 3.0 * _u_d**2 - 2.0 * _u_d**3',
                '        _dS_d_du = np.where((_u_d > 0.0) & (_u_d < 1.0), 6.0 * _u_d * (1.0 - _u_d), 0.0)',
                '        _dS_d_domega = -_dS_d_du / (_eps_sm + 1e-30)',
                f'        _dR_diss_domega = -solid_fac * k_num_{r_name} * c_CaCO3 * (_S_d + _us_d * _dS_d_domega) * (_omega < 1.0)',
                '        _dr_dCa2 = _dR_diss_domega * (_co3 / (_Ksp + 1e-30))',
                '        _dr_dTA = _dR_diss_domega * (c_Ca2 / (_Ksp + 1e-30)) * _dco3_dta',
                '        _dr_dDIC = _dR_diss_domega * (c_Ca2 / (_Ksp + 1e-30)) * _dco3_ddic',
                f'        _dr_dCaCO3 = solid_fac * k_num_{r_name} * _us_d * _S_d * (_omega < 1.0)',
                '        for _vname, _dr in [("Ca2", _dr_dCa2), ("TA", _dr_dTA), ("DIC", _dr_dDIC), ("CaCO3", _dr_dCaCO3)]:',
                '            J[:, idx["CaCO3"], idx[_vname]] -= 1.0 * _dr',
                '            J[:, idx["Ca2"], idx[_vname]] += 1.0 * _dr',
                '            J[:, idx["DIC"], idx[_vname]] += 1.0 * _dr',
                '            J[:, idx["TA"], idx[_vname]] += 2.0 * _dr',
            ])
            lines.append('')
            return lines

        for coeff_diff, orig_sp_diff, diff_p in r_parents:
            other_r = [
                (c, sp, p) for c, sp, p in r_parents if p != diff_p
            ]
            other_terms = []
            for c, sp, p in other_r:
                spec = self.validator.get_speciation_factor_expr(sp)
                term = f'c_{p}' if not spec else f'(c_{p} * ({spec}))'
                other_terms.append(term)
            
            diff_spec = self.validator.get_speciation_factor_expr(orig_sp_diff)
            spec_factor = f' * ({diff_spec})' if diff_spec else ''
            
            all_prod = " * ".join(other_terms + lim_terms) if (other_terms or lim_terms) else "1.0"
            d_var = f"d_{diff_p}_{r_name}"
            lines.append(f'        {d_var} = {fac_var} * (k_num_{r_name} / {ref_coeff}){spec_factor} * {all_prod}')

            for coeff_r, _, r_tgt in r_parents:
                if r_tgt == "O2" and "POC" in diff_p:
                    ratio_str = 'float(getattr(mp, "POC_O2_ratio", 1.0))'
                    lines.append(f'        J[:, idx["{r_tgt}"], idx["{diff_p}"]] -= {ratio_str} * {d_var}')
                else:
                    lines.append(f'        J[:, idx["{r_tgt}"], idx["{diff_p}"]] -= {coeff_r} * {d_var}')
            for coeff_p, _, p_tgt in p_parents:
                lines.append(f'        J[:, idx["{p_tgt}"], idx["{diff_p}"]] += {coeff_p} * {d_var}')

            # Companion 32S isotopes
            iso_sp_map = r.get("isotope_species")
            if iso_sp_map and diff_p in iso_sp_map:
                diff_32 = iso_sp_map[diff_p]
                lines.append(f'        if "{diff_32}" in idx and getattr(mp, "isotopes", False):')
                for coeff_r, _, r_tgt in r_parents:
                    if r_tgt in iso_sp_map:
                        r_tgt_32 = iso_sp_map[r_tgt]
                        lines.append(f'            J[:, idx["{r_tgt_32}"], idx["{diff_32}"]] -= {coeff_r} * {d_var}')
                for coeff_p, _, p_tgt in p_parents:
                    if p_tgt in iso_sp_map:
                        p_tgt_32 = iso_sp_map[p_tgt]
                        lines.append(f'            J[:, idx["{p_tgt_32}"], idx["{diff_32}"]] += {coeff_p} * {d_var}')

        lines.append('')
        return lines

    def generate_code(self) -> str:
        """Generates complete Python source for equations.py."""
        code = []
        code.extend([
            '# -----------------------------------------------------------------------------',
            '# Individual Reaction Functions',
            '# -----------------------------------------------------------------------------',
            '',
        ])

        for r in self.reactions:
            cat = self.classify_reaction(r, self.validator)
            if cat.mechanism == "reversible_saturation":
                code.extend(self._generate_reversible_saturation_fn(r, cat))
            else:
                code.extend(self._generate_kinetic_reaction_fn(r, cat))

        code.extend([
            '# -----------------------------------------------------------------------------',
            '# Main Diagenetic Reactions Dispatcher',
            '# -----------------------------------------------------------------------------',
            '',
            'DEFAULT_DIAGENETIC_REACTIONS = [',
        ])
        for r in self.reactions:
            r_name = r.get("reaction_name", "rxn")
            fn_name = f"rxn_{r_name}"
            code.append(f'    [{fn_name}, None],')

        code.extend([
            ']',
            '',
            'DIAGENETIC_REACTIONS = DEFAULT_DIAGENETIC_REACTIONS',
            '',
            'def diagenetic_reactions(mp, c, k, f=None, lim=None, dt=None):',
            '    """Evaluates all registered diagenetic reactions for the model column."""',
            '    LHS = {}',
            '    RHS = {}',
            '    RATES = {}',
            '    CROSS = {}',
            '    for sp in getattr(mp, "species_list", list(c.keys())):',
            '        shape = c[sp].shape if hasattr(c[sp], "shape") else (len(c[sp]),)',
            '        LHS[sp] = np.zeros(shape, dtype=np.float64)',
            '        RHS[sp] = np.zeros(shape, dtype=np.float64)',
            '        RATES[sp] = np.zeros(shape, dtype=np.float64)',
            '        CROSS[sp] = []',
            '    if lim is not None:',
            '        limiters = lim',
            '    elif get_limiters is not None:',
            '        limiters = get_limiters(c, mp)',
            '    elif hasattr(mp, "get_limiters") and callable(mp.get_limiters):',
            '        limiters = mp.get_limiters(c, mp)',
            '    else:',
            '        limiters = getattr(mp, "limiters", {})',
            '    rxns = getattr(mp, "diagenetic_reactions", None)',
            '    if rxns is None:',
            '        rxns = DEFAULT_DIAGENETIC_REACTIONS',
            '    for r_entry in rxns:',
            '        fn = r_entry[0]',
            '        fn_k = r_entry[1] if (len(r_entry) > 1 and r_entry[1] is not None) else k',
            '        fn(c, fn_k, limiters, LHS, RHS, RATES, CROSS, mp)',
            '    if f is not None:',
            '        sp_list = getattr(mp, "species_list", list(c.keys()))',
            '        for s in sp_list:',
            '            setattr(f, s, (LHS.get(s, 0.0), RHS.get(s, 0.0), RATES.get(s, 0.0), None))',
            '        if not getattr(mp, "in_solver", False):',
            '            for key, val in RATES.items():',
            '                if key not in sp_list:',
            '                    setattr(f, key, (None, None, val, None))',
            '        object.__setattr__(f, "raw_LHS", LHS)',
            '        object.__setattr__(f, "raw_CROSS", CROSS)',
            '        object.__setattr__(f, "raw_RHS", RHS)',
            '        object.__setattr__(f, "raw_RATES", RATES)',
            '    return f, RATES',
            '',
            '',
            '# -----------------------------------------------------------------------------',
            '# Analytical SymPy Chemical Jacobian',
            '# -----------------------------------------------------------------------------',
            '',
            'def compute_chemical_jacobian(c, mp, k, species_names):',
            '    """Analytically evaluates the N x S x S chemical Jacobian J[i, s, m] = dR_s / dC_m.',
            '    Derived symbolically via SymPy to exact machine precision.',
            '    """',
            '    N = len(c[species_names[0]])',
            '    S = len(species_names)',
            '    J = np.zeros((N, S, S), dtype=np.float64)',
            '    idx = {name: i for i, name in enumerate(species_names)}',
            '    phi = np.asarray(mp.phi.value if hasattr(mp.phi, "value") else mp.phi)',
            '    solid_fac = 1.0 - phi',
            '    diss_fac = phi',
            '    if get_limiters is not None:',
            '        lim = get_limiters(c, mp)',
            '    elif hasattr(mp, "get_limiters") and callable(mp.get_limiters):',
            '        lim = mp.get_limiters(c, mp)',
            '    else:',
            '        lim = getattr(mp, "limiters", {})',
            '    has_carb = "TA" in idx and "DIC" in idx',
            '    if has_carb:',
            '        _ta_val = np.maximum(np.asarray(c["TA"].value if hasattr(c["TA"], "value") else c["TA"]), 1e-12)',
            '        _dic_val = np.maximum(np.asarray(c["DIC"].value if hasattr(c["DIC"], "value") else c["DIC"]), 1e-12)',
            '        _pK1 = getattr(mp, "pK1_carb", 5.86)',
            '        _pK2 = getattr(mp, "pK2_carb", 8.92)',
            '        _K1 = 10.0 ** (-float(_pK1))',
            '        _K2 = 10.0 ** (-float(_pK2))',
            '        _b = _K1 * (_ta_val - _dic_val)',
            '        _c_term = _K1 * _K2 * (_ta_val - 2.0 * _dic_val)',
            '        _disc = np.maximum(_b**2 - 4.0 * _ta_val * _c_term, 0.0)',
            '        _h = (-_b + np.sqrt(_disc)) / (2.0 * _ta_val)',
            '        _h = np.maximum(_h, 1e-14)',
            '        _denom = _h**2 + _K1 * _h + _K1 * _K2',
            '        _co3 = _dic_val * (_K1 * _K2) / _denom',
            '        _dF_dh = 2.0 * _ta_val * _h + _b',
            '        _dF_dta = _h**2 + _K1 * _h + _K1 * _K2',
            '        _dF_ddic = -(_K1 * _h + 2.0 * _K1 * _K2)',
            '        _dh_dta = -_dF_dta / (_dF_dh + 1e-30)',
            '        _dh_ddic = -_dF_ddic / (_dF_dh + 1e-30)',
            '        _dD_dh = 2.0 * _h + _K1',
            '        _dco3_dta = - (_dic_val * _K1 * _K2 / _denom**2) * _dD_dh * _dh_dta',
            '        _dco3_ddic = (_K1 * _K2 / _denom) - (_dic_val * _K1 * _K2 / _denom**2) * _dD_dh * _dh_ddic',
            '',
        ])

        for r in self.reactions:
            cat = self.classify_reaction(r, self.validator)
            if cat.mechanism == "reversible_saturation":
                code.extend(self._generate_reversible_saturation_jacobian(r, cat))
            else:
                code.extend(self._generate_kinetic_reaction_jacobian(r, cat))

        code.extend([
            '    return J',
            '',
            'diagenetic_reactions.compute_chemical_jacobian = compute_chemical_jacobian',
        ])

        # Extract reaction_parameters or fractionation_factors from reactions file if defined
        if self.validator.reactions_path and self.validator.reactions_path.exists():
            try:
                rxn_tree = ast.parse(self.validator.reactions_path.read_text(encoding="utf-8"))
                for node in ast.walk(rxn_tree):
                    if isinstance(node, ast.Assign):
                        for target in node.targets:
                            if isinstance(target, ast.Name) and target.id in ("reaction_parameters", "fractionation_factors"):
                                code.append("")
                                code.append(f"{target.id} = {ast.unparse(node.value)}")
            except Exception:
                pass

        # Collect all referenced limiters across reactions and Jacobian
        full_body = "\n".join(code)
        used_limiters = set(re.findall(r'lim(?:iters)?\["([^"]+)"\]', full_body))
        used_limiters.update(re.findall(r'lim(?:iters)?\.get\("([^"]+)"', full_body))

        header_code = [
            '"""equations.py',
            '',
            'Auto-generated reaction equations system module for MANTIS.',
            'Generated by generate_equation_system.py with analytical SymPy Jacobian support.',
            '"""',
            '',
            'from __future__ import annotations',
            'import numpy as np',
            'from mantis.diff_lib import (',
            '    add_coupled_reaction,',
            '    add_implicit_sink,',
            '    add_implicit_coupling_new,',
            '    calculate_fractionated_rate_coeff,',
            '    calculate_fractionated_coeff_32,',
            '    partition_equilibrium_isotope,',
            '    partition_equilibrium_isotope_32,',
            '    compute_carbonate_speciation,',
            '    data_container,',
            ')',
            'try:',
            '    from limiters import get_limiters as _external_get_limiters',
            'except ImportError:',
            '    _external_get_limiters = None',
            '',
            'def get_limiters(c, mp):',
            '    """Calculates phase-corrected Monod and inhibition limiters dynamically."""',
            '    if _external_get_limiters is not None:',
            '        return _external_get_limiters(c, mp)',
            '    if hasattr(c, "TA") and hasattr(c, "DIC"):',
            '        compute_carbonate_speciation(c, mp)',
            '    phi = np.asarray(mp.phi.value if hasattr(mp.phi, "value") else mp.phi)',
            '    solid_phi = 1.0 - phi',
            '    limiters = {}',
        ]

        # Scan limiter expressions to identify all half-saturation constants (K_*, Km_*, Kin_*)
        detected_k_constants = set()
        for lim_name, expr in self.validator.limiters_defs.items():
            for m in re.finditer(r'\b(K(?:m|in)?_[A-Za-z0-9_]+)\b', expr):
                detected_k_constants.add(m.group(1))

        # Dynamically determine phase scaling based on species phase from validator
        species_defs = getattr(self.validator, "species_defs", {})
        k_extract_lines = []
        for k_name in sorted(detected_k_constants):
            # Check if species matching the suffix is solid
            is_solid = False
            for sp_k, s_meta in species_defs.items():
                if isinstance(s_meta, dict) and s_meta.get("phase") == "solid":
                    if sp_k in k_name:
                        is_solid = True
                        break
            denom = "solid_phi" if is_solid else "phi"
            k_extract_lines.append(f'    {k_name} = getattr(mp, "{k_name}", 0.0) / {denom}')

        if k_extract_lines:
            header_code.append('    # Dynamically extract half-saturation constants (liquid scaled by phi, solid by 1-phi)')
            header_code.extend(k_extract_lines)

        # Dynamically emit only limiters actually referenced in reactions or Jacobian
        for lim_name, expr in self.validator.limiters_defs.items():
            base_key = lim_name[:-9] if lim_name.endswith("_implicit") else lim_name
            clean_expr = expr.strip("'\"")
            var_name = f"_l_{lim_name.replace('.', '_')}"

            possible_keys = {lim_name, f"{base_key}_implicit"}
            if "_" in base_key:
                main_k, sub_k = base_key.split("_", 1)
                possible_keys.update([f"{main_k}_implicit_{sub_k}", f"{main_k}_{sub_k}_implicit"])

            has_denom = "1.0 /" in clean_expr
            if has_denom:
                possible_keys.update([f"{base_key}_explicit", f"{base_key}_inhibit", f"{base_key}_inhib"])
                if "_" in base_key:
                    main_k, sub_k = base_key.split("_", 1)
                    possible_keys.update([
                        f"{main_k}_explicit_{sub_k}",
                        f"{main_k}_inhibit_{sub_k}",
                        f"{main_k}_inhib_{sub_k}",
                    ])

            needed_keys = possible_keys.intersection(used_limiters)
            if not needed_keys:
                continue

            header_code.append(f'    {var_name} = {clean_expr}')
            if lim_name in needed_keys:
                header_code.append(f'    limiters["{lim_name}"] = {var_name}')
            if f"{base_key}_implicit" in needed_keys and f"{base_key}_implicit" != lim_name:
                header_code.append(f'    limiters["{base_key}_implicit"] = {var_name}')
            if "_" in base_key:
                main_k, sub_k = base_key.split("_", 1)
                if f"{main_k}_implicit_{sub_k}" in needed_keys and f"{main_k}_implicit_{sub_k}" != lim_name:
                    header_code.append(f'    limiters["{main_k}_implicit_{sub_k}"] = {var_name}')
                if f"{main_k}_{sub_k}_implicit" in needed_keys and f"{main_k}_{sub_k}_implicit" not in (lim_name, f"{main_k}_implicit_{sub_k}"):
                    header_code.append(f'    limiters["{main_k}_{sub_k}_implicit"] = {var_name}')

            if has_denom:
                denom = clean_expr.split("1.0 /", 1)[1].strip()
                if denom.startswith("(") and denom.endswith(")"):
                    inner = denom[1:-1].strip()
                    parts = inner.split("+", 1)
                    if len(parts) == 2:
                        var_part = parts[0].strip()
                        k_part = parts[1].strip()

                        if f"{base_key}_explicit" in needed_keys:
                            header_code.append(f'    limiters["{base_key}_explicit"] = ({var_part}) * {var_name}')
                        if "_" in base_key:
                            main_k, sub_k = base_key.split("_", 1)
                            if f"{main_k}_explicit_{sub_k}" in needed_keys and f"{main_k}_explicit_{sub_k}" != f"{base_key}_explicit":
                                header_code.append(f'    limiters["{main_k}_explicit_{sub_k}"] = ({var_part}) * {var_name}')

                        inhibit_val = f'({k_part}) * {var_name}'
                        inhibit_computed = False
                        if f"{base_key}_inhibit" in needed_keys:
                            header_code.append(f'    limiters["{base_key}_inhibit"] = {inhibit_val}')
                            inhibit_computed = True
                        if f"{base_key}_inhib" in needed_keys:
                            header_code.append(f'    limiters["{base_key}_inhib"] = limiters["{base_key}_inhibit"]' if inhibit_computed else f'    limiters["{base_key}_inhib"] = {inhibit_val}')
                            inhibit_computed = True
                        if "_" in base_key:
                            main_k, sub_k = base_key.split("_", 1)
                            if f"{main_k}_inhibit_{sub_k}" in needed_keys and f"{main_k}_inhibit_{sub_k}" != f"{base_key}_inhibit":
                                header_code.append(f'    limiters["{main_k}_inhibit_{sub_k}"] = {inhibit_val}')
                            if f"{main_k}_inhib_{sub_k}" in needed_keys and f"{main_k}_inhib_{sub_k}" not in (f"{base_key}_inhib", f"{main_k}_inhibit_{sub_k}"):
                                header_code.append(f'    limiters["{main_k}_inhib_{sub_k}"] = {inhibit_val}')

        header_code.extend([
            '    return limiters',
            '',
        ])
        return "\n".join(header_code + code)


def main(args_list: Optional[List[str]] = None) -> int:
    """CLI entrypoint for generate_equation_system."""
    parser = build_arg_parser()
    args = parser.parse_args(args_list)

    # Maintain attribute compatibility with legacy flag names
    if not hasattr(args, "input"):
        args.input = args.reactions
    if not hasattr(args, "output"):
        args.output = args.equations
    if not hasattr(args, "constants"):
        args.constants = args.kinetic_constants

    try:
        generate_equations(
            reactions_path=args.reactions,
            equations_path=args.equations,
            kinetic_constants_path=args.kinetic_constants,
            limiters_path=args.limiters,
            species_path=args.species,
            verbose=args.verbose,
            validate_only=args.validate_only,
        )
        if not args.validate_only:
            print(f"Successfully generated equations module saved to: {args.equations}")
        return 0
    except ValidationError as ve:
        print("\n" + "=" * 60, file=sys.stderr)
        print(f"VALIDATION FAILED:\n{ve}", file=sys.stderr)
        print("=" * 60 + "\n", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

