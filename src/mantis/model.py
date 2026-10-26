"""model.py

Unified, declarative DiageneticModel interface for MANTIS reaction-transport simulations.

Provides a clean object-oriented entry point that:
- Automatically resolves paths relative to the calling script.
- Centralizes species, phases, dynamic speciation, and default boundary conditions.
- Extracts pH from species.py (Hplus['bc']['pH']) as a single source of truth.
- Sets bioirrigation explicitly to zero by default.
- Provides built-in CLI argument parsing (--help, --solver, --video).
"""

from __future__ import annotations

import argparse
import ast
import copy
import importlib.util
import inspect
import os
from pathlib import Path
import sys
import re
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import pint

from mantis.diff_lib import (
    VariableArray as CellVariable,
    check_peclet_numbers,
    compute_bio_irrigation_alpha,
    compute_exponential_db,
    compute_sigmoidal_db,
    data_container,
    diff_coeff,
    get_l_mass,
    get_total_delta,
    make_grid_from_sections,
    make_grid2,
    read_state,
    save_data,
    save_state,
)
from mantis.generate_equation_system import resolve_and_prepare_equations
from mantis.live_plot_lib import LivePlotter
from mantis.monolithic_solver import run_non_steady_state_solver_monolithic as run_solver


class DiageneticModelResult:
    """Encapsulates the output of a completed DiageneticModel simulation."""

    def __init__(
        self,
        mp: data_container,
        c: data_container,
        k: data_container,
        species_list: List[str],
        z: np.ndarray,
        D_mol: data_container,
        diagenetic_reactions: Any,
        converged: str,
        step: int,
        total_time: float,
        df: Any = None,
        fqfn: Optional[str] = None,
    ):
        self.mp = mp
        self.c = c
        self.k = k
        self.species_list = species_list
        self.z = z
        self.D_mol = D_mol
        self.diagenetic_reactions = diagenetic_reactions
        self.converged = converged
        self.step = step
        self.total_time = total_time
        self.df = df
        self.fqfn = fqfn

    def save_csv(self) -> Tuple[Any, str]:
        """Saves tabular concentration and rate data to CSV."""
        df, fqfn = save_data(
            self.mp, self.c, self.k, self.species_list, self.z, self.D_mol, self.diagenetic_reactions
        )
        self.df = df
        self.fqfn = fqfn
        return df, fqfn

    def save_state(self, filepath: Optional[Union[str, Path]] = None) -> str:
        """Saves complete CellVariable array state to a compressed .npz archive."""
        if filepath is None:
            caller_dir = getattr(self.mp, "base_dir", None) or getattr(self.mp, "_caller_dir", None) or Path.cwd()
            filepath = Path(caller_dir) / f"{self.mp.plot_name}_state.npz"
        filepath_str = str(filepath)
        save_state(self.c, filepath_str, mp=self.mp)
        return filepath_str

    def plot_summary(self, layout_file: Optional[Union[str, Path]] = None, outfile: Optional[str] = None) -> Optional[str]:
        """Renders the final summary plot PDF."""
        import mantis.plot_data_new as pdn
        layout = layout_file or getattr(self.mp, "layout_file", "plot_layout.py")
        p_layout = Path(layout)
        if not p_layout.is_absolute():
            if (Path.cwd() / p_layout).exists():
                p_layout = Path.cwd() / p_layout
            elif hasattr(self.mp, "_caller_dir") and (Path(self.mp._caller_dir) / p_layout.name).exists():
                p_layout = Path(self.mp._caller_dir) / p_layout.name
            elif (Path.cwd() / p_layout.name).exists():
                p_layout = Path.cwd() / p_layout.name
            elif hasattr(self.mp, "_caller_dir") and (Path(self.mp._caller_dir) / p_layout).exists():
                p_layout = Path(self.mp._caller_dir) / p_layout

        if not p_layout.exists():
            caller_dir = getattr(self.mp, "_caller_dir", Path.cwd())
            print(
                f"\n[MANTIS ERROR] Plot layout file '{layout}' not found.\n"
                f"Searched in:\n"
                f"  • Current working directory: {Path.cwd()}\n"
                f"  • Model directory:           {caller_dir}\n"
                f"Plot summary rendering skipped.\n",
                file=sys.stderr,
            )
            return None

        if self.df is not None:
            caller_dir = getattr(self.mp, "base_dir", None) or getattr(self.mp, "_caller_dir", None) or Path.cwd()
            pdf_out = outfile or str(Path(caller_dir) / f"{self.mp.plot_name}.pdf")
            try:
                plt_desc = pdn.load_layout_from_file(self.df, str(p_layout))
                display_len = getattr(self.mp, "display_length", None) or getattr(self.mp, "max_depth", 1.0)
                pdn.plot(self.df, display_len, outfile=pdf_out, show=False, plot_description=plt_desc)
                return pdf_out
            except Exception as e:
                print(f"[DiageneticModel] Plot rendering skipped or failed: {e}")
                return None
        return None


class DiageneticModel:
    """Declarative Diagenetic Model Runner."""

    def __init__(
        self,
        species: Union[str, Path] = "species.py",
        reactions: Optional[Union[str, Path]] = "reactions.py",
        constants: Optional[Union[str, Path]] = "reaction_constants.py",
        equations: Optional[Union[str, Path]] = None,
        limiters: Optional[Union[str, Path]] = None,
        layout: Union[str, Path] = "plot_layout.py",
        experiment: Optional[str] = None,
        base_dir: Optional[Union[str, Path]] = None,
    ):
        # 1. Determine caller directory and caller frame automatically
        caller_frame = None
        for frame_info in inspect.stack()[1:]:
            frame_file = frame_info.filename
            frame_name = Path(frame_file).name
            if (
                frame_name not in ("model.py", "__init__.py")
                and "_pytest" not in frame_file
                and "pluggy" not in frame_file
                and frame_name not in ("pytest", "pytest.py")
            ):
                caller_frame = frame_info
                break

        if base_dir is not None:
            self.base_dir = Path(base_dir).resolve()
        elif caller_frame is not None and caller_frame.filename != "<stdin>":
            self.base_dir = Path(caller_frame.filename).resolve().parent
        else:
            self.base_dir = Path.cwd()

        # Derive experiment name automatically from caller filename stem
        if experiment is not None:
            self.experiment = experiment
        elif caller_frame is not None and caller_frame.filename not in ("<stdin>", "") and not Path(caller_frame.filename).name.startswith("pytest"):
            self.experiment = Path(caller_frame.filename).stem
        elif sys.argv and sys.argv[0] and not Path(sys.argv[0]).name.startswith("pytest"):
            self.experiment = Path(sys.argv[0]).stem
        else:
            self.experiment = "diagenetic_run"

        # 2. Store original user inputs for clear diagnostics
        self._orig_species = species
        self._orig_reactions = reactions
        self._orig_constants = constants
        self._orig_equations = equations
        self._orig_limiters = limiters
        self._orig_layout = layout

        # 3. Resolve filenames relative to cwd or base_dir (wherever the model resides)
        self.species_file = self._resolve(species)
        rxn_fallbacks = [f"{self.experiment}_reactions.py"] if str(reactions) == "reactions.py" else None
        self.reactions_file = self._resolve(reactions, fallbacks=rxn_fallbacks) if reactions else None
        self.constants_file = self._resolve(constants) if constants else None
        self.equations_file = self._resolve(equations) if equations else None
        self.limiters_file = self._resolve(limiters) if limiters else None
        self.layout_file = self._resolve(layout)

        # 4. If running as a script (not querying CLI help/options or in pytest), validate files immediately
        is_cli_query = any(arg in sys.argv for arg in ("-h", "--help", "--solver", "--video", "--run"))
        in_pytest = "pytest" in sys.modules or (sys.argv and Path(sys.argv[0]).name.startswith("pytest"))
        if not is_cli_query and not in_pytest:
            self.validate_simulation_files()

        self.ureg = pint.UnitRegistry()
        self.Q_ = self.ureg.Quantity

        # Master parameter container
        self.mp = data_container()
        self.p_dict: Dict[str, Any] = {}
        self._user_configured_keys: set = set()

        # Default parameters
        self._init_defaults()

    def _exit_with_missing_file_error(self, filename: Union[str, Path], file_description: str) -> None:
        """Prints a user-friendly error message and exits when a simulation file cannot be found."""
        print(f"\n[MANTIS ERROR] Could not find {file_description}: '{filename}'", file=sys.stderr)
        print("Searched locations:", file=sys.stderr)
        print(f"  • Current working directory: {Path.cwd()}", file=sys.stderr)
        if self.base_dir.resolve() != Path.cwd().resolve():
            print(f"  • Model directory:           {self.base_dir}", file=sys.stderr)
        print(
            f"\nPlease ensure that '{filename}' is located in the current working directory "
            f"or the directory where the model script resides ({self.base_dir}).\n",
            file=sys.stderr,
        )
        sys.exit(1)

    def validate_simulation_files(self) -> None:
        """Validates that all required files describing the simulation exist in cwd or model directory."""
        if not self.species_file.exists():
            self._exit_with_missing_file_error(self._orig_species, "species definitions file")
        if self._orig_reactions is not None and self.reactions_file is not None and not self.reactions_file.exists():
            self._exit_with_missing_file_error(self._orig_reactions, "declarative reactions file")
        if self._orig_constants is not None and self.constants_file is not None and not self.constants_file.exists():
            if not self.is_diffusion_only:
                self._exit_with_missing_file_error(self._orig_constants, "reaction constants file")
        if self._orig_limiters is not None and self.limiters_file is not None and not self.limiters_file.exists():
            self._exit_with_missing_file_error(self._orig_limiters, "limiters file")

    def _resolve(self, path_like: Union[str, Path], fallbacks: Optional[List[str]] = None) -> Path:
        p = Path(path_like)
        # Check absolute path
        if p.is_absolute() and p.exists():
            return p.resolve()

        # Check in current working directory
        if (Path.cwd() / p).exists():
            return (Path.cwd() / p).resolve()

        # Check in model directory (base_dir)
        if (self.base_dir / p).exists():
            return (self.base_dir / p).resolve()

        # Check fallbacks if provided
        if fallbacks:
            for fb in fallbacks:
                fb_p = Path(fb)
                if (Path.cwd() / fb_p).exists():
                    return (Path.cwd() / fb_p).resolve()
                if (self.base_dir / fb_p).exists():
                    return (self.base_dir / fb_p).resolve()

        # Fallback to model directory or cwd
        if self.base_dir != Path.cwd():
            return (self.base_dir / p).resolve()
        return (Path.cwd() / p).resolve()

    def _init_defaults(self) -> None:
        """Initializes default numerical solver, video output, and experiment metadata parameters."""
        Q_ = self.Q_
        self.p_dict.update({
            # Experiment metadata
            "experiment": self.experiment,
            "plot_name": self.experiment,
            "state_data": None,
            "layout_file": str(self.layout_file),
            "start_time": 0,

            # Video recording options
            "video_dt": None,
            "video_dpi": 120,
            "video_queue_size": 100,
            "video_resume_threshold": 50,
            "video_workers": None,

            # Solver parameters
            "solver": "monolithic",
            "t_end": Q_("1 kyr").to("seconds").magnitude,
            "dt_min": Q_("1 minute").to("seconds").magnitude,
            "dt_init": Q_("1 day").to("seconds").magnitude,
            "dt_max": Q_("1 year").to("seconds").magnitude,
            "tolerance": 1e-12,
            "dt_tolerance": 1e-5,

            # Monolithic Newton Settings
            "max_newton_iters": 12,
            "newton_tol": 1e-3,
            "newton_norm": "linf",
            "newton_target_optimal": 4,
            "newton_max_acceptable": 8,
            "use_numerical_jacobian": False,
        })

    def _convert_value(self, val: Any) -> Any:
        Q_ = self.Q_
        if isinstance(val, (list, tuple)):
            converted = [self._convert_value(item) for item in val]
            return tuple(converted) if isinstance(val, tuple) else converted
        if isinstance(val, dict):
            return {k: self._convert_value(v) for k, v in val.items()}
        if isinstance(val, str):
            try:
                q = Q_(val)
                if hasattr(q, "dimensionality"):
                    dim = q.dimensionality
                    if dim == Q_("1 second").dimensionality:
                        return q.to("seconds").magnitude
                    elif dim == Q_("1 meter").dimensionality:
                        return q.to("meter").magnitude
                    elif dim == Q_("1 m/s").dimensionality:
                        return q.to("m/s").magnitude
                    elif dim == Q_("1 m^2/s").dimensionality:
                        return q.to("m^2/second").magnitude
                    elif dim == Q_("1 / second").dimensionality:
                        return q.to("1/second").magnitude
                    elif dim == Q_("1 mol/m^3").dimensionality:
                        return q.to("mol/m^3").magnitude
                    elif dim == Q_("1 mol/(m^2 * second)").dimensionality:
                        return q.to("mol/(m^2 * second)").magnitude
                    elif q.dimensionless and not isinstance(q.magnitude, str):
                        return q.magnitude
            except Exception:
                pass
        elif hasattr(val, "dimensionality"):
            dim = val.dimensionality
            if dim == Q_("1 second").dimensionality:
                return val.to("seconds").magnitude
            elif dim == Q_("1 meter").dimensionality:
                return val.to("meter").magnitude
            elif dim == Q_("1 m/s").dimensionality:
                return val.to("m/s").magnitude
            elif dim == Q_("1 m^2/s").dimensionality:
                return val.to("m^2/second").magnitude
            elif dim == Q_("1 / second").dimensionality:
                return val.to("1/second").magnitude
            elif dim == Q_("1 mol/m^3").dimensionality:
                return val.to("mol/m^3").magnitude
            elif dim == Q_("1 mol/(m^2 * second)").dimensionality:
                return val.to("mol/(m^2 * second)").magnitude
            elif val.dimensionless:
                return val.to_base_units().magnitude
        return val

    def configure(self, **kwargs) -> DiageneticModel:
        """Sets or overrides model configuration parameters with automatic unit handling.

        Any parameter passed with units (as a string or Pint Quantity) is automatically
        detected and converted to standard SI simulation units based on its physical dimensionality:
        - [time] -> seconds (e.g. t_end='2 kyr', dt_max='1 year', dt_min='1 second')
        - [length] -> meters (e.g. max_depth='0.5 m', initial_spacing='0.01 mm')
        - [length] / [time] -> m/s (e.g. w='0.2 cm/yr')
        - [length]^2 / [time] -> m^2/s (e.g. BT0='4 cm^2/year')
        - 1 / [time] -> 1/s (e.g. BI0='0.1 /year')
        - [substance] / [length]^3 -> mol/m^3 (e.g. concentration='0.001 umol/cm^3')
        - [substance] / ([length]^2 * [time]) -> mol/(m^2*s) (e.g. flux='1.2 mol/(m^2*yr)')

        Nested structures (dicts, lists, tuples) are handled recursively.
        """
        self._user_configured_keys.update(kwargs.keys())
        for key, val in kwargs.items():
            if key in ("video", "video_dt"):
                self.p_dict["video_dt"] = None if str(val).lower() in ("none", "null", "false") else val
            elif key == "pH":
                self.p_dict["pH"] = float(val)
            elif key in ("steady_state_tolerance", "dt_tolerance"):
                f_val = float(val)
                self.p_dict["steady_state_tolerance"] = f_val
                self.p_dict["dt_tolerance"] = f_val
            elif key == "check_reactant_coexistence":
                self.p_dict["check_reactant_coexistence"] = bool(val)
            elif key in ("diffusive_boundary_layer", "dbl"):
                val_m = self._convert_value(val)
                self.p_dict["diffusive_boundary_layer"] = float(val_m)
                self.p_dict["dbl"] = float(val_m)
            elif key == "bioturbation" and isinstance(val, dict):
                c_val = self._convert_value(val)
                if "BT0" in c_val:
                    self.p_dict["BT0"] = float(c_val["BT0"])
                if "depth" in c_val:
                    self.p_dict["BT_depth"] = float(c_val["depth"])
                if "attenuation" in c_val:
                    self.p_dict["BT_attenuation"] = float(c_val["attenuation"])
                if "decay" in c_val:
                    self.p_dict["BT_attenuation"] = float(c_val["decay"])
                if "profile" in c_val:
                    self.p_dict["BT_profile"] = str(c_val["profile"]).lower()
                self.p_dict["bioturbation"] = c_val
            elif key == "bioirrigation" and isinstance(val, dict):
                c_val = self._convert_value(val)
                if "BI0" in c_val:
                    self.p_dict["BI0"] = float(c_val["BI0"])
                if "depth" in c_val:
                    self.p_dict["BI_depth"] = float(c_val["depth"])
                if "attenuation" in c_val:
                    self.p_dict["BI_attenuation"] = float(c_val["attenuation"])
                if "decay" in c_val:
                    self.p_dict["BI_attenuation"] = float(c_val["decay"])
                if "profile" in c_val:
                    self.p_dict["BI_profile"] = str(c_val["profile"]).lower()
                self.p_dict["bioirrigation"] = c_val
            else:
                self.p_dict[key] = self._convert_value(val)

        self.mp.update(self.p_dict)
        return self

    def parse_cli_args(self, args: Optional[List[str]] = None) -> None:
        """Parses CLI arguments if requested, supporting category inspection and option overrides."""
        cur_solver = self.p_dict.get("solver", "monolithic")
        cur_tol = self.p_dict.get("newton_tol", 1e-3)
        cur_jac = self.p_dict.get("use_numerical_jacobian", False)
        cur_max_iters = self.p_dict.get("max_newton_iters", 8)
        cur_norm = self.p_dict.get("newton_norm", "linf")
        cur_opt_iters = self.p_dict.get("newton_target_optimal", 4)
        cur_max_acc_iters = self.p_dict.get("newton_max_acceptable", 6)

        cur_video = self.p_dict.get("video_dt", "None")
        cur_dpi = self.p_dict.get("video_dpi", 120)

        cur_tend = self.p_dict.get("t_end", "None")
        cur_dtmax = self.p_dict.get("dt_max", "None")
        cur_dtinit = self.p_dict.get("dt_init", "None")
        cur_dtmin = self.p_dict.get("dt_min", "None")
        cur_max_depth = self.p_dict.get("max_depth", "None")
        cur_disp_len = self.p_dict.get("display_length", "None")

        cur_dt_tol = self.p_dict.get("dt_tolerance", 1e-5)
        cur_vqueue = self.p_dict.get("video_queue_size", 100)
        cur_vresume = self.p_dict.get("video_resume_threshold", 50)
        cur_vworkers = self.p_dict.get("video_workers", None)

        solver_section = [
            ("--solver <backend>", cur_solver, "Backend solver choice ('monolithic' or 'coupled')"),
            ("--solver-tol <float>", cur_tol, "Newton convergence tolerance (e.g. 1e-3)"),
            ("--numerical-jacobian", cur_jac, "Use finite-difference Jacobian instead of SymPy analytical"),
            ("--max-newton-iters <int>", cur_max_iters, "Maximum Newton iterations allowed per timestep"),
            ("--newton-norm <choice>", cur_norm, "Residual norm formulation ('linf' or 'wrms')"),
            ("--newton-target-optimal <int>", cur_opt_iters, "Target Newton iterations for timestep expansion"),
            ("--newton-max-acceptable <int>", cur_max_acc_iters, "Max acceptable iterations before timestep reduction"),
            ("--dt-tolerance <float>", cur_dt_tol, "Adaptive timestep convergence tolerance (e.g. 1e-5)"),
            ("--force-regen", False, "Force regeneration of equation system and analytical Jacobian"),
        ]

        video_section = [
            ("--video <interval>", cur_video, "Live video recording frame interval (e.g. '10 year', 'none')"),
            ("--video-dpi <int>", cur_dpi, "Resolution (DPI) of recorded video frames"),
            ("--video-queue-size <int>", cur_vqueue, "Maximum video frame render queue size"),
            ("--video-resume-threshold <int>", cur_vresume, "Video frame resume queue threshold"),
            ("--video-workers <int>", cur_vworkers, "Dedicated background worker count for video rendering"),
        ]

        run_section = [
            ("--t-end <duration>", f"{cur_tend} s" if isinstance(cur_tend, (int, float)) else str(cur_tend), "Total duration of simulation (e.g. '10 kyr')"),
            ("--dt-max <timestep>", f"{cur_dtmax} s" if isinstance(cur_dtmax, (int, float)) else str(cur_dtmax), "Maximum adaptive timestep cap (e.g. '1 year')"),
            ("--dt-init <timestep>", f"{cur_dtinit} s" if isinstance(cur_dtinit, (int, float)) else str(cur_dtinit), "Initial timestep size (e.g. '1 day')"),
            ("--dt-min <timestep>", f"{cur_dtmin} s" if isinstance(cur_dtmin, (int, float)) else str(cur_dtmin), "Minimum allowable timestep before aborting"),
            ("--max-depth <length>", f"{cur_max_depth} m" if isinstance(cur_max_depth, (int, float)) else str(cur_max_depth), "Total depth of sediment model column (e.g. '10 m')"),
            ("--display-length <length>", f"{cur_disp_len} m" if isinstance(cur_disp_len, (int, float)) else str(cur_disp_len), "Plotted column depth for video and PDF summaries"),
        ]

        parser = argparse.ArgumentParser(
            description=(
                f"MANTIS Diagenetic Model: {self.experiment}\n\n"
                "Interactive Category Query:\n"
                "  Passing category flags without arguments displays all configurable options and current defaults:\n"
                "    %(prog)s --solver        Show all solver options and current settings\n"
                "    %(prog)s --video         Show all video output options\n"
                "    %(prog)s --run           Show all time, grid, and domain options\n"
            ),
            formatter_class=argparse.RawDescriptionHelpFormatter,
            add_help=True,
        )

        # Solver options
        parser.add_argument(
            "--solver",
            nargs="?",
            const="__SHOW__",
            metavar="BACKEND",
            help=f"Solver backend: 'monolithic' or 'coupled' (or query section: pass flag with no value)",
        )
        parser.add_argument(
            "--solver-tol",
            "--newton-tol",
            nargs="?",
            const="__SHOW__",
            metavar="TOL",
            help=f"Newton convergence tolerance (default: {cur_tol})",
        )
        parser.add_argument(
            "--numerical-jacobian",
            action="store_true",
            help=f"Use finite-difference Jacobian instead of SymPy analytical (default: {cur_jac})",
        )
        parser.add_argument(
            "--max-newton-iters",
            type=int,
            metavar="ITERS",
            help=f"Max Newton iterations per step (default: {cur_max_iters})",
        )
        parser.add_argument(
            "--newton-norm",
            choices=["linf", "wrms"],
            metavar="NORM",
            help=f"Residual norm metric: 'linf' or 'wrms' (default: {cur_norm})",
        )
        parser.add_argument(
            "--newton-target-optimal",
            type=int,
            metavar="ITERS",
            help=f"Target iterations for timestep expansion (default: {cur_opt_iters})",
        )
        parser.add_argument(
            "--newton-max-acceptable",
            type=int,
            metavar="ITERS",
            help=f"Max iterations before timestep reduction (default: {cur_max_acc_iters})",
        )
        parser.add_argument(
            "--dt-tolerance",
            type=float,
            metavar="TOL",
            help=f"Adaptive timestep tolerance (default: {cur_dt_tol})",
        )
        parser.add_argument(
            "--force-regen",
            "--rebuild-equations",
            action="store_true",
            help="Force regeneration of equations and analytical Jacobian",
        )

        # Video options
        parser.add_argument(
            "--video",
            nargs="?",
            const="__SHOW__",
            metavar="INTERVAL",
            help=f"Live video recording interval, e.g. '10 year' or 'none' (or query section: pass flag with no value)",
        )
        parser.add_argument(
            "--video-dpi",
            nargs="?",
            const="__SHOW__",
            metavar="DPI",
            help=f"Video DPI resolution (default: {cur_dpi})",
        )
        parser.add_argument(
            "--video-queue-size",
            type=int,
            metavar="SIZE",
            help=f"Max video queue size (default: {cur_vqueue})",
        )
        parser.add_argument(
            "--video-resume-threshold",
            type=int,
            metavar="SIZE",
            help=f"Video queue resume threshold (default: {cur_vresume})",
        )
        parser.add_argument(
            "--video-workers",
            type=int,
            metavar="COUNT",
            help="Worker count for video encoding",
        )

        # Run / Domain options
        parser.add_argument(
            "--run",
            nargs="?",
            const="__SHOW__",
            metavar="SETTINGS",
            help="Query time and domain run settings (pass flag with no value)",
        )
        parser.add_argument(
            "--t-end",
            metavar="DURATION",
            help=f"Simulation duration, e.g. '10 kyr' (default: {cur_tend} s)",
        )
        parser.add_argument(
            "--dt-max",
            metavar="TIMESTEP",
            help=f"Maximum timestep, e.g. '1 year' (default: {cur_dtmax} s)",
        )
        parser.add_argument(
            "--dt-init",
            metavar="TIMESTEP",
            help=f"Initial timestep, e.g. '1 day' (default: {cur_dtinit} s)",
        )
        parser.add_argument(
            "--dt-min",
            metavar="TIMESTEP",
            help=f"Minimum timestep, e.g. '1 second' (default: {cur_dtmin} s)",
        )
        parser.add_argument(
            "--max-depth",
            metavar="LENGTH",
            help=f"Sediment column depth, e.g. '10 m' (default: {cur_max_depth} m)",
        )
        parser.add_argument(
            "--display-length",
            metavar="LENGTH",
            help=f"Plotted display depth, e.g. '0.5 m' (default: {cur_disp_len} m)",
        )

        parsed, _ = parser.parse_known_args(args)

        # Check section queries first
        def _print_section(title: str, rows: List[Tuple[str, Any, str]]) -> None:
            print(f"\n[{self.experiment}] {title}:")
            for flag, cur_val, desc in rows:
                print(f"  {flag:<30} : {str(cur_val):<12} ({desc})")
            print()

        show_solver = parsed.solver in ("__SHOW__", "help", "?", "show") or parsed.solver_tol == "__SHOW__"
        show_video = parsed.video in ("__SHOW__", "help", "?", "show") or parsed.video_dpi == "__SHOW__"
        show_run = parsed.run in ("__SHOW__", "help", "?", "show")

        if show_solver or show_video or show_run:
            if show_solver:
                _print_section("Solver Configuration & Configurable Options", solver_section)
            if show_video:
                _print_section("Video Output Configuration & Configurable Options", video_section)
            if show_run:
                _print_section("Simulation Time & Domain Configuration & Configurable Options", run_section)
            sys.exit(0)

        # Apply solver overrides
        if parsed.solver is not None:
            if parsed.solver in ("monolithic", "coupled"):
                self.p_dict["solver"] = parsed.solver
            else:
                parser.error(f"argument --solver: invalid choice: '{parsed.solver}' (choose from 'monolithic', 'coupled')")

        if parsed.solver_tol is not None:
            self.p_dict["newton_tol"] = float(parsed.solver_tol)
        if parsed.numerical_jacobian:
            self.p_dict["use_numerical_jacobian"] = True
            self._user_configured_keys.add("use_numerical_jacobian")
        if parsed.max_newton_iters is not None:
            self.p_dict["max_newton_iters"] = parsed.max_newton_iters
        if parsed.newton_norm is not None:
            self.p_dict["newton_norm"] = parsed.newton_norm
        if parsed.newton_target_optimal is not None:
            self.p_dict["newton_target_optimal"] = parsed.newton_target_optimal
        if parsed.newton_max_acceptable is not None:
            self.p_dict["newton_max_acceptable"] = parsed.newton_max_acceptable
        if parsed.dt_tolerance is not None:
            self.configure(dt_tolerance=parsed.dt_tolerance)
        if parsed.force_regen:
            self.p_dict["force_regen"] = True

        # Apply video overrides
        if parsed.video is not None:
            self.p_dict["video_dt"] = None if parsed.video.lower() in ("none", "null", "false") else parsed.video
        if parsed.video_dpi is not None:
            self.p_dict["video_dpi"] = int(parsed.video_dpi)
        if parsed.video_queue_size is not None:
            self.p_dict["video_queue_size"] = int(parsed.video_queue_size)
        if parsed.video_resume_threshold is not None:
            self.p_dict["video_resume_threshold"] = int(parsed.video_resume_threshold)
        if parsed.video_workers is not None:
            self.p_dict["video_workers"] = int(parsed.video_workers)

        # Apply run / domain overrides
        if parsed.t_end is not None:
            self.configure(t_end=parsed.t_end)
        if parsed.dt_max is not None:
            self.configure(dt_max=parsed.dt_max)
        if parsed.dt_init is not None:
            self.configure(dt_init=parsed.dt_init)
        if parsed.dt_min is not None:
            self.configure(dt_min=parsed.dt_min)
        if parsed.max_depth is not None:
            self.configure(max_depth=parsed.max_depth)
        if parsed.display_length is not None:
            self.configure(display_length=parsed.display_length)

    def _load_species_definitions(self) -> Dict[str, Any]:
        """Loads and parses the species definitions dictionary from species.py."""
        if not self.species_file.exists():
            self._exit_with_missing_file_error(self._orig_species, "species definitions file")
        spec = importlib.util.spec_from_file_location("species_mod", self.species_file)
        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load species file: {self.species_file}")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        if not hasattr(mod, "species") or not isinstance(mod.species, dict):
            raise ValueError(f"Species file {self.species_file} must define a 'species' dictionary.")
        self._species_mod = mod
        return mod.species

    @property
    def is_diffusion_only(self) -> bool:
        """Determines if the simulation is purely diffusive without diagenetic reactions."""
        if self._orig_reactions is None or self.reactions_file is None:
            return True
        if hasattr(self, "_raw_reactions") and len(self._raw_reactions) == 0:
            return True
        if hasattr(self, "reactions_file") and self.reactions_file and self.reactions_file.exists():
            try:
                with open(self.reactions_file, "r") as f:
                    rxn_tree = ast.parse(f.read())
                for node in rxn_tree.body:
                    if isinstance(node, ast.Assign):
                        for target in node.targets:
                            if isinstance(target, ast.Name) and target.id == "reactions":
                                if isinstance(node.value, (ast.List, ast.Tuple)) and len(node.value.elts) == 0:
                                    return True
            except Exception:
                pass
        return False

    def check_requires_ph(self, species_defs: Dict[str, Any], speciation_defs: Dict[str, Any]) -> Tuple[bool, str]:
        """Detects whether the chemical model (species, speciation, reactions, or kinetic constants) requires pH."""
        ph_re = re.compile(r"\b(pH|Hplus|H\+)\b")

        # 1. Species definitions
        for sp_name, meta in species_defs.items():
            if sp_name in ("Hplus", "H+", "pH"):
                return True, f"species '{sp_name}' is defined"
            if isinstance(meta, dict):
                formula = meta.get("formula", "")
                if formula in ("H+", "Hplus"):
                    return True, f"species '{sp_name}' has formula '{formula}'"
                for field in ("fraction", "partition"):
                    val = meta.get(field)
                    if isinstance(val, str) and ph_re.search(val):
                        return True, f"species '{sp_name}' {field} expression '{val}'"
                    elif isinstance(val, dict):
                        for pk, pv in val.items():
                            if isinstance(pv, str) and ph_re.search(pv):
                                return True, f"species '{sp_name}' {field}.{pk} expression '{pv}'"

        # 2. Speciation definitions
        if isinstance(speciation_defs, dict):
            for sys_name, sys_cfg in speciation_defs.items():
                if not isinstance(sys_cfg, dict):
                    continue
                for sub_k, sub_v in sys_cfg.get("species", {}).items():
                    if isinstance(sub_v, dict) and isinstance(sub_v.get("fraction"), str):
                        if ph_re.search(sub_v["fraction"]):
                            return True, f"speciation '{sys_name}.{sub_k}' fraction expression"
                for f_k, f_expr in sys_cfg.get("fractions", {}).items():
                    if isinstance(f_expr, str) and ph_re.search(f_expr):
                        return True, f"speciation '{sys_name}' partition fraction '{f_k}'"

        # 3. Kinetic constants file
        if hasattr(self, "constants_file") and self.constants_file and self.constants_file.exists():
            try:
                with open(self.constants_file, "r") as f:
                    const_tree = ast.parse(f.read())
                for node in const_tree.body:
                    if isinstance(node, ast.FunctionDef) and node.name == "get_reaction_constants":
                        for stmt in node.body:
                            for sub in ast.walk(stmt):
                                if isinstance(sub, ast.Name) and isinstance(sub.ctx, ast.Load) and sub.id in ("pH", "Hplus"):
                                    return True, f"get_reaction_constants uses variable '{sub.id}'"
                                if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
                                    if ph_re.search(sub.value):
                                        return True, f"get_reaction_constants references '{sub.value}'"
            except Exception:
                pass

        # 4. Reactions file
        if hasattr(self, "reactions_file") and self.reactions_file and self.reactions_file.exists():
            try:
                with open(self.reactions_file, "r") as f:
                    rxn_tree = ast.parse(f.read())
                for node in rxn_tree.body:
                    if isinstance(node, ast.Assign):
                        for target in node.targets:
                            if isinstance(target, ast.Name) and target.id in ("reactions", "limiters"):
                                for sub in ast.walk(node.value):
                                    if isinstance(sub, ast.Name) and sub.id in ("pH", "Hplus"):
                                        return True, f"reactions file references '{sub.id}'"
                                    if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
                                        if ph_re.search(sub.value):
                                            return True, f"reactions file references '{sub.value}'"
            except Exception:
                pass

        return False, "no reactions or species depend on pH or Hplus"

    def classify_jacobian_strategy(self) -> Tuple[bool, str]:
        """Classifies whether the model requires the vectorized numerical Jacobian
        or can safely and optimally use the analytical SymPy Jacobian.

        Classification Heuristic:
          - Tier 1 (Static / Uncoupled Speciation, e.g. Velde 2016 with fixed pH):
            All speciation fractions are constant scalars or depend only on fixed constants.
            -> Algebraic SymPy Jacobian is mathematically exact, fastest, 0% FD error.
          - Tier 2 (Implicit / Iterative / Coupled Speciation, e.g. SimpleFe, dynamic carbonate):
            Speciation fractions dynamically couple multiple primary state variables,
            or dynamic carbonate equilibrium is coupled to mineral saturation.
            -> Vectorized Numerical Jacobian provides exact total derivatives, preventing
               Newton divergence, dt cutbacks, and solver failure.

        Returns:
            (requires_numerical: bool, rationale: str)
        """
        species_defs = getattr(self.mp, "species_defs", {})
        speciation_defs = getattr(self.mp, "speciation", {})
        transport_species = {
            name for name, meta in species_defs.items()
            if isinstance(meta, dict) and meta.get("include", True)
        }

        # 1. Procedural speciation hooks (e.g. Du 2023 iron speciation)
        if self.p_dict.get("enable_iron_speciation", False) or "beta_FeHS" in self.p_dict or hasattr(self.mp, "beta_FeHS"):
            return True, "coupled non-linear iron speciation (enable_iron_speciation/beta_FeHS)"

        # 2. Check speciation systems for cross-species dependencies
        for sys_name, sys_cfg in speciation_defs.items():
            primary = sys_cfg.get("primary") or sys_cfg.get("parent")
            parent_names = set()
            if isinstance(primary, str):
                parent_names.add(primary)
            elif isinstance(primary, dict):
                parent_names.update(primary.values())
            elif isinstance(primary, (list, tuple)):
                parent_names.update(primary)

            # Check sub-species fractions
            sub_species = sys_cfg.get("species", {})
            if isinstance(sub_species, dict):
                for sp_k, sp_v in sub_species.items():
                    if isinstance(sp_v, dict):
                        frac_expr = sp_v.get("fraction", "")
                        if isinstance(frac_expr, str):
                            for tsp in transport_species:
                                if tsp not in parent_names and re.search(r'\b' + re.escape(tsp) + r'\b', frac_expr):
                                    return True, f"speciation '{sys_name}.{sp_k}' dynamically couples transport species '{tsp}'"

            # Check phase partition fractions
            fractions = sys_cfg.get("fractions", {})
            if isinstance(fractions, dict):
                for f_k, f_expr in fractions.items():
                    if isinstance(f_expr, str):
                        for tsp in transport_species:
                            if tsp not in parent_names and re.search(r'\b' + re.escape(tsp) + r'\b', f_expr):
                                return True, f"phase partition '{sys_name}.{f_k}' dynamically couples transport species '{tsp}'"

        # 3. Dynamic carbonate chemistry coupling to mineral saturation
        is_carbonate = getattr(self.mp, "is_carbonate", False) or "carbonate_equilibrium_chemistry" in self.p_dict
        has_ta_dic = "TA" in transport_species and "DIC" in transport_species
        if is_carbonate and has_ta_dic:
            has_non_carb_ta = False
            non_carb_species = []
            for sp_k, sp_v in species_defs.items():
                if isinstance(sp_v, dict) and sp_v.get("include_in_TA") is not None:
                    has_non_carb_ta = True
                    non_carb_species.append(sp_k)

            # Check if reactions contain mineral saturation
            has_saturation = False
            sat_rxn_name = ""
            for r in getattr(self, "_raw_reactions", []):
                rtype = r.get("type", "")
                rname = r.get("reaction_name", "")
                if rtype == "reversible_saturation" or "saturation" in rname.lower() or "caco3" in rname.lower() or "fes" in rname.lower():
                    has_saturation = True
                    sat_rxn_name = rname
                    break

            if has_non_carb_ta and has_saturation:
                return True, f"non-carbonate alkalinity ({', '.join(non_carb_species)}) coupled to carbonate equilibrium and mineral saturation ({sat_rxn_name})"
            if has_saturation:
                return True, f"dynamic carbonate equilibrium coupled to mineral saturation ({sat_rxn_name})"

        # 4. Otherwise, all speciation fractions are static scalars (fixed pH, constant Kd)
        return False, "all speciation fractions are static scalars / uncoupled (analytical SymPy Jacobian is mathematically exact)"

    def setup(self, plot_queue=None) -> Tuple[data_container, data_container, data_container, List[str], np.ndarray, data_container, Any]:
        """Assembles grid, species variables, transport coefficients, and compiles reaction equations."""
        self.validate_simulation_files()
        self.mp.update(self.p_dict)
        self.mp._caller_dir = str(self.base_dir)

        # Fallback defaults for geometry if not specified by user configuration
        if not hasattr(self.mp, "grid_setup") or self.mp.grid_setup is None:
            self.mp.grid_setup = [
                [0.0, 0.5, 0.0001, 0.01],
            ]

        # Fallback defaults for environmental and transport parameters
        if not hasattr(self.mp, "temp"):
            self.mp.temp = [10.0, 10.1]
        if not hasattr(self.mp, "w"):
            self.mp.w = self.Q_("0.2 cm/yr").to("m/s").magnitude
        if not hasattr(self.mp, "advection"):
            self.mp.advection = 0.0
        if not hasattr(self.mp, "phi"):
            self.mp.phi = 0.8
        if not hasattr(self.mp, "BT0"):
            self.mp.BT0 = 0.0
        if not hasattr(self.mp, "BT_depth"):
            self.mp.BT_depth = self.Q_("7.6 cm").to("meter").magnitude
        if not hasattr(self.mp, "BT_attenuation"):
            self.mp.BT_attenuation = self.Q_("2 cm").to("meter").magnitude
        if not hasattr(self.mp, "BI0"):
            self.mp.BI0 = 0.0
        if not hasattr(self.mp, "BI_depth"):
            self.mp.BI_depth = 0.0
        if not hasattr(self.mp, "isotopes"):
            self.mp.isotopes = True

        # 1. Parse species definitions & isotope standards
        raw_species = self._load_species_definitions()
        raw_stds = getattr(self._species_mod, "isotope_standards", {})
        raw_speciation = getattr(self._species_mod, "speciation", None)

        from mantis.species_parser import parse_and_validate_species_dict, validate_and_normalize_speciation
        species_defs, self._parsed_bc_map, iso_stds = parse_and_validate_species_dict(
            raw_species, raw_stds, ureg=self.Q_
        )
        speciation_defs = validate_and_normalize_speciation(raw_speciation, species_defs, ureg=self.Q_)
        self.mp.species_defs = species_defs
        self.mp.speciation = speciation_defs
        self.mp.isotope_standards = iso_stds

        # Expose standards on mp
        for std_k, std_m in iso_stds.items():
            setattr(self.mp, f"{std_k}_ratio", float(std_m["ratio"]))
            if std_m.get("element") == "S":
                self.mp.VCDT = float(std_m["ratio"])
            elif std_m.get("element") == "O":
                self.mp.VSMOW = float(std_m["ratio"])

        # 2. Extract pH and check carbonate system boundary consistency
        carb_cfg = self.p_dict.get("carbonate_equilibrium_chemistry")
        is_carbonate = carb_cfg is not None
        self.mp.is_carbonate = is_carbonate

        ph_val = None
        if "pH" in self.p_dict and self.p_dict["pH"] is not None:
            ph_val = float(self.p_dict["pH"])
        elif "Hplus" in species_defs and "bc_top" in species_defs["Hplus"]:
            h_bc = species_defs["Hplus"]["bc_top"]
            if isinstance(h_bc, dict) and "Dirichlet" in h_bc and h_bc["Dirichlet"] is not None:
                ph_val = float(h_bc["Dirichlet"])

        if is_carbonate:
            from mantis.seawater_carbonate import (
                calc_carbonate_constants,
                adjust_for_calcite_equilibrium,
                solve_carbonate_speciation_from_ta_dic,
                solve_carbonate_speciation_from_dic_ph,
                solve_carbonate_speciation_from_ta_ph,
            )

            # 1. Environmental parameters (P, S, T)
            pressure_bar = float(carb_cfg.get("pressure", 1.0))
            salinity = float(carb_cfg.get("salinity", 34.0))
            adjust_omega = bool(carb_cfg.get("adjust_omega", True))
            ph_scale = str(carb_cfg.get("ph_scale", carb_cfg.get("pH_scale", "free"))).lower().strip()

            temp_val = self.mp.temp
            temp_c = float(temp_val[0] if isinstance(temp_val, (list, tuple, np.ndarray)) else temp_val)

            carb_constants = calc_carbonate_constants(
                temp_c=temp_c, salinity=salinity, pressure_bar=pressure_bar, ph_scale=ph_scale
            )
            self.mp.carb_constants = carb_constants
            self.mp.ph_scale = carb_constants.ph_scale
            self.mp.scale_factor = carb_constants.scale_factor
            self.mp.pK1_carb = carb_constants.pK1
            self.mp.pK2_carb = carb_constants.pK2
            self.mp.K1_carb = carb_constants.K1
            self.mp.K2_carb = carb_constants.K2
            self.mp.Ksp_calcite = carb_constants.Ksp_calcite

            # 2. Check Ca2: must be defined and > 0
            ca2_top = None
            if "Ca2" in carb_cfg and carb_cfg["Ca2"] is not None:
                ca2_top = carb_cfg["Ca2"]
            elif "Ca2" in self.mp.get("bc_overrides", {}):
                ca2_top = self.mp["bc_overrides"]["Ca2"]
            elif hasattr(self, "_parsed_bc_map") and "Ca2" in self._parsed_bc_map:
                ca2_top = self._parsed_bc_map["Ca2"].get("top")
            elif "Ca2" in species_defs:
                ca2_meta = species_defs["Ca2"].get("bc_top", species_defs["Ca2"].get("bc", {}))
                if isinstance(ca2_meta, dict):
                    ca2_top = ca2_meta.get("Dirichlet", ca2_meta.get("top"))
                else:
                    ca2_top = ca2_meta

            if isinstance(ca2_top, dict):
                ca2_top = ca2_top.get("top")
            if isinstance(ca2_top, str):
                ca2_top = self.Q_(ca2_top).to("mol/m^3").magnitude

            if ca2_top is None or float(ca2_top) <= 0.0:
                print(
                    "\n[MANTIS ERROR] Carbonate equilibrium chemistry requires species 'Ca2' "
                    "with a valid, positive concentration at the upper boundary to compute calcite saturation state Omega.\n"
                    "Please define 'Ca2' in species.py or specify it in carbonate_equilibrium_chemistry.\n",
                    file=sys.stderr,
                )
                sys.exit(1)
            ca2_top = float(ca2_top)

            # 3. Extract TA, DIC, pH solely from carbonate_equilibrium_chemistry
            # When carbonate_equilibrium_chemistry is provided, species.py boundary values
            # for carbonate species are disregarded and decisions are based solely on carb_cfg.
            ta_top = carb_cfg.get("TA")
            if isinstance(ta_top, dict):
                ta_top = ta_top.get("top")
            if isinstance(ta_top, str):
                ta_top = self.Q_(ta_top).to("mol/m^3").magnitude

            dic_top = carb_cfg.get("DIC")
            if isinstance(dic_top, dict):
                dic_top = dic_top.get("top")
            if isinstance(dic_top, str):
                dic_top = self.Q_(dic_top).to("mol/m^3").magnitude

            ph_top = carb_cfg.get("pH")
            if ph_top is not None:
                ph_top = float(ph_top)

            has_ph = ph_top is not None
            has_ta = ta_top is not None
            has_dic = dic_top is not None

            if has_ph and has_ta and has_dic:
                print(
                    "\n[MANTIS ERROR] Overdefined carbonate system in carbonate_equilibrium_chemistry.\n"
                    "You have specified boundary values for all three: pH, DIC, and TA.\n"
                    "Because pH, DIC, and TA are thermodynamically linked in seawater, only two of the three\n"
                    "can be independently specified as boundary controls.\n"
                    "Please leave one of them unspecified (or None) in carbonate_equilibrium_chemistry.\n",
                    file=sys.stderr,
                )
                sys.exit(1)

            num_specified = sum([has_ph, has_ta, has_dic])
            if num_specified < 2:
                print(
                    "\n[MANTIS ERROR] Underdefined carbonate system in carbonate_equilibrium_chemistry.\n"
                    "At least two of (pH, DIC, TA) must be specified.\n",
                    file=sys.stderr,
                )
                sys.exit(1)

            if adjust_omega:
                carb_spec, adjust_msg = adjust_for_calcite_equilibrium(
                    ca2_top=ca2_top,
                    dic=dic_top,
                    ta=ta_top,
                    ph=ph_top,
                    const=carb_constants,
                )
                self.mp.pH = carb_spec["pH"]
                self.mp.carb_spec_top = carb_spec
                print(
                    f"\n[MANTIS CARBONATE] Carbonate equilibrium chemistry enabled:\n"
                    f"         Temperature: {temp_c:.1f} °C, Salinity: {salinity:.1f} PSU, Pressure: {pressure_bar:.1f} bar\n"
                    f"         pH scale: {carb_constants.ph_scale.capitalize()} (scale factor = {carb_constants.scale_factor:.4f})\n"
                    f"         Thermodynamic constants (Lueker et al. 2000, Mucci 1983, Millero 1995, Dickson 1990):\n"
                    f"           pK1 = {carb_constants.pK1:.3f}, pK2 = {carb_constants.pK2:.3f}, K'sp(calcite) = {carb_constants.Ksp_calcite:.4f} (mol/m^3)^2\n"
                    f"         adjust_omega = True:\n"
                    f"           {adjust_msg}\n"
                )
            else:
                if has_ta and has_dic:
                    carb_spec = solve_carbonate_speciation_from_ta_dic(ta_top, dic_top, carb_constants)
                elif has_ph and has_dic:
                    carb_spec = solve_carbonate_speciation_from_dic_ph(dic_top, ph_top, carb_constants)
                else:  # has_ph and has_ta
                    carb_spec = solve_carbonate_speciation_from_ta_ph(ta_top, ph_top, carb_constants)

                omega = (ca2_top * carb_spec["CO3"]) / carb_constants.Ksp_calcite
                self.mp.pH = carb_spec["pH"]
                self.mp.carb_spec_top = carb_spec
                omega_state = "supersaturated" if omega > 1.05 else ("undersaturated" if omega < 0.95 else "near equilibrium")
                print(
                    f"\n[MANTIS CARBONATE] Carbonate equilibrium chemistry enabled:\n"
                    f"         Temperature: {temp_c:.1f} °C, Salinity: {salinity:.1f} PSU, Pressure: {pressure_bar:.1f} bar\n"
                    f"         pH scale: {carb_constants.ph_scale.capitalize()} (scale factor = {carb_constants.scale_factor:.4f})\n"
                    f"         Thermodynamic constants (Lueker et al. 2000, Mucci 1983, Millero 1995, Dickson 1990):\n"
                    f"           pK1 = {carb_constants.pK1:.3f}, pK2 = {carb_constants.pK2:.3f}, K'sp(calcite) = {carb_constants.Ksp_calcite:.4f} (mol/m^3)^2\n"
                    f"         adjust_omega = False: Using boundary conditions as is.\n"
                    f"           Overlying water calcite saturation state: Omega = {omega:.2f} ({omega_state}).\n"
                    f"           Top boundary speciation: DIC = {carb_spec['DIC']:.2f}, TA = {carb_spec['TA']:.2f} mol/m^3, pH = {carb_spec['pH']:.2f}.\n"
                )
        else:
            if ph_val is not None:
                self.mp.pH = ph_val
            else:
                needs_ph, reason = self.check_requires_ph(species_defs, speciation_defs)
                if needs_ph:
                    print(
                        f"\n[MANTIS ERROR] pH is not defined in the model, but is required ({reason}).\n"
                        "Please specify pH in species.py (e.g. Hplus: {'bc': {'pH': 7.5}}) "
                        "or pass model.configure(pH=7.5).\n",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                else:
                    self.mp.pH = None

        # 3. Compile or load reaction equations
        if self.is_diffusion_only:
            from mantis.generate_equation_system import make_null_reactions
            eq_mod = make_null_reactions()
            self.mp.reactions_module = eq_mod
            self._raw_reactions = []
            get_rc_fn = None
            self.p_dict["use_numerical_jacobian"] = False
            self.mp.use_numerical_jacobian = False
        else:
            prep_dict = dict(self.mp)
            prep_dict["reactions"] = str(self.reactions_file) if self.reactions_file else None
            prep_dict["species"] = str(self.species_file)
            prep_dict["kinetic_constants"] = str(self.constants_file) if self.constants_file else None
            if self.limiters_file and self.limiters_file.exists():
                prep_dict["limiters"] = str(self.limiters_file)
            if self.equations_file:
                prep_dict["equations"] = str(self.equations_file)

            eq_mod, _rxn_p, _eq_p = resolve_and_prepare_equations(
                prep_dict, experiment_name=self.experiment, base_dir=self.base_dir
            )
            self.mp.reactions_module = eq_mod
            if hasattr(eq_mod, "get_limiters"):
                self.mp.get_limiters = eq_mod.get_limiters

            # Extract reaction_parameters, fractionation_factors, and raw reaction definitions
            rxn_params = getattr(eq_mod, "reaction_parameters", None) or getattr(eq_mod, "fractionation_factors", None)
            if hasattr(self, "reactions_file") and self.reactions_file and self.reactions_file.exists():
                spec_rxn = importlib.util.spec_from_file_location("rxn_param_mod", self.reactions_file)
                if spec_rxn and spec_rxn.loader:
                    mod_rxn = importlib.util.module_from_spec(spec_rxn)
                    spec_rxn.loader.exec_module(mod_rxn)
                    self._raw_reactions = getattr(mod_rxn, "reactions", [])
                    if rxn_params is None:
                        rxn_params = getattr(mod_rxn, "reaction_parameters", None) or getattr(mod_rxn, "fractionation_factors", None)
            else:
                self._raw_reactions = []

            if isinstance(rxn_params, dict):
                for rp_k, rp_v in rxn_params.items():
                    if not hasattr(self.mp, rp_k):
                        setattr(self.mp, rp_k, rp_v)

            # 4. Load reaction constants
            spec_rc = importlib.util.spec_from_file_location("rc_mod", self.constants_file)
            if spec_rc is None or spec_rc.loader is None:
                raise ImportError(f"Cannot load reaction constants file: {self.constants_file}")
            mod_rc = importlib.util.module_from_spec(spec_rc)
            spec_rc.loader.exec_module(mod_rc)

            if hasattr(mod_rc, "get_reaction_constants"):
                get_rc_fn = mod_rc.get_reaction_constants
            else:
                raise AttributeError(f"{self.constants_file} must define get_reaction_constants(pH, phi, k_values).")

            # 5. Chemical Jacobian strategy selection (heuristic switching between analytical and numerical)
            if "use_numerical_jacobian" not in getattr(self, "_user_configured_keys", set()):
                requires_num, reason = self.classify_jacobian_strategy()
                if requires_num:
                    self.p_dict["use_numerical_jacobian"] = True
                    self.mp.use_numerical_jacobian = True
                    print(
                        f"\n[MANTIS] Coupled non-linear speciation detected ({reason}).\n"
                        "         Automatically enabled vectorized finite-difference chemical Jacobian for exact cross-coupling.\n"
                    )
                else:
                    self.p_dict["use_numerical_jacobian"] = False
                    self.mp.use_numerical_jacobian = False

        # 6. Grid generation
        mesh, z = make_grid_from_sections(self.mp.grid_setup)
        self.mp.mesh = mesh
        self.mp.z = z
        self.mp.grid_points = len(z)
        self.mp.max_depth = float(mesh.faceCoordinates[-1])
        if not hasattr(self.mp, "display_length") or self.mp.display_length is None:
            self.mp.display_length = self.mp.max_depth
        if isinstance(self.mp.phi, dict):
            from mantis.diff_lib import compute_porosity_profile
            phi_top = float(self.mp.phi.get("top", 0.9))
            phi_bottom = float(self.mp.phi.get("bottom", self.mp.phi.get("inf", 0.8)))
            z_half = self.mp.phi.get("half_depth", self.mp.phi.get("z_half", None))
            z_decay = self.mp.phi.get("decay", self.mp.phi.get("z_decay", None))
            phi_arr = compute_porosity_profile(
                z,
                phi_top=phi_top,
                phi_bottom=phi_bottom,
                z_half=z_half,
                z_decay=z_decay,
            )
            self.mp.phi = CellVariable(name="porosity", mesh=mesh, value=phi_arr)
        elif callable(self.mp.phi):
            self.mp.phi = CellVariable(name="porosity", mesh=mesh, value=self.mp.phi(z))
        else:
            self.mp.phi = CellVariable(name="porosity", mesh=mesh, value=self.mp.phi)

        # 6. Evaluate k values
        k = data_container()
        if get_rc_fn is not None:
            _k1, k = get_rc_fn(self.mp.pH, self.mp.phi, k_values=k)
        self.mp.k = k

        # Propagate reaction and equilibrium constants from k to self.mp if not already configured
        k_items = k.items() if hasattr(k, "items") else k.__dict__.items()
        for k_name, k_val in k_items:
            if not hasattr(self.mp, k_name):
                setattr(self.mp, k_name, k_val)

        # 7. Species lists & phase classification
        species_transport: List[str] = []
        species_all: List[str] = []
        bc_map: Dict[str, Dict[str, Any]] = {}
        D_mol = data_container()
        c = data_container()

        # Dynamic speciation & phase partition calculations from species_defs
        Hplus = 10.0 ** (-float(self.mp.pH)) if self.mp.pH is not None else None
        eval_ctx: Dict[str, Any] = {
            "k": k,
            "mp": self.mp,
            "np": np,
            "Hplus": Hplus,
            "pH": float(self.mp.pH) if self.mp.pH is not None else None,
        }

        for sp_name, meta in species_defs.items():
            # Check partition for mixed species (e.g. Fe2_total: dissolved, sorbed)
            if "partition" in meta and isinstance(meta["partition"], dict):
                for p_key, p_expr in meta["partition"].items():
                    raw_val = eval(p_expr, eval_ctx)
                    val = float(raw_val) if np.ndim(raw_val) == 0 else np.asarray(raw_val, dtype=np.float64)
                    setattr(self.mp, f"{sp_name}_{p_key}", val)
                    eval_ctx[f"{sp_name}_{p_key}"] = val
                    # Alias for Fe2_total.dissolved / sorbed
                    if sp_name == "Fe2_total" and p_key == "dissolved":
                        self.mp.Fe2_diss = val
                        eval_ctx["Fe2_diss"] = val
                    elif sp_name == "Fe2_total" and p_key == "sorbed":
                        self.mp.Fe2_sorb = val
                        eval_ctx["Fe2_sorb"] = val

            # Check fraction for equilibrium speciation (e.g. HS: Ka1 / (Hplus + Ka1))
            if meta.get("fraction") is not None:
                frac_ctx = dict(eval_ctx)
                if "constants" in meta and isinstance(meta["constants"], dict):
                    for c_name, c_val in meta["constants"].items():
                        frac_ctx[c_name] = float(c_val)
                        if c_name.startswith("pK"):
                            base_c = c_name[1:]
                            frac_ctx[base_c] = 10.0 ** (-float(c_val))
                frac_raw = meta["fraction"]
                try:
                    if isinstance(frac_raw, (int, float)):
                        f_val = float(frac_raw)
                    else:
                        f_val = float(eval(str(frac_raw), frac_ctx))
                    setattr(self.mp, f"{sp_name.lower()}_frac", f_val)
                    eval_ctx[f"{sp_name.lower()}_frac"] = f_val
                except Exception:
                    # Dynamic fraction depends on spatially varying species not yet initialized;
                    # will be updated across grid by update_all_species during setup and solver steps.
                    pass

        # Register equilibrium isotope fractionation factors from speciation definitions
        if isinstance(speciation_defs, dict):
            for sys_name, sys_cfg in speciation_defs.items():
                eq_frac = sys_cfg.get("equilibrium_fractionation", sys_cfg.get("isotope_fractionation"))
                if eq_frac is not None:
                    alpha_entry = eq_frac.get("alpha") if isinstance(eq_frac, dict) else eq_frac
                    sub_sp = sys_cfg.get("species", {})
                    sp_keys = list(sub_sp.keys()) if isinstance(sub_sp, dict) else []

                    if isinstance(alpha_entry, dict):
                        for a_name, a_val in alpha_entry.items():
                            val = float(a_val)
                            clean_name = str(a_name).strip()
                            # Check for A_to_B or A/B
                            tokens = None
                            if "_to_" in clean_name:
                                tokens = clean_name.split("_to_")
                            elif " / " in clean_name:
                                tokens = clean_name.split(" / ")
                            elif "/" in clean_name:
                                tokens = clean_name.split("/")

                            if tokens and len(tokens) == 2:
                                s1, s2 = tokens[0].strip(), tokens[1].strip()
                                # Set specific key e.g. h2s_s32_to_hs_s32_alpha
                                setattr(self.mp, f"{s1.lower()}_to_{s2.lower()}_alpha", val)
                                setattr(self.mp, f"{s2.lower()}_to_{s1.lower()}_alpha", 1.0 / val)
                                # Also derive base species names (strip _32 or _S32)
                                s1_base = species_defs.get(s1, {}).get("isotope_of", s1)
                                s2_base = species_defs.get(s2, {}).get("isotope_of", s2)
                                setattr(self.mp, f"{s1_base.lower()}_{s2_base.lower()}_alpha", val)
                                setattr(self.mp, f"{s2_base.lower()}_{s1_base.lower()}_alpha", 1.0 / val)
                            else:
                                if not clean_name.endswith("_alpha"):
                                    clean_name = f"{clean_name}_alpha"
                                setattr(self.mp, clean_name.lower(), val)
                    elif isinstance(alpha_entry, (int, float)):
                        val = float(alpha_entry)
                        setattr(self.mp, f"{sys_name.lower()}_alpha", val)
                        if len(sp_keys) == 2:
                            setattr(self.mp, f"{sp_keys[1].lower()}_{sp_keys[0].lower()}_alpha", val)
                            setattr(self.mp, f"{sp_keys[0].lower()}_{sp_keys[1].lower()}_alpha", 1.0 / val)


        # Initialize CellVariables
        zeros = np.zeros(self.mp.grid_points)
        T_profile = np.linspace(self.mp.temp[0], self.mp.temp[1], self.mp.grid_points)

        # Sort species: included in transport vs diagnostic
        for sp_name, meta in species_defs.items():
            species_all.append(sp_name)
            setattr(D_mol, sp_name, zeros.copy())
            setattr(c, sp_name, CellVariable(name=sp_name, mesh=mesh, value=0.0, hasOld=True))

            is_inc = meta.get("include", True)
            is_iso = "isotope_of" in meta

            if is_inc:
                if not is_iso or self.mp.isotopes:
                    species_transport.append(sp_name)

        # Build boundary conditions from strictly parsed bc_map
        bc_map = copy.deepcopy(self._parsed_bc_map)
        dbl_val = float(self.p_dict.get("diffusive_boundary_layer", self.p_dict.get("dbl", 0.0)))
        if dbl_val > 0.0:
            for sp_name, b_entry in bc_map.items():
                if b_entry.get("phase") == "liquid" and "dbl" not in b_entry:
                    b_entry["dbl"] = dbl_val

        # Apply any explicit user overrides to boundary conditions
        bc_overrides = self.mp.get("bc_overrides", {})
        for sp_k, ov in bc_overrides.items():
            if sp_k in bc_map:
                if isinstance(ov, dict):
                    bc_map[sp_k].update(ov)
                else:
                    bc_map[sp_k]["top"] = float(ov)

        # Carbonate system boundary condition consistency:
        if is_carbonate and hasattr(self.mp, "carb_spec_top"):
            carb_spec = self.mp.carb_spec_top
            if "TA" in bc_map:
                bc_map["TA"]["top"] = carb_spec["TA"]
            if "DIC" in bc_map:
                bc_map["DIC"]["top"] = carb_spec["DIC"]
            if "CO3" in bc_map:
                bc_map["CO3"]["top"] = carb_spec["CO3"]
            if "HCO3" in bc_map:
                bc_map["HCO3"]["top"] = carb_spec["HCO3"]

        # 8. Compute molecular diffusion coefficients from species recipes
        # Import diff_lib module for diff_coeff functions
        import mantis.diff_lib as diff_lib_mod

        for sp_name, meta in species_defs.items():
            diff_spec = meta.get("diffusion")
            fn_name = meta.get("diff_coeff_function")

            # Check if this species inherits diffusion from another species
            if isinstance(diff_spec, str):
                if hasattr(D_mol, diff_spec):
                    setattr(D_mol, sp_name, getattr(D_mol, diff_spec))
                continue

            # Determine function name from diff_coeff_function or dict 'type'
            if fn_name is None and isinstance(diff_spec, dict):
                fn_name = diff_spec.get("type", "diff_coeff")

            if fn_name is not None:
                # 1. Look up function in species.py module first, then diff_lib
                diff_fn = None
                if hasattr(self, "_species_mod") and hasattr(self._species_mod, fn_name):
                    diff_fn = getattr(self._species_mod, fn_name)
                elif hasattr(diff_lib_mod, fn_name):
                    diff_fn = getattr(diff_lib_mod, fn_name)

                if diff_fn is not None:
                    sig = inspect.signature(diff_fn)
                    kwargs = {}
                    if "T" in sig.parameters:
                        kwargs["T"] = T_profile
                    if "T_profile" in sig.parameters:
                        kwargs["T_profile"] = T_profile
                    if "temp" in sig.parameters:
                        kwargs["temp"] = T_profile
                    if "phi" in sig.parameters:
                        kwargs["phi"] = self.mp.phi
                    if isinstance(diff_spec, dict):
                        if "m0" in sig.parameters and "m0" in diff_spec:
                            kwargs["m0"] = diff_spec["m0"]
                        if "m1" in sig.parameters and "m1" in diff_spec:
                            kwargs["m1"] = diff_spec["m1"]
                    try:
                        d_val = diff_fn(**kwargs)
                        setattr(D_mol, sp_name, d_val)
                    except Exception as e:
                        raise RuntimeError(f"Failed to compute diffusion for '{sp_name}' with {fn_name}(): {e}")
                else:
                    raise AttributeError(f"Diffusion function '{fn_name}' for species '{sp_name}' not found in species.py or diff_lib.")

        # Lumped mixed-phase diffusion scaling (porewater retardation)
        for sp_name, meta in species_defs.items():
            if meta.get("phase") == "mixed" and hasattr(D_mol, sp_name):
                diss_frac = 1.0
                part = meta.get("partition", {})
                if "dissolved" in part:
                    try:
                        raw_diss = eval(part["dissolved"], eval_ctx)
                        diss_frac = float(raw_diss) if np.ndim(raw_diss) == 0 else np.asarray(raw_diss, dtype=np.float64)
                    except Exception:
                        diss_frac = getattr(self.mp, f"{sp_name}_dissolved", getattr(self.mp, "Fe2_diss", 1.0))
                else:
                    diss_frac = getattr(self.mp, f"{sp_name}_dissolved", getattr(self.mp, "Fe2_diss", 1.0))
                setattr(D_mol, sp_name, getattr(D_mol, sp_name) * diss_frac)

        # Bioturbation & Bio-irrigation profiles
        bi_atten = getattr(self.mp, "BI_attenuation", self.mp.BI_depth)
        D_mol.D_irr = compute_bio_irrigation_alpha(z, self.mp.BI0, bi_atten)

        bt_profile = getattr(self.mp, "BT_profile", "sigmoidal")
        if bt_profile == "exponential":
            bt_atten = getattr(self.mp, "BT_attenuation", getattr(self.mp, "BT_depth", 0.03))
            D_mol.D_bio = compute_exponential_db(z, self.mp.BT0, bt_atten)
        else:
            D_mol.D_bio = compute_sigmoidal_db(z, self.mp.BT0, self.mp.BT_depth, self.mp.BT_attenuation)

        # 9. Register boundary values on mp
        for sp_k, b_entry in bc_map.items():
            setattr(self.mp, f"bc_{sp_k}", b_entry.get("top", 0.0))

        # 10. Initial conditions and boundary constraints
        phi_val = self.mp.phi.value if hasattr(self.mp.phi, "value") else self.mp.phi
        phi_top = phi_val[0] if hasattr(phi_val, "__getitem__") else phi_val
        phi_bottom = phi_val[-1] if hasattr(phi_val, "__getitem__") else phi_val
        w_top = self.mp.w
        if hasattr(phi_val, "__getitem__") and len(phi_val) > 1 and np.ndim(self.mp.w) == 0:
            w_top = float(self.mp.w) * (1.0 - phi_bottom) / np.maximum(1.0 - phi_top, 1e-6)
        for sp_name in species_transport:
            var = getattr(c, sp_name)
            props = bc_map.get(sp_name, {"top": 0.0, "type": "Dirichlet", "phase": "liquid"})
            meta = species_defs.get(sp_name, {})

            # Determine prefill flag (default: False)
            if "isotope_of" in meta:
                parent_name = meta["isotope_of"]
                parent_meta = species_defs.get(parent_name, {})
                prefill = meta.get("prefill", parent_meta.get("prefill", False))
            else:
                prefill = meta.get("prefill", False)

            bc_type = props.get("type", "Dirichlet")
            is_solid = (props.get("phase") == "solid")

            if bc_type == "Robin":
                eff_phi_top = (1.0 - phi_top) if is_solid else phi_top
                J_flux = props["top"] / np.maximum(eff_phi_top, 1e-6)
                val = J_flux / w_top if w_top > 0 else 0.0

                if "init" in meta and meta["init"] is not None:
                    init_val = meta["init"]
                    if isinstance(init_val, str):
                        init_val = self.Q_(init_val).to("mol/m^3").magnitude
                    var.setValue(init_val)
                elif prefill:
                    var.setValue(val)
                else:
                    var.setValue(0.0)

                D_total = getattr(D_mol, sp_name, 0.0) + (D_mol.D_bio if is_solid else 0.0)
                if not isinstance(D_total, CellVariable):
                    D_total = CellVariable(mesh=mesh, value=D_total)

                d_left = D_total.faceValue[mesh.facesLeft.value][0]
                if d_left > 1e-20:
                    var.faceGrad.constrain(
                        [(self.mp.w * var.faceValue - J_flux) / D_total.faceValue],
                        mesh.facesLeft,
                    )
                else:
                    var.constrain(val, mesh.facesLeft)
            elif bc_type == "Dirichlet":
                val = props["top"]
                if "init" in meta and meta["init"] is not None:
                    init_val = meta["init"]
                    if isinstance(init_val, str):
                        init_val = self.Q_(init_val).to("mol/m^3").magnitude
                    var.setValue(init_val)
                elif prefill:
                    var.setValue(val)
                else:
                    var.setValue(0.0)
                var.constrain(val, mesh.facesLeft)
            elif bc_type == "Neumann":
                grad = props["top"]
                if "init" in meta and meta["init"] is not None:
                    init_val = meta["init"]
                    if isinstance(init_val, str):
                        init_val = self.Q_(init_val).to("mol/m^3").magnitude
                    var.setValue(init_val)
                elif prefill:
                    var.setValue(0.0)
                else:
                    var.setValue(0.0)
                var.faceGrad.constrain([grad], mesh.facesLeft)

            bot_val = float(props.get("bot", 0.0))
            var.faceGrad.constrain([bot_val], mesh.facesRight)

        # Initialize carbonate speciation diagnostics across cells if carbonate system is active
        if getattr(self.mp, "is_carbonate", False) and hasattr(self.mp, "carb_spec_top"):
            from mantis.diff_lib import compute_carbonate_speciation

            carb_spec = self.mp.carb_spec_top
            if hasattr(c, "TA"):
                c.TA.setValue(carb_spec["TA"])
            if hasattr(c, "DIC"):
                c.DIC.setValue(carb_spec["DIC"])
            if hasattr(c, "CO3"):
                c.CO3.setValue(carb_spec["CO3"])
            if hasattr(c, "HCO3"):
                c.HCO3.setValue(carb_spec["HCO3"])
            if hasattr(c, "CO2"):
                c.CO2.setValue(carb_spec["CO2"])
            if hasattr(c, "pH"):
                c.pH.setValue(carb_spec["pH"])
            if hasattr(c, "Hplus"):
                c.Hplus.setValue(carb_spec["Hplus"] * 1000.0)

            compute_carbonate_speciation(c, self.mp)

        # Restore state if specified
        if self.mp.state_data:
            read_state(c, self.mp.state_data)

        check_peclet_numbers(mesh, self.mp, D_mol, species_transport, bc_map)

        self._c = c
        self._k = k
        self._species_transport = species_transport
        self._species_all = species_all
        self._D_mol = D_mol
        self._bc_map = bc_map
        self._mesh = mesh

        return self.mp, c, k, species_transport, z, D_mol, eq_mod.diagenetic_reactions

    def run(self, plot_queue=None, **override_kwargs) -> DiageneticModelResult:
        """Runs the simulation and returns a DiageneticModelResult object."""
        if override_kwargs:
            self.configure(**override_kwargs)

        self.setup(plot_queue=plot_queue)

        plotter = None
        video_dt = getattr(self.mp, "video_dt", None)
        if plot_queue is None and video_dt is not None:
            if not self.layout_file.exists():
                self._exit_with_missing_file_error(self._orig_layout, "plot layout file (required for video recording)")
            output_path = f"{self.mp.plot_name}.pdf"
            video_path = os.path.abspath(f"{self.mp.plot_name}.mp4")
            plotter = LivePlotter(
                layout_path=self.mp.layout_file,
                display_length=self.mp.display_length,
                output_path=output_path,
                video_path=video_path,
                gui=False,
                max_queue_size=getattr(self.mp, "video_queue_size", 100),
                resume_threshold=getattr(self.mp, "video_resume_threshold", 50),
                video_dpi=getattr(self.mp, "video_dpi", 120),
                video_workers=getattr(self.mp, "video_workers", None),
            )
            plotter.start()
            plot_queue = plotter.queue

        try:
            diag_rxn_func = self.mp.reactions_module.diagenetic_reactions

            def default_equilibrium_reactions(mp, c, k, f, RATES, dt):
                inst_rxns = getattr(mp, "instantenous_reactions", [])
                for r in inst_rxns:
                    sig = inspect.signature(r[0])
                    if "f" in sig.parameters:
                        r[0](c, r[1], mp, dt, RATES, f=f)
                    else:
                        r[0](c, r[1], mp, dt, RATES)
                return f, RATES

            step, max_change = run_solver(
                self.mp,
                self._c,
                self._species_all,
                self._species_transport,
                self._k,
                diag_rxn_func,
                default_equilibrium_reactions,
                self._mesh,
                self._D_mol,
                self._bc_map,
                self.mp.z,
                plot_queue=plot_queue,
            )
        finally:
            if plotter:
                plotter.stop()

        # Run post-simulation reactant coexistence diagnostic by default
        if getattr(self.mp, "check_reactant_coexistence", True):
            try:
                from .diff_lib import check_reactant_coexistence
                check_reactant_coexistence(
                    mesh_or_z=self._mesh if self._mesh is not None else self.mp.z,
                    c=self._c,
                    k=self._k,
                    D_mol=self._D_mol,
                    verbose=True,
                )
            except Exception as e:
                print(f"[Warning] Reactant coexistence diagnostic failed: {e}")

        max_steps = getattr(self.mp, "max_steps", None)
        converged = "Yes" if (max_steps is None or step < max_steps) else "No"
        total_time = getattr(self.mp, "current_time", getattr(self.mp, "t_end", 0.0))

        result = DiageneticModelResult(
            mp=self.mp,
            c=self._c,
            k=self._k,
            species_list=self._species_all if hasattr(self, "_species_all") else self._species_transport,
            z=self.mp.z,
            D_mol=self._D_mol,
            diagenetic_reactions=diag_rxn_func,
            converged=converged,
            step=step,
            total_time=total_time,
        )
        return result
