# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-10-19

### Added

- **Marine Carbonate Equilibrium Chemistry (`seawater_carbonate.py`)**:
  - Implemented consensus oceanographic equilibrium formulations (Lueker et al. 2000, Mucci 1983, Millero 1995, Ingle 1975, UNESCO 1981 density equation) to compute stoichiometric constants ($K_1^*, K_2^*, K_{\text{sp}}'$) and calcite saturation state ($\Omega$) as functions of in-situ temperature, salinity, and hydrostatic pressure.
  - Carbonate speciation routines calculating $[\text{CO}_2^*]$, $[\text{HCO}_3^-]$, and $[\text{CO}_3^{2-}]$ from any valid parameter pair: pH & DIC, pH & TA, or DIC & TA.
  - Boundary condition validation and automatic equilibrium adjustment (`adjust_omega=True`), preventing overdefinition of the carbonate system and initializing sediment porewater profiles in thermodynamic equilibrium with solid $\text{CaCO}_3$.
  - Gated all carbonate equilibrium speciation behind the `carbonate_equilibrium_chemistry` parameter dictionary in `model.configure()`, isolating non-carbonate models from side effects.
  - Added scientific documentation with thermodynamic derivations and literature citations in `docs/carbonate_chemistry.org`.
- **Species-Agnostic Equation System & SymPy Analytical Jacobian Generator**:
  - Completely generalized `generate_equation_system.py` to be species-agnostic, replacing hardcoded species checks (`HS`, `H2S`, `Fe2`, `FeS`, `TS2`, `VCDT`) with phase-permutation detection (`liquid_to_liquid`, `solid_to_liquid`, `liquid_to_solid`, `solid_transformation`) and reaction mechanism classification.
  - Automated analytical chemical Jacobian assembly via SymPy symbolic differentiation.
  - Dynamic derivation of speciation factors and equilibrium partitioning from `species.py` metadata (`fraction`, `equilibrium_of`, `alias_of`, `isotope_of`, `isotope_standards`).
  - Added generic helper functions `calculate_fractionated_rate_coeff` and `partition_equilibrium_isotope` in `diff_lib.py`.
- **Multi-Isotope Cross-Derivatives & Coupling**:
  - Populated analytical cross-derivatives for companion isotope species when differentiating with respect to non-isotope reactants, scaling rates by reactant isotope ratios to eliminate numerical front ratio divergence.
  - Generalized cross-derivative logic across all isotope systems (S, C, N, O), dynamically resolving standard ratios from species definitions.
  - Added `get_total_delta` function in `diff_lib.py` for multi-species pooled isotope delta calculations with test suite in `tests/test_get_total_delta.py`.
- **Calcite Reactions & Diagnostics**:
  - Added dedicated calcite precipitation and dissolution reaction handling with smooth transitions across saturation state ($\Omega$).
  - Coupled reversible saturation isotope reactions to Total Alkalinity (TA) and Dissolved Inorganic Carbon (DIC) in the analytical chemical Jacobian.
  - Added carbonate system diagnostics and live plotting support.
  - Added `prefill=True` keyword in `species.py` to allow initializing species profiles to boundary values instead of zeros.
- **Timestamp-Aware Automatic Equation Regeneration**:
  - The equation generator now inspects modification timestamps of all core Mantis source files to automatically rebuild equation systems when the library is updated.

### Fixed

- Corrected units in reversible saturation reactions: unified $[\text{H}^+]$ scaling in $\partial\Omega/\partial h$ for mineral reactions and ensured consistent units for $k_{H^+}$.
- Simplified dynamic $H^+$ parameter handling in generated equation systems.
- Added `__future__.annotations` and `Sequence` import in `diff_lib.py` to resolve Python 3.12 CI compatibility.
- Excluded non-carbonate TA from kinetic isotope rate expressions.

## [0.2.1] - 2026-10-15
- fixed github CI errors

## [0.2.0] - 2026-10-15

### Added

- **Generalized Isotope Modeling & Standards Framework**:
  - Dynamically handle all isotope companion species across `diff_lib.py` and `live_plot_lib.py` instead of relying on hardcoded mappings.
  - Added support for explicit isotope reference standards (e.g., `VCDT` for sulfur) with custom standard ratios and metadata-defined `delta_name`.
  - Added support for species without explicit stoichiometry (defaulting to 1.0) and maintained backwards compatibility for legacy `d_<base>` column names alongside `d34S_<species>`.
  - Updated isotope unit representations from per mil (`permil` / `‰`) to milli Urey (`mUr`) per IUPAC conventions across documentation, examples, and plot exports.
 
- **Equation System Generation & Dynamic Limiters**:
  - Dynamic extraction of Monod half-saturation constants from rate expressions with phase-aware (solid vs. fluid phase) scaling based on sediment porosity.
  - Implemented dynamic limiter inclusion based on actual reactant and inhibition usage, eliminating unused terms and redundant code in generated equation scripts.
  - Added `force_regen` and `rebuild_equations` options in `generate_equation_system` and added the `--force-regen` CLI argument to bypass cached equations and rebuild.
  - Added automatic equation regeneration checks comparing file modification timestamps of generator and reaction scripts.
- **Model Configuration & Video Controls**:
  - Added strict pH validation during `DiageneticModel` initialization with descriptive error guidance.
  - Enhanced `configure()` with generic Pint-based unit conversions for geometric, environmental, transport, and grid parameters.
- **CI/CD & Documentation**:
  - Added Python 3.14 to GitHub Actions test matrix across Linux environments.
  - Added status, PyPI version, Python versions (3.12–3.14), and AGPL license badges to `README.md`.
  - Expanded tutorial documentation (`docs/tutorial.md`) with isotope delta boundary condition examples, updated geochemical constants, and equation regeneration options.
  - Added comprehensive test coverage in `tests/test_diagenetic_model.py` for isotope standards, delta columns, and missing file diagnostics.

### Changed

- **Examples**:
  - Streamlined `examples/velde_et_al_2016/run_velde_fast.py` by removing unnecessary manual time stepping parameters (`dt_init`, `dt_min`, `dt_tolerance`).
  - Simplified limiter calculation logic in `velde_et_al_2016` by consolidating phase scaling into a unified block.
- **Plotting**:
  - Refactored `plot_data_new.py` to dynamically detect and plot isotope delta columns.

### Removed

- Removed obsolete `docs/tutorial.org` file.
- Removed deprecated default parameter assignments for `Fe2_diss`, `Fe2_sorb`, `hs_frac`, and `h2s_frac`.

## [0.1.2] - 2026-10-12

### Added

- **File Path Resolution & Diagnostic Error Handling**:
  - Enhanced file path resolution in `DiageneticModel` (`_resolve`), searching the current working directory, the model script directory (`base_dir`), and fallback naming patterns (such as `<experiment>_reactions.py`).
  - Implemented early file validation (`validate_simulation_files()`) during `DiageneticModel` initialization and setup to catch missing simulation files before simulation execution.
  - Added informative error diagnostics (`[MANTIS ERROR]`) detailing exact search locations and resolution guidance when species, reactions, constants, limiters, or plot layout files are missing.
  - Improved plot layout file resolution in `plot_data_new.load_layout_from_file` and `plot_summary` with clear error reporting and search path feedback.
  - Added unit test suite in `tests/test_diagenetic_model.py` testing missing file scenarios, exit behavior, and directory resolution from both CWD and model base directories.

## [0.1.1] - 2026-10-12

### Added

- **Steady-State Detection & Convergence**:
  - Implemented a grid-independent spatial $L_2$ steady-state convergence check in `run_non_steady_state_solver_monolithic`. Computes a volume-weighted root-mean-square (RMS) change across grid cells per species and checks against `steady_state_tolerance` (or falls back to `dt_tolerance` if set).
  - Added a minimum step requirement (step > 5) to prevent premature steady-state termination.
- **Model Parameters & Unit Conversions**:
  - Added support for grid spacing parameters with unit conversion in `DiageneticModel`: `initial_spacing`, `reaction_zone_spacing`, and `max_spacing`.
  - Added support for the `reaction_zone` tuple/list interval parameter with unit conversion.
  - Added support for `steady_state_tolerance` and `dt_tolerance` solver tolerance parameters.
  - Added `check_reactant_coexistence` boolean option and enabled the post-simulation reactant coexistence diagnostic by default.
- **Examples**:
  - Added fast and slow simulation setups for the Velde et al. (2016) sediment diagenesis benchmark:
    - `examples/velde_et_al_2016/run_velde_fast.py` & `kinetic_constants_fast.py`: fast configuration with pyrite precipitation.
    - `examples/velde_et_al_2016/run_velde_slow.py` & `kinetic_constants_slow.py`: slow configuration tracking sulfur isotopes.
- **CI/CD & Automation**:
  - Added GitHub Actions workflow `.github/workflows/tests.yml` to run automated test suites across Python 3.12, 3.13, and 3.14.
  - Added GitHub Actions workflow `.github/workflows/publish-pypi.yml` for automated PyPI package publishing upon GitHub releases.
- **Documentation**:
  - Added comprehensive documentation on physical units (`pint` integration) and non-uniform grid generation to `docs/tutorial.md` and `docs/tutorial.org`.

### Changed

- **Rate Laws & Equation System Generation**:
  - Simplified particulate organic carbon (POC) respiration rate equations by removing redundant coefficient calculations.
  - Simplified FeS precipitation and dissolution rate formulation by removing unnecessary Michaelis-Menten terms.
  - Adjusted default smoothing parameter `fes_smooth_epsilon` from `0.05` to `0.01` for sharper, smoother reaction transitions.
  - Dynamically extract and handle active POC species and kinetic parameter pairs instead of using hardcoded POC fast/slow loops.
- **Differential Operators Library (`diff_lib`)**:
  - Added type annotations across numerical and differential operator functions.
  - Improved formatting, line wrapping, indentation, and visual structure.
- **Project Packaging & Metadata**:
  - Updated package distribution name to `mantis-rtm` and enabled dynamic versioning using `setuptools_scm`.
  - Added `pint` as an explicit core dependency in `pyproject.toml` and `requirements.txt`.
  - Added AGPL-3.0-or-later license configuration and classifiers.
  - Configured project URLs (Homepage, Repository, Issues) and updated build ignore patterns in `.gitignore`.
- **Examples**:
  - Corrected FeS dissolution rate constant in `examples/velde_et_al_2016/kinetic_constants_slow.py`.
  - Adjusted broken axis limits in `examples/velde_et_al_2016/plot_layout.py` for clearer concentration profile visualization.

### Removed

- Removed obsolete `examples/velde_et_al_2016/run_velde_model.py` script.
- Removed legacy table of contents from tutorial documentation.

## [0.1.0] - 2026-10-11

- Initial tagged release with core MANTIS reactive transport solver and tutorial documentation.
