# User facing functions and methods

This guide provides a concise yet context-rich reference for the primary user-facing functions, methods, and parameters in Mantis. It bridges the gap between high-level walkthrough tutorials and the raw API index, explaining what each parameter controls, typical values, default behaviors, and when or why you should adjust it.


## Model Interface


### `model.DiageneticModel`

The entry class for setting up, compiling, and executing a diagenetic reactive-transport simulation.

```python
from mantis import DiageneticModel

model = DiageneticModel(
    species="species.py",
    reactions="reactions.py",
    constants="kinetic_constants.py",
    equations="equations.py",
    layout="plot_layout.py",
    experiment="my_run",
)
```

-   `species` (str | Path, default: `"species.py"`):
    -   **Context**: Defines the chemical species, phases (`"liquid"`, `"solid"`, `"mixed"`), molecular diffusion coefficients, boundary conditions, and isotope companions.
    -   **When to set**: Point to a custom filename or shared definition if organizing multi-experiment directories.

-   `reactions` (str | Path | None, default: `"reactions.py"`):
    -   **Context**: Declares primary chemical reaction stoichiometries, Monod-style substrate limiters, and inhibition terms.
    -   **When to set**: Set to `None` if running a conservative or purely physical transport benchmark without diagenetic reactions.

-   `constants` (str | Path | None, default: `"reaction_constants.py"`):
    -   **Context**: Defines the rate constants, half-saturation constants, and equilibrium coefficients via `get_reaction_constants(pH, phi)`.
    -   **When to set**: Use when switching between alternative rate parameter sets (e.g., fast vs. slow organic matter degradation).

-   `equations` (str | Path | None, default: `None`):
    -   **Context**: Name of the pre-compiled Python equation file assembled by Mantis from declarative reactions.
    -   **When to set**: Leave as `None` or `"equations.py"` for automatic code generation and updates. Mantis inspects timestamps and regenerates equations when reaction definitions change.

-   `limiters` (str | Path | None, default: `None`):
    -   **Context**: Custom limitation and inhibition functions if overriding the standard Monod/hyperbolic forms.
    -   **When to set**: Only needed when using non-standard or user-written rate formulation expressions.

-   `layout` (str | Path, default: `"plot_layout.py"`):
    -   **Context**: Declarative configuration for live simulation monitoring and post-run summary plots.
    -   **When to set**: Specify a different layout file if generating specialized multi-panel figures for different sub-systems.

-   `experiment` (str | None, default: `None`):
    -   **Context**: Output prefix for generated artifacts (`<experiment>.pdf`, `<experiment>.csv`, `<experiment>_state.npz`, `<experiment>.mp4`).
    -   **When to set**: Defaults to the caller script&rsquo;s stem name (e.g., `run_velde_fast`). Set explicitly when looping over parameter sweeps or batch runs to avoid overwriting outputs.


### `model.configure`

Configures physical dimensions, numerical mesh, environmental conditions, transport dynamics, time-stepping, and solver tolerances. All physical quantities can be specified as unit strings (e.g., `"10 cm"`, `"2 kyr"`, `"0.2 cm/yr"`); Mantis handles automatic conversion into SI units.

```python
model.configure(
    grid_setup=[
        ["0 cm", "0.1 cm", "0.1 mm", "0.05 mm"],
        ["0.1 cm", "4.0 cm", "0.05 mm", "0.05 mm"],
        ["4.0 cm", "50 cm", "0.05 mm", "1.0 cm"],
    ],
    t_end="2 kyr",
    dt_max="1 year",
    w="0.2 cm/yr",
    phi=0.8,
    temp=[10.0, 10.1],
    bioturbation={"BT0": "4 cm^2/year", "depth": "6 cm", "attenuation": "2 cm"},
)
```


#### Grid & Geometry

-   `grid_setup` (list of 4-element lists):
    -   **Format**: `[[start, stop, start_dz, stop_dz], ...]`
    -   **Context**: Defines piecewise grid domains with geometric or spline transitions.
    -   **When to set**: Critical for resolving thin, sharp reaction fronts (such as oxic/suboxic boundaries or the SMTZ) with sub-millimeter spacing while smoothly expanding cell size at depth to keep computation fast.
-   `max_depth` (str | float, default: `"1.0 m"`):
    -   **Context**: Total sediment depth modeled.
    -   **When to set**: Used when `grid_setup` is omitted to generate a standard default grid.
-   `display_length` (str | float, optional):
    -   **Context**: Depth window shown in summary plots and videos.
    -   **When to set**: Set when the computational column extends to 5 m to capture deep boundary effects, but the geochemical action of interest is concentrated in the top 20 cm.


#### Time & Stepping

-   `t_end` (str | float, default: `"1 kyr"`):
    -   **Context**: Total physical time to integrate.
    -   **When to set**: Ensure `t_end` exceeds the transit time ($L/w$) or diffusive equilibration timescale ($L^2 / D$) if seeking steady-state profiles.
-   `dt_init` (str | float, default: `"1 day"`):
    -   **Context**: Initial time step size.
    -   **When to set**: Decrease if the simulation encounters convergence difficulties during early transients.
-   `dt_max` (str | float, default: `"1 year"`):
    -   **Context**: Maximum cap on time-step expansion.
    -   **When to set**: Set smaller for dynamic, rapidly changing systems; relax to larger values (e.g., `"10 year"`) when marching toward deep diagenetic steady state.
-   `dt_tolerance` / `steady_state_tolerance` (float, default: `1e-5`):
    -   **Context**: Relative tolerance threshold used to determine whether the column has reached steady state, triggering early completion.


#### Transport & Environment

-   `w` (str | float, default: `0.0`):
    -   **Context**: Sediment burial accumulation rate (e.g., `"0.2 cm/yr"`). Advects solids downward and drives porewater advection.
    -   **When to set**: Always set for accumulating marine or lacustrine sediment columns.
-   `advection` (str | float, default: `0.0`):
    -   **Context**: Additional porewater fluid advection velocity (e.g., `"10 cm/yr"` upwards in methane seep or cold vent systems).
    -   **When to set**: Essential for cold vents, groundwater discharge, or compaction-driven upward advection.
-   `phi` (float | dict | callable | list | ndarray, default: `0.8`):
    -   **Context**: Sediment porosity. May be a constant float, a `dict` defining an exponential compaction profile (keys: `"top"`, `"bottom"` / `"inf"`, and `"half_depth"` / `"z_half"` or `"decay"` / `"z_decay"`, which automatically calls [`compute_porosity_profile`](#org0ef8063)), a callable `phi(z)` accepting depth coordinates `z` [m], or a pre-computed `list` / `ndarray` matching the grid length. Note that passing a `list` or `ndarray` uses the values directly without calling `compute_porosity_profile`; Mantis provides [`compute_porosity_profile`](#org0ef8063) to compute exponential depth profiles directly.
    -   **When to set**: Set based on core measurements. Porosity scales solute inventories and governs tortuosity corrections ($D_{\text{eff}} = \phi^2 D_{\text{mol}}$).
-   `temp` (float | list, default: `10.0`):
    -   **Context**: Sediment temperature [°C]. Can be a constant or a two-element list `[T_top, T_bottom]` to impose a linear geothermal gradient.
    -   **When to set**: Affects temperature-dependent molecular diffusion coefficients.
-   `pH` (float, default: `7.5`):
    -   **Context**: Master sediment pH. Passed to `get_reaction_constants` and used in speciation/solubility calculations.
-   `diffusive_boundary_layer` / `dbl` (str | float, default: `0.0`):
    -   **Context**: Stagnant boundary layer thickness above the sediment-water interface (e.g., `"1 mm"`).
    -   **When to set**: Include when modeling micro-sensor benthic chamber experiments or low-energy deep-sea environments where DBL controls solute influx.
-   `bioturbation` (dict):
    -   **Keys**: `"BT0"` (biodiffusion coefficient, e.g., `"4 cm^2/yr"`), `"depth"` (mixed layer depth, e.g., `"6 cm"`), `"attenuation"` / `"decay"` (transition decay width, e.g., `"2 cm"`), `"profile"` (`"exponential"` or `"linear"`).
    -   **When to set**: Set whenever benthic macrofaunal particle mixing is present in upper sediment layers.
-   `bioirrigation` (dict):
    -   **Keys**: `"BI0"` (irrigation rate constant, e.g., `"182.6 /year"`), `"depth"`, `"attenuation"`, `"profile"`.
    -   **When to set**: Set when burrow-ventilating infauna pump bottom water into pore fluids, enhancing solute exchange independently of particle mixing.
-   `isotopes` (bool, default: `False`):
    -   **Context**: Enables automatic synthesis and tracking of isotopic companions and fractionation matrices.
    -   **When to set**: Set to `True` whenever modeling kinetic isotope systems (e.g., $^{34}\text{S}/^{32}\text{S}$, $^{13}\text{C}/^{12}\text{C}$).


#### Solver & Convergence

-   `solver` (str, default: `"monolithic"`):
    -   **Context**: Numerical solver engine. Mantis uses a monolithic coupled block Newton method.
-   `newton_tol` (float, default: `1e-3`):
    -   **Context**: Convergence tolerance for block Newton iterations.
-   `newton_norm` (str, default: `"linf"`):
    -   **Context**: Residual norm metric (`"linf"` for maximum cell error or `"wrms"` for weighted root-mean-square).
-   `use_numerical_jacobian` (bool, default: `False`):
    -   **Context**: Fallback to finite-difference numerical Jacobian if analytical SymPy Jacobian derivation is skipped or disabled.
-   `check_reactant_coexistence` (bool, default: `False`):
    -   **Context**: Diagnostic check verifying mutual exclusivity of incompatible chemical species (e.g., dissolved $\text{O}_2$ and $\text{H}_2\text{S}$).


#### Live Video Recording

-   `video_dt` (str | float | None, default: `None`):
    -   **Context**: Frequency at which frames are rendered to an MP4 video (e.g., `"10 yr"`). Set to `None` to disable.
-   `video_dpi` (int, default: `120`):
    -   **Context**: Video resolution DPI.
-   `video_workers` (int | None, optional):
    -   **Context**: Worker process pool size for background video frame encoding.


### `model.run`

Executes the simulation until `t_end` or steady-state convergence, returning a `DiageneticModelResult` object.

```python
result = model.run()
# Optional: pass keyword overrides directly
result = model.run(t_end="5 kyr", dt_max="0.5 year")
```

-   `plot_queue` (multiprocessing.Queue | None, optional):
    -   **Context**: Frame queue for live streaming during video generation. Handled automatically when `video_dt` is set.
-   `**override_kwargs`:
    -   **Context**: Any parameter accepted by `model.configure` can be passed directly to `model.run` as an override without modifying previous configuration calls.


### `result.save_csv`

Exports full depth-resolved numerical simulation results to a CSV table.

```python
df, fqfn = result.save_csv()
```

-   **Returns**:
    -   `df` (pandas.DataFrame): DataFrame containing depth (`z`, in meters), porosity (`phi`), temperature (`T`), concentrations (`c_<species>` in $\text{mol/m}^3$), bulk reaction rates (`k_<reaction>` in $\text{mol/(m}^3\cdot\text{s)}$), and calculated isotopic delta values (e.g., `d34S_SO4` in mUr/‰).
    -   `fqfn` (str): Fully qualified filesystem path to the saved CSV file (`<experiment>.csv`).
-   **When to use**: Call to export data for post-processing in external tools, plotting with Seaborn/Matplotlib, or archive generation.


### `result.save_state`

Saves the complete numerical state arrays to a compressed NumPy `.npz` archive.

```python
archive_path = result.save_state()  # saves to <experiment>_state.npz
# Or specify explicit path:
archive_path = result.save_state("checkpoints/spinup_state.npz")
```

-   `filepath` (str | Path | None, default: `None`):
    -   **Context**: Destination file path. If `None`, defaults to `<experiment>_state.npz` in the model script directory.
-   **When to use**: Essential for multi-stage spinups, restarting prolonged simulations, or performing scenario analyses where one steady-state provides the initial condition for transient pertubations.


### `result.save_report`


### `result.plot_summary()`

Renders a multi-panel publication-grade PDF plot based on the layout specification.

```python
pdf_path = result.plot_summary()
# Or override layout / output name:
pdf_path = result.plot_summary(layout_file="detailed_layout.py", outfile="final_figure.pdf")
```

-   `layout_file` (str | Path | None, default: `None`):
    -   **Context**: Layout script declaring the arrangement of panels, axis scalings, and units. If omitted, uses the layout passed to `DiageneticModel`.
-   `outfile` (str | None, default: `None`):
    -   **Context**: Destination PDF path. Defaults to `<experiment>.pdf`.
-   **When to use**: Call at the end of the driver script to immediately inspect the simulated profiles.


### `get_reaction_constants`

User-defined function placed in `kinetic_constants.py` that supplies rate constants and thresholds as unit-aware dictionaries.

```python
def get_reaction_constants(pH: float, phi: float) -> dict:
    return {
        "POC_fast": "10 1/year",
        "K_SO4": "0.9 umol/cm^3",
        "FeS_sp": "10**-3.5 mol/l",
        "Hplus": f"{10 ** (-pH)} mol/l",
    }
```

-   `pH` (float): Column pH passed from model configuration; allows dynamic calculation of $\text{H}^+$ concentration and pH-dependent mineral solubilities.
-   `phi` (float | ndarray): Sediment porosity; enables porosity-dependent rate adjustments.
-   **When to modify**: Adjust this file whenever calibrating reaction rates against empirical concentration profiles, testing parameter sensitivities, or altering chemical equilibrium constants.


### `get_total_delta`

Calculates bulk, phase-weighted isotopic delta values (e.g., total buried sulfur $\delta^{34}\text{S}$) across combined liquid and solid pools.

```python
from mantis.diff_lib import get_total_delta

d34S_burial = get_total_delta(
    result.c,
    result.mp,
    species="SO4, TS2, FeS",
    isotope_standard="S",
    index=-1,
)
print(f"Total buried d34S: {d34S_burial:.2f} mUr")
```

-   `c`: Simulation concentration state container (`result.c`).
-   `mp`: Model parameters container (`result.mp`).
-   `species` (str | Sequence[str], default: `"SO4, TS2, FeS"`): Comma-separated list or sequence of species to aggregate into the elemental isotope budget.
-   `isotope_standard` (str | float, default: look up in `mp.isotope_standards`): Element standard name (e.g., `"S"` for VCDT, `"C"` for VPDB) or direct reference ratio.
-   `index` (int | None, default: `-1`): Spatial index where the budget is computed. `-1` evaluates burial at the bottom boundary; `0` evaluates influx at the surface; `None` computes across the entire depth column.
-   **When to use**: Use to verify column-wide isotopic mass conservation or compute net isotopic signatures of buried sediments.


## Plot Interface

Mantis uses a declarative plot system configured in `plot_layout.py`. Individual panels are created using specification dictionaries or convenience wrapper classes (`liquid_species`, `solid_species`, `reaction_rates`, `isotopes`), which are then arranged in a 2D grid matrix via `plot_panels`.

```python
# Arranging panels: 2 rows of subplots
plot_panels = [
    [liquid_panel, solid_panel],   # Row 1: porewater solutes & solid minerals
    [rates_panel, isotopes_panel], # Row 2: reaction rates & isotope delta signatures
]
```

Global layout options include:

-   `colormap` (str, default: `"default"`): Color palette choice.
-   `colormap_override` (dict): Custom styling for individual species (e.g., `{"D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True}}`).
-   `z_scale` (str, default: `"lin"` or `"log"`): Vertical depth axis scale.
-   `z_lim` (tuple, default: `(0, 1)`): Depth axis limits `(z_min, z_max)` in meters.


### `liquid_species`

Configures subplots displaying dissolved porewater solute concentrations.

```python
panel_solutes = liquid_species({
    "species": ["SO4", "TS2", "O2", "Fe2_liq"],
    "units": "mmol/l",
    "scale": "lin",          # 'lin' or 'log'
    "lim": (0, 30),          # axis min and max in specified units
})
```

-   `species` (list of str): List of dissolved species names corresponding to entries in `species.py`.
-   `units` (str, default: internal $\text{mol/m}^3$): Display units (e.g., `"mmol/l"`, `"umol/l"`, `"uM"`, `"mM"`). Mantis automatically performs unit conversion.
-   `scale` (str, default: `"linear"`): Axis scale: `"lin"` or `"log"`. Logarithmic scaling is ideal when trace species (e.g., nanomolar oxygen) coexist on the same panel as millimolar sulfate.
-   `lim` (tuple, optional): Explicit limits `(min, max)` in designated units.
-   `broken_axis` / `lim` with 4 elements (tuple, optional): Generates a split axis (e.g., `(0, 2, 20, 30)`) to display low and high concentrations without wasting space on empty intervals.


### `solid_species`

Configures subplots displaying particulate and solid mineral phases.

```python
panel_minerals = solid_species([
    {
        "species": ["POC"],
        "units": ["wt"],
        "relative_to": "C",
        "scale": "lin",
        "lim": (0, 3),
    },
    {
        "species": ["FeS", "Fe3"],
        "units": ["wt"],
        "relative_to": "Fe",
        "scale": "lin",
        "lim": (0, 1),
    },
])
```

-   `species` (list of str): Solid species names (e.g., `["POC"]`, `["FeS", "Fe3"]`).
-   `units` (str | list of str): Display units:
    -   `"wt"` or `"percent"`: Weight percent of dry sediment (standard for sediment geochemical records).
    -   `"umol/g"`: Micromoles per gram dry sediment.
    -   `"mmol/l"`: Concentrations per volume of sediment.
-   `relative_to` (str, optional): Target element for weight percent conversion (e.g., `"C"` for POC wt% carbon; `"Fe"` for iron mineral wt% iron).
-   `scale` (str, default: `"linear"`): `"lin"` or `"log"`.
-   `lim` (tuple, optional): Axis limits `(min, max)`.


### `reaction_rates`

Configures subplots displaying rates of chemical transformations.

```python
panel_rates = reaction_rates({
    "species": ["POC", "SO4", "TS2", "FeS"],
    "units": "umol/(cm^3*yr)",  # optional display unit
    "scale": "symlog",           # standard for rates spanning orders of magnitude
    "lim": (-1e-2, 1e-2),
})
```

-   `species` (list of str): Chemical species whose net production or consumption rates should be plotted, or individual reaction names. You can also include physical transport profiles such as `"D_bio"` (biodiffusion) or `"bioirrigation"` to visualize biological mixing alongside reaction kinetics.
-   `units` (str, optional): Unit for reaction rates (e.g., `"umol/(cm^3*yr)"`, `"mol/(m^3*s)"`). Defaults to $\text{mol}&sdot;\text{m}<sup>-3</sup>&sdot;\text{s}<sup>-1</sup>$.
-   `scale` (str, default: `"symlog"`): Symmetrical logarithmic scaling, accommodating both positive (net production) and negative (net consumption) rates with a smooth linear region around zero.


### `isotopes`

Configures subplots displaying isotopic compositions ($&delta;$-values in mUr or ‰).

```python
panel_isotopes = isotopes({
    "species": ["d34S_SO4", "d34S_TS2", "d34S_FeS"],
    "scale": "lin",
    "lim": (-45, -35, 15, 25),  # 4-element tuple creates broken axis
})
```

-   `species` (list of str): Names of isotopic delta variables generated during simulation (e.g., `"d34S_SO4"`, `"d13C_DIC"`).
-   `lim` (tuple, optional):
    -   2 elements `(min, max)`: Standard continuous axis limits.
    -   4 elements `(min_lower, max_lower, min_upper, max_upper)`: Broken axis. In sulfur isotope systems where residual sulfate is enriched (+20‰) while biogenic sulfide and pyrite are depleted (-40‰), this cleanly excludes the large intermediate unpopulated zone.
-   Multiple isotope systems: If species belonging to multiple isotope systems (e.g., $\delta^{34}\text{S}$ and $\delta^{13}\text{C}$) are included in a single panel, Mantis automatically places the first group on the primary axis and the second group on a twin right axis.


## Utility Functions


### `compute_porosity_profile`

Computes an exponential sediment porosity profile with depth:

$$\phi(z) = \phi_{\text{bottom}} + (\phi_{\text{top}} - \phi_{\text{bottom}}) \exp\left(-\frac{z}{z_{\text{decay}}}\right)$$

```python
from mantis.diff_lib import compute_porosity_profile

phi_profile = compute_porosity_profile(
    z,
    phi_top=0.9,
    phi_bottom=0.75,
    z_half="5 cm",
)
```

-   `z` (array-like):
    -   **Context**: Sediment depth coordinates in meters.
-   `phi_top` (float):
    -   **Context**: Porosity at the sediment-water interface ($z = 0$).
-   `phi_bottom` (float):
    -   **Context**: Asymptotic porosity at depth ($z \to \infty$).
-   `z_half` (float | str | None, default: `None`):
    -   **Context**: Depth at which porosity has decreased halfway toward `phi_bottom` (accepts unit strings such as `"5 cm"` or meters as float). Related to `z_decay` by $z<sub>\text{decay}</sub> = z<sub>\text{half}</sub> / ln 2$.
-   `z_decay` (float | str | None, default: `None`):
    -   **Context**: $e$-folding attenuation depth ($1/\lambda$; accepts unit strings or meters as float). Either `z_half` or `z_decay` must be specified.
-   **Returns**:
    -   `np.ndarray`: Porosity profile evaluated at each depth coordinate in `z`.
-   **When to use**: Call directly when you want to construct, inspect, or customize an exponential porosity array before passing it to `model.configure(phi=...)`, or pass a configuration `dict` to `model.configure(phi={"top": 0.9, "bottom": 0.75, "half_depth": "5 cm"})` to have Mantis invoke it automatically.


### `compute_sediment_age`

Calculates the steady-state mean sediment age profile $A(z)$ [yr] by solving the 1D advective-diffusive aging equation for solid sediment particles:

$$\frac{d}{dz}\left[-(1 - \phi)\,D_b\,\frac{dA}{dz} + (1 - \phi)\,v_s\,A\right] = 1 - \phi$$

```python
from mantis.diff_lib import compute_sediment_age

age_yr = compute_sediment_age(
    z,
    w=0.002,       # m/yr (consistent with z in m)
    phi=phi_arr,
    Db=Db_arr,     # m^2/yr
)
```

-   `z` (array-like):
    -   **Context**: Sediment depth coordinates (e.g., in meters or centimeters). Units of `w` and `Db` must be consistent with `z`.
-   `w` (float | array-like):
    -   **Context**: Solid sediment burial velocity $v_s(z)$ in $\text{depth\_unit}/\text{yr}$.
-   `phi` (float | array-like | None, default: `None`):
    -   **Context**: Sediment porosity profile. Compaction is incorporated through the solid volume fraction $(1 - &phi;)$. If `None`, porosity is assumed uniform (no compaction).
-   `Db` (float | array-like | None, default: `None`):
    -   **Context**: Bioturbation biodiffusion coefficient profile in $\text{depth\_unit}<sup>2</sup>/\text{yr}$. If `None` or all zero, evaluates the pure unmixed advective burial integral. Bioirrigation affects only dissolved solutes and does not alter solid particle age.
-   **Returns**:
    -   `np.ndarray`: Mean sediment age profile $A(z)$ in years.
-   **When to use**: Use in post-processing or analysis scripts to convert sediment depth into mean depositional age, accounting for burial rate, compaction, and bioturbative mixing.
