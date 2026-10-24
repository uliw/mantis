# MANTIS quickstart guide 

## Overall Model Architecture: What is Where?

A typical MANTIS project directory contains five primary configuration files, alongside the master run script. An working example is provided in the ~./examples~ directory 

```
./velde_et_al_2016/
├── species.py                 <- Who lives in the sediment? (chemical species, phases, diffusion, BCs)
├── fes_reactions.py           <- What reactions happen? (stoichiometry, kinetics, limiters)
├── kinetic_constants_fast.py  <- How fast do reactions happen? (rate constants, half-saturation values)
├── plot_layout.py             <- How should the profiles look? (panels, colors, geochemist units like wt%)
└── run_velde_model_fast.py         <- Master driver script (sediment accumulation rate, depth, duration)
```

When you execute `python run_velde_model_fast.py`, MANTIS automatically detects these files, ensures the underlying equation system is up to date, solves the time evolution of the sediment column, and produces your outputs:

-   `<model_name>.pdf`: A multi-panel summary plot of porewater and solid-phase profiles.
-   `<model_name>.csv`: A tabular export of concentrations, reaction rates, and isotope ratios across depth.
-   `<model_name>_state.npz`: The complete numerical state of the simulation for restarting or continuing runs.
-   `<model_name>.mp4` (optional): An animated video showing the transient evolution


## Step 1: Defining Chemical Species (`species.py`)

The `species.py` file is your geochemical inventory. It declares every chemical entity that exists in your model domain.

### Understanding Phases: Dissolved vs. Solid vs. Mixed

In a sediment column, species exist in different physical compartments:

1.  **Liquid (Porewater)**: Solutes dissolved in the pore fluid (e.g., $\text{SO}_4^{2-}$, $\text{O}_2$, $\text{HS}^-$). Their concentrations are expressed per unit volume of porewater ($\text{mmol/L}_{\text{pw}}$) and their inventory is scaled by the sediment porosity ($\phi$). They move by molecular diffusion, porewater advection, and bio-irrigation.
2.  **Solid (Particulate)**: Minerals and particulate organic carbon bound in the sediment framework (e.g., $\text{FeS}_2$, reactive iron $\text{Fe}^{3+}$, organic carbon $\text{POC}$). Their inventory is scaled by the solid fraction ($1 - \phi$). They move solely via sediment accumulation (burial) and macrofaunal bioturbation.
3.  **Mixed Pools**: Lumped entities like total ferrous iron ($\text{Fe}^{2+}_{\text{total}}$) that partition rapidly into a dissolved fraction and a surface-adsorbed fraction.


### Example: Standard Species Configuration

Here is how common species are declared inside the `species` dictionary in `species.py`:

```python
# -----------------------------------------------------------------------------
# Isotope Reference Standards (milli Urey, mUr)
# -----------------------------------------------------------------------------
# Standard ratios must be explicitly defined under 'ratio'. MANTIS strictly looks
# up ratios from isotope_standards without arbitrary numerical fallbacks.
isotope_standards: dict = {
    "S": {
        "standard": "VCDT",
        "ratio": 0.044162589,  # 34S / 32S
        "symbol": "d34S",
        "unit": "mUr",
    },
    "C": {
        "standard": "VPDB",
        "ratio": 0.0112372,    # 13C / 12C
        "symbol": "d13C",
        "unit": "mUr",
    },
}

species = {
    # -------------------------------------------------------------------------
    # Dissolved Porewater Solutes
    # -------------------------------------------------------------------------
    "SO4": {
        "phase": "liquid",
        "include": True,               # True = solve transport PDE; False = diagnostic/fixed
        "formula": "SO4",
        "diff_coeff_function": "diff_coeff",
        "diffusion": {"m0": 4.88, "m1": 0.232},  # Boudreau (1996) diffusion parameters
        "bc": {"top": 28.2, "type": "dissolved", "delta": 21.0}, # 28.2 mmol/L seawater sulfate, d34S = +21 mUr
        "isotope": "SO4_32",           # Companion light-sulfur isotope species
        "isotope_standard": "S",       # Reference standard key in isotope_standards
        "delta_name": "d34S_SO4",      # Explicit output column name in CSV / df
        "prefill": True,               # True: prefill initial profile with top boundary concentration (28.2 mM)
    },

    "SO4_32": {
        "phase": "liquid",
        "include": True,
        "formula": "SO4_32",
        "diff_coeff_function": "diff_coeff",
        "diffusion": "SO4",            # Inherits the molecular diffusion of parent SO4
        "isotope_of": "SO4",           # Automatically inherits prefill=True from parent SO4
    },

    # -------------------------------------------------------------------------
    # Solid Minerals and Particulate Matter
    # -------------------------------------------------------------------------
    "Fe3": {
        "phase": "solid",
        "include": True,
        "formula": "Fe3",
        "diffusion": None,             # Solids do not have molecular diffusion!
        "bc": {"top": "12 umol/(cm^2 * yr)", "type": "particulate"}, # Solid depositional flux
        "prefill": False,              # False (default): initial column starts at 0.0 mmol/L
    },

    "FeS2": {
        "phase": "solid",
        "include": True,
        "formula": "FeS2",
        "diffusion": None,
        "bc": {"top": 0.0, "type": "particulate"},   # Formed authigenically in sediment
        "isotope": "FeS2_32",
        "stoichiometry": 2.0,          # 2 moles of S per mole of FeS2
        "isotope_standard": "S",
        "delta_name": "d34S_FeS2",
        "prefill": False,
    },

    "FeS2_32": {
        "phase": "solid",
        "include": True,
        "formula": "FeS2",
        "diffusion": None,
        "isotope_of": "FeS2",
    },

    # -------------------------------------------------------------------------
    # System Master Variables (e.g., pH)
    # -------------------------------------------------------------------------
    "Hplus": {
        "phase": "liquid",
        "include": False,              # Fixed background parameter, not an active transport PDE
        "formula": "H",
        "bc": {"pH": 7.5},             # Set porewater pH as the single source of truth!
        "prefill": True,
    },
}
```


### Initial Conditions and the `prefill` Keyword

When starting a simulation, each active transport variable must have an initial spatial concentration vector across the sediment grid:

- **`prefill: False` (Default)**:
  The concentration vector is initialized with zeros (`0.0`). This is typical for authigenic minerals (e.g. $\text{FeS}_2$, $\text{FeS}$), metabolic intermediates (e.g. $\text{HS}^-$, $\text{Fe}^{2+}$), and solid phases deposited incrementally at the sediment-water interface.
- **`prefill: True`**:
  The entire column is prefilled with the concentration at the upper boundary:
  - **For dissolved solutes (`type: "dissolved"`)**: The initial concentration is set directly to the top Dirichlet boundary value: $C(z, t=0) = C_{\text{top}}$ (e.g. $28.2\ \text{mmol/L}$ for seawater $\text{SO}_4^{2-}$).
  - **For particulate species (`type: "particulate"`)**: The upper boundary flux $J_{\text{bulk}}$ is converted to the solid-phase advective equilibrium concentration:
    $$C_{\text{top}} = \frac{J_{\text{solid}}}{w} = \frac{J_{\text{bulk}}}{(1 - \phi_{\text{top}}) \cdot w}$$
    If `prefill: True`, $C(z, t=0)$ is set to $C_{\text{top}}$.
- **Isotope Companions**:
  Isotope species (defined with `isotope_of`) automatically inherit the `prefill` setting of their parent species unless explicitly overridden. If `prefill: True`, the companion isotope vector is prefilled with the calculated isotopic mass abundance at the upper boundary.


### Carbonate System Modeling: Boundary Conditions & Porewater Equilibrium

Modeling the dissolved inorganic carbon (DIC) and total alkalinity (TA) system in porewaters involves coupled acid-base equilibria:
$$\text{CO}_2 + \text{H}_2\text{O} \rightleftharpoons \text{HCO}_3^- + \text{H}^+ \quad (K_1^*)$$
$$\text{HCO}_3^- \rightleftharpoons \text{CO}_3^{2-} + \text{H}^+ \quad (K_2^*)$$

#### 1. Transport Variables vs. Local Algebraic Diagnostics
In reactive transport modeling (e.g. Boudreau 1996, Soetaert et al. 1996), **$\text{TA}$ and $\text{DIC}$ are conservative transport scalars**—they satisfy standard advection-diffusion-reaction PDEs with stoichiometric reaction sources and sinks. In contrast, free proton concentration $[\text{H}^+]$ is not conservative because rapid acid-base buffering constantly exchanges protons among species. Therefore, porewater $\text{pH}$ and speciation are solved algebraically at each grid cell from the local $(\text{TA}, \text{DIC})$ balance.

#### 2. Upper Boundary Controls: Exactly Two of Three
Because porewater $\text{pH}$, $\text{DIC}$, and $\text{TA}$ are thermodynamically interdependent:
$$\text{DIC} = [\text{CO}_2] + [\text{HCO}_3^-] + [\text{CO}_3^{2-}]$$
$$\text{TA} = [\text{HCO}_3^-] + 2[\text{CO}_3^{2-}]$$
specifying all three independently at the upper boundary overdefines the chemical system.

MANTIS requires that **only two of the three** ($\text{pH}$, $\text{DIC}$, $\text{TA}$) have non-`None` boundary conditions. If all three are specified with numerical values, MANTIS halts execution with an explicit error prompting you to set `'top': None` (or `'pH': None`) for the dependent variable:
- **Case A: Define $(\text{pH}, \text{DIC})$** $\rightarrow$ Set `TA: {'bc': {'top': None, 'type': 'dissolved'}}`. MANTIS calculates the corresponding bottom water $\text{TA}_{\text{top}}$ and $[\text{CO}_3^{2-}]_{\text{top}}$.
- **Case B: Define $(\text{TA}, \text{DIC})$** $\rightarrow$ Set `Hplus: {'bc': {'pH': None}}`. MANTIS calculates bottom water $\text{pH}_{\text{top}}$ and $[\text{H}^+]$ from the alkalinity quadratic.
- **Case C: Define $(\text{pH}, \text{TA})$** $\rightarrow$ Set `DIC: {'bc': {'top': None, 'type': 'dissolved'}}`. MANTIS calculates bottom water $\text{DIC}_{\text{top}}$.

#### 3. Overlying Water vs. Sediment Porewater Initial Equilibrium
A critical distinction in marine diagenesis is between the **overlying water boundary condition** and the **sediment porewater initial condition**:
- **Overlying Bottom Water ($z = 0$)**: In the open ocean, bottom water is typically supersaturated with respect to calcite ($\Omega_{\text{calcite}} \sim 4\text{--}8$, $\text{pH} \approx 8.0\text{--}8.2$) because high dissolved $\text{Mg}^{2+}$ kinetically inhibits abiotic calcite precipitation in the water column.
- **Sediment Porewater ($z > 0$)**: Inside the sediment column, abundant mineral calcite surfaces ($\sim 50,000\ \text{mol/m}^3$ in carbonate-bearing sediments) rapidly buffer porewater to thermodynamic equilibrium ($\Omega_{\text{calcite}} = 1.0$, $\text{pH} \approx 7.22$ for typical seawater $\text{DIC} = 2.2\ \text{mM}$ and $\text{Ca}^{2+} = 10.3\ \text{mM}$).

If a simulation blindly copies the supersaturated bottom water ($\text{pH} = 8.2$) throughout the entire sediment column at $t = 0$ via `prefill: True`, every grid cell undergoes an immediate, artificial precipitation burst ($\sim 23,500\ \text{mol}/(\text{m}^3\cdot\text{yr})$). The adaptive Newton solver must then take small timesteps ($\Delta t \sim \text{hours}$) for thousands of iterations until the column relaxes to equilibrium.

To prevent this artificial transient, MANTIS automatically initializes the sediment porewater profile in thermodynamic equilibrium with solid $\text{CaCO}_3$:
$$[\text{CO}_3^{2-}]_{\text{init}} = \frac{K_{\text{sp, calcite}}}{[\text{Ca}^{2+}]}$$
$$\alpha_{2,\text{eq}} = \frac{[\text{CO}_3^{2-}]_{\text{init}}}{\text{DIC}_{\text{init}}} = \frac{K_1^* K_2^*}{[\text{H}^+]^2 + K_1^* [\text{H}^+] + K_1^* K_2^*}$$
$$[\text{TA}]_{\text{init}} = \text{DIC}_{\text{init}} \left( \frac{K_1^* [\text{H}^+] + 2 K_1^* K_2^*}{[\text{H}^+]^2 + K_1^* [\text{H}^+] + K_1^* K_2^*} \right)$$
while maintaining the user's specified overlying water concentration at $z = 0$. MANTIS announces this initialization with a clear console notice upon model startup:
```text
[MANTIS] Carbonate system detected (DIC, TA, Ca2, CaCO3):
         Sediment porewater initialized in equilibrium with solid CaCO3
         (Omega = 1.00, pH_init = 7.22, TA_init = 2.15 mol/m^3 based on DIC and Ca2+)
         to avoid artificial column-wide precipitation shock.
         Overlying water boundary condition remains TA_top = 2.40 mol/m^3 (pH_top = 7.96).
```
*(If you wish to override this equilibrium prefill with custom concentrations, you can explicitly set `'init': <value>` in `species.py`).*

#### 4. Vectorized Speciation & Dynamic pH Evolution
During the time-stepping loop, reactions like organic matter degradation, iron reduction, and calcite dissolution/precipitation produce or consume $\text{DIC}$ and $\text{TA}$. At each timestep, MANTIS solves the exact quadratic relationship for carbonate alkalinity in a fully vectorized manner across all grid cells via `compute_carbonate_speciation`, updating porewater $[\text{H}^+]$, dynamic $\text{pH}$, $[\text{CO}_3^{2-}]$, and $[\text{HCO}_3^-]$.

#### 5. Total Alkalinity vs. Carbonate Alkalinity: The `include_in_TA` Keyword

In marine geochemistry, **Total Alkalinity ($\text{TA}$)** is defined relative to reference protonation levels for conservative ions and weak acid-base pairs (**Dickson, 1981**):
$$\text{TA} = \underbrace{[\text{HCO}_3^-] + 2[\text{CO}_3^{2-}]}_{\text{Carbonate Alkalinity (CA)}} + [\text{B(OH)}_4^-] + [\text{OH}^-] + [\text{HPO}_4^{2-}] + 2[\text{PO}_4^{3-}] + [\text{HS}^-] + [\text{NH}_3] + \dots - [\text{H}^+] - [\text{HSO}_4^-] - [\text{HF}]$$

In MANTIS, **$\text{TA}$ is transported as the conservative master state variable**. However, carbonate mineral saturation ($\Omega_{\text{calcite}}$) and $\text{pH}$ depend specifically on **Carbonate Alkalinity ($\text{CA}$)**:
$$\text{CA} = \text{TA} - \text{TA}_{\text{non-carb}}$$
where $\text{TA}_{\text{non-carb}}$ represents the sum of all non-carbonate species that buffer protons:
$$\text{TA}_{\text{non-carb}} = \sum_i \nu_i [S_i]$$

To tell MANTIS which species contribute to alkalinity during speciation and $\text{pH}$ calculations, specify the **`include_in_TA`** keyword in `species.py`:

- **Positive Proton Acceptors ($\nu_i > 0$)**: Species that carry excess alkalinity relative to their zero-alkalinity reference level. For example, bisulfide ($\text{HS}^-$, $\nu = +1.0$), ammonia ($\text{NH}_3$, $\nu = +1.0$), borate ($\text{B(OH)}_4^-$, $\nu = +1.0$), or monohydrogen phosphate ($\text{HPO}_4^{2-}$, $\nu = +1.0$). When these species accumulate, they consume part of the total alkalinity pool, reducing the available carbonate alkalinity ($\text{CA} = \text{TA} - \sum \nu_i [S_i]$).
- **Negative Proton Donors ($\nu_i < 0$)**: Species that donate protons below their zero-alkalinity reference level. For example, excess hydronium ion ($[\text{H}^+]$, $\nu = -1.0$) or bisulfate ($\text{HSO}_4^-$, $\nu = -1.0$). Subtracting a negative term adds back to carbonate alkalinity.

##### Declarative Syntax in `species.py`

You can attach `include_in_TA` to fast equilibrium sub-species or directly transported species:

```python
species = {
    # 1. Equilibrium sub-species (HS- dynamically partitioned from TS2):
    "HS": {
        "phase": "liquid",
        "include": False,
        "formula": "HS",
        "equilibrium_of": "TS2",
        "fraction": "Ka1 / (Hplus + Ka1)",
        "constants": {"pKa1": 7.0},
        "include_in_TA": 1.0,          # Bisulfide contributes +1.0 * [HS-] to TA
    },

    # 2. Negative proton donor (e.g. bisulfate HSO4- partitioned from sulfate):
    "HSO4": {
        "phase": "liquid",
        "include": False,
        "formula": "HSO4",
        "equilibrium_of": "SO4",
        "fraction": "Hplus / (Hplus + Ka2_SO4)",
        "constants": {"pKa2_SO4": 1.5},
        "include_in_TA": -1.0,         # Bisulfate contributes -1.0 * [HSO4-] to TA
    },

    # 3. Transported weak base (e.g. unprotonated ammonia NH3):
    "NH3": {
        "phase": "liquid",
        "include": True,
        "formula": "NH3",
        "include_in_TA": 1.0,          # Contributes +1.0 * [NH3] to TA
        # ...
    },
}
```

##### Why Carbonate Species ($\text{HCO}_3^-$, $\text{CO}_3^{2-}$) Do NOT Take `include_in_TA`

A common question is: *if bisulfide ($\text{HS}^-$) takes `include_in_TA: 1.0`, why don't we specify `include_in_TA` for bicarbonate ($\text{HCO}_3^-$) or carbonate ($\text{CO}_3^{2-}$)?*

The answer lies in how the speciation engine calculates porewater $\text{pH}$:

1. **Carbonate Species Are Endogenous (The Solution, Not the Input)**:
   The carbonate alkalinity equation is:
   $$\text{CA} = [\text{HCO}_3^-] + 2[\text{CO}_3^{2-}]$$
   Using the thermodynamic equilibrium relations $[\text{HCO}_3^-] = \text{DIC} \frac{K_1^* [\text{H}^+]}{[\text{H}^+]^2 + K_1^* [\text{H}^+] + K_1^* K_2^*}$ and $[\text{CO}_3^{2-}] = \text{DIC} \frac{K_1^* K_2^*}{[\text{H}^+]^2 + K_1^* [\text{H}^+] + K_1^* K_2^*}$, MANTIS sets up the exact quadratic equation:
   $$\text{CA} [\text{H}^+]^2 + K_1^*(\text{CA} - \text{DIC}) [\text{H}^+] + K_1^* K_2^* (\text{CA} - 2\,\text{DIC}) = 0$$
   Notice that the $+1$ weighting of $\text{HCO}_3^-$ and the $+2$ weighting of $\text{CO}_3^{2-}$ are **already mathematically embedded in the coefficients of this quadratic equation**. They are the unknowns being solved for.

2. **What `include_in_TA` Actually Does**:
   The `include_in_TA` keyword exists exclusively for **non-carbonate** species that compete for the transported total alkalinity pool. MANTIS deducts non-carbonate alkalinity before solving the quadratic:
   $$\text{CA} = \text{TA} - \text{TA}_{\text{non-carb}} = \text{TA} - \sum \nu_i [S_i]$$
   If you were to attach `"include_in_TA": 1.0` to `HCO3` and `"include_in_TA": 2.0` to `CO3`, MANTIS would subtract carbonate alkalinity from itself:
   $$\text{CA} = \text{TA} - ([\text{HCO}_3^-] + 2[\text{CO}_3^{2-}]) = \text{TA} - \text{CA} = 0$$
   This would artificially collapse the carbonate alkalinity to zero and crash the speciation solver!

##### Why $\text{Ca}^{2+}$ and $\text{SO}_4^{2-}$ Do NOT Take `include_in_TA`

A frequent question is whether conservative ions like $\text{Ca}^{2+}$ and $\text{SO}_4^{2-}$ should be given an `include_in_TA` entry. **They should not**:

1. **Reaction Stoichiometry Covers Their Effect Completely**: The production or consumption of alkalinity associated with reactions involving $\text{Ca}^{2+}$ and $\text{SO}_4^{2-}$ is already written directly into the reaction stoichiometry in `reactions.py`:
   - *Sulfate Reduction*: $2\ \text{POC} + \text{SO}_4^{2-} \to 2\ \text{DIC} + \text{TS2} + 2\ \text{TA}$ (consuming 1 mole $\text{SO}_4^{2-}$ explicitly generates $+2\ \text{TA}$).
   - *Calcite Precipitation*: $\text{Ca}^{2+} + \text{DIC} + 2\ \text{TA} \to \text{CaCO}_3$ (consuming 1 mole $\text{Ca}^{2+}$ explicitly consumes $-2\ \text{TA}$).
2. **Conservative Reference Levels at Porewater pH**: In the neutral to mildly acidic pH range of marine porewaters ($\text{pH} \approx 6\text{--}9$), $\text{Ca}^{2+}$ is a conservative divalent cation with no acid-base equilibria, and sulfate is virtually $100\%$ deprotonated as $\text{SO}_4^{2-}$ (its second dissociation constant $pK_a \approx 1.0\text{--}1.5$ is well below ambient pH). Because $\text{SO}_4^{2-}$ is the zero-alkalinity reference species in Dickson's definition, unprotonated $\text{SO}_4^{2-}$ has zero contribution to speciation alkalinity. Adding `include_in_TA` to $\text{Ca}^{2+}$ or $\text{SO}_4^{2-}$ would double-count their effects and corrupt porewater $\text{pH}$ calculations.

##### Summary Rule of Thumb

| Chemical Category | Examples | `include_in_TA` in `species.py` | Mechanism |
| :--- | :--- | :--- | :--- |
| **Carbonate System** | $\text{HCO}_3^-$, $\text{CO}_3^{2-}$, $\text{CO}_2$ | **Omit** (Do NOT specify) | Endogenous to the carbonate quadratic engine; algebraic weights ($+1, +2$) are built into the solver. |
| **Non-Carbonate Weak Bases / Acids** | $\text{HS}^-$, $\text{NH}_3$, $\text{B(OH)}_4^-$, $\text{HPO}_4^{2-}$, $[\text{H}^+]$, $\text{HSO}_4^-$ | **Specify** (`1.0`, `-1.0`, etc.) | Deducted from $\text{TA}$ to yield the true carbonate alkalinity ($\text{CA} = \text{TA} - \sum \nu_i [S_i]$). |
| **Conservative Ions** | $\text{Ca}^{2+}$, $\text{SO}_4^{2-}$, $\text{Na}^+$, $\text{Cl}^-$ | **Omit** (Do NOT specify) | Zero-alkalinity reference species at porewater $\text{pH}$ (6–9). Their alkalinity effects are handled 100% by reaction stoichiometry in `reactions.py`. |


### Handling Diffusion Coefficients

For dissolved ions, molecular diffusion coefficients in sediment are calculated using the empirical linear temperature formulation from **Boudreau (1996)** corrected for tortuosity via Archie's law: $$ D_0(T) = (m_0 + m_1 \cdot T) \times 10^{-10} \quad [\text{m}^2/\text{s}] $$ $$ D_{\text{sed}} = \frac{D_0(T)}{1 - \ln(\phi^2)} $$

For ions following this relationship, you simply set:

```python
"diff_coeff_function": "diff_coeff",
"diffusion": {"m0": 4.88, "m1": 0.232},
```



#### Custom Diffusion Functions (e.g., Dissolved Oxygen)

Some species, such as dissolved $\text{O}_2$, require a non-standard formulation. You can declare your own custom diffusion function directly at the top of `species.py` and point to it by name:

```python
def boudreau_o2(T_profile, phi):
    """Temperature and Archie porosity corrected O2 molecular diffusion."""
    phi_val = phi.value if hasattr(phi, "value") else phi
    return (
        (0.2604 + 0.006363 * (T_profile + 273.15))
        * 1e-9
        / (1.0 - np.log(phi_val**2))
    )

species = {
    "O2": {
        "phase": "liquid",
        "include": True,
        "formula": "O2",
        "diff_coeff_function": "boudreau_o2",  # Calls the custom function defined above!
        "diffusion": None,
        "bc": {"top": 0.28, "type": "dissolved"}, # 280 uM typical bottom water O2
    },
    # ...
}
```

## Concentration Units & Conversions: Inputs, Internal Engine, and Outputs

A common source of confusion in diagenetic modeling is the distinction between porewater solute concentrations, solid sediment mineral concentrations, and geochemical laboratory units (such as $\text{wt}\%$ or $\mu\text{mol/g}$). MANTIS handles physical units transparently while maintaining strict mass conservation across phases.

### 1. Default Input Units

When defining your model setup and boundary conditions:

- **Driver Script Parameters (`run_velde_model_fast.py`)**:
  MANTIS supports human-readable strings parsed via Python's `pint` physical unit library:
  - Lengths and depths: `"0.5 m"`, `"50 cm"`
  - Times and time steps: `"0.5 kyr"`, `"1 year"`, `"24 hours"`
  - Sedimentation rates: `"0.2 cm/yr"`, `"1 mm/a"`
  - Bioturbation coefficients: `"4 cm^2/year"`, `"1e-9 m^2/s"`

- **Boundary Conditions (`species.py`)**:
  - **Dissolved Solutes (`type: "dissolved"`)**:
    Upper boundary concentrations can be specified as a numeric float or a unit string. The default numeric input is **$\text{mmol/L}$** (which is identical to $\text{mol/m}^3$ or $\text{mM}$).
    ```text
    "bc": {"top": 28.2, "type": "dissolved"}         # 28.2 mmol/L (seawater sulfate)
    "bc": {"top": "280 umol/L", "type": "dissolved"} # parsed automatically to 0.28 mmol/L
    "bc": {"top": 28.2, "type": "dissolved", "dbl": "0.5 mm"} # Optional Diffusive Boundary Layer
    ```
  - **Particulate Solids (`type: "particulate"`)**:
    Deposition at the sediment-water interface is specified as a depositional surface flux:
    ```text
    "bc": {"top": "12 umol/(cm^2 * yr)", "type": "particulate"}
    ```
    MANTIS automatically converts this to internal flux units ($\text{mol}/(\text{m}^2 \cdot \text{s})$) during initialization.

### 2. How Units Are Handled Internally by the Solver

Under the hood, MANTIS standardizes all transport and reaction PDEs into SI-derived volumetric units:
$$\text{Internal Concentration Unit} = \mathbf{mol/m^3} \quad (\equiv 1\ \text{mmol/L} \equiv 1\ \text{mM} \equiv 1000\ \mu\text{mol/L} \equiv 1\ \mu\text{mol/cm}^3)$$

Crucially, concentrations are defined **per unit volume of their respective phase**:
1. **Liquid Solutes**: Expressed per volume of porewater ($\text{mol / m}^3_{\text{liquid}} = \text{mmol / L}_{\text{porewater}}$). The bulk sediment inventory is $\phi \cdot C_{\text{liquid}}$.
2. **Solid Minerals**: Expressed per volume of pure solid sediment ($\text{mol / m}^3_{\text{solid}} = \text{mmol / L}_{\text{solid}}$). The bulk sediment inventory is $(1 - \phi) \cdot C_{\text{solid}}$.
3. **Phase Transfer & Reactions**: Whenever a reaction couples a solid phase and a liquid phase (for example, $\text{Fe}^{2+} + \text{HS}^- \to \text{FeS}$), MANTIS scales the source and sink terms with the appropriate phase ratio ($(1-\phi)/\phi$ or $\phi/(1-\phi)$) so that mass is strictly conserved.
4. **Reaction Rates**: Solved internally in $\text{mol} / (\text{m}^3_{\text{phase}} \cdot \text{s})$.
5. **Depths and Transport**: Depth $z$ is in meters ($\text{m}$), time $t$ in seconds ($\text{s}$), and diffusion coefficients in $\text{m}^2/\text{s}$.

### 3. Output Units in Result Files (`.csv`)

When the simulation completes or exports data (`<model_name>.csv`):
- `c_<species>`: Reported directly in **$\text{mmol/L}$** (numerically equivalent to $\text{mol/m}^3$):
  - Solutes (e.g., `c_SO4`, `c_O2`, `c_TS2`) are in **$\text{mmol / L of porewater}$**.
  - Solids (e.g., `c_FeS`, `c_Fe3`, `c_POC_slow`) are in **$\text{mmol / L of pure solid}$**.
- `f_<species>`: Net production/consumption reaction rates in **$\text{mol} / (\text{m}^3 \cdot \text{s})$**.
- `d_<species>` and `d34S_<species>`: Isotope delta values (e.g., `d_SO4`, `d34S_SO4`, `d34S_FeS2`) in **$\text{milli Urey}$** ($\text{mUr}$, relative to reference standards such as VCDT or VPDB).
- `z`: Depth below sediment-water interface in **meters** ($\text{m}$).
- `phi` and `w`: Porosity (fraction, dimensionless) and burial velocity ($\text{m/s}$).

### 4. Converting Solid Concentrations to Standard Geochemical Units

Geochemists typically report solid mineral abundances as dry-sediment mass fractions ($\text{wt}\%$) or mass-normalized molar concentrations ($\mu\text{mol/g}$). MANTIS includes built-in conversion utilities in `mantis.diff_lib`:

#### A. Converting to $\mu\text{mol / g}$ (micromoles per gram of dry sediment)
To convert a solid concentration $C_{\text{solid}}$ ($\text{mmol/L of solid} = \mu\text{mol/cm}^3_{\text{solid}}$) to $\mu\text{mol/g}$ of dry sediment, divide by the mineral grain density $\rho_{\text{sed}}$ (typically $2.6\ \text{g/cm}^3$ for siliciclastic sediments):
$$C_{\mu\text{mol/g}} = \frac{C_{\text{solid}}}{\rho_{\text{sed}}}$$

In Python:
```python
from mantis.diff_lib import solid_conc_to_umol_per_g

# df.c_FeS is in mmol/L of solid; density in g/cm^3
c_fes_umol_g = solid_conc_to_umol_per_g(df.c_FeS, d=2.6)
```

#### B. Converting to Dry Weight Percent ($\text{wt}\%$)
To convert $C_{\text{solid}}$ to weight percent of dry sediment matrix given molecular weight $M_w$ ($\text{g/mol}$):
$$\text{wt}\% = \frac{C_{\text{solid}} \cdot M_w}{\rho_{\text{sed}} \cdot 10^4}$$

In Python:
```python
from mantis.diff_lib import solid_conc_to_wt_percent

# Example: FeS (Mw = 87.91 g/mol, density = 2.6 g/cm^3)
fes_wt_pct = solid_conc_to_wt_percent(df.c_FeS, mw=87.91, d=2.6, phi=df.phi)
```


## Step 2: Defining Kinetic Constants (`kinetic_constants_*.py`)

The kinetic constants file defines how fast reactions occur and the affinity thresholds of microbial communities.

In diagenesis, we typically distinguish:

1.  **Rate Constants ($k$)**: The maximum rate of a reaction under substrate-saturated conditions (e.g., $k_{\text{POC}}$ in $\text{yr}^{-1}$ or $\text{s}^{-1}$, or mineral precipitation rates $k_{\text{prec}}$).
2.  **Half-Saturation Constants ($K_m$)**: The concentration of an electron acceptor at which the microbial degradation rate reaches $50\%$ of its maximum. Lower $K_m$ means the microbe is efficient at scavenging low concentrations.
3.  **Inhibition Constants**: The threshold concentration of a more energetically favorable oxidant that shuts down a less favorable pathway (e.g., $O_2$ inhibiting iron reduction or sulfate reduction).


### Structure of `kinetic_constants_slow.py`

Constants are encapsulated in a callable function `get_reaction_constants(pH, phi)` so that pH-dependent speciation and porosity-dependent conversions can be calculated automatically:

```python
def get_reaction_constants(pH, phi):
    """Returns reaction rate constants and half-saturation values for the diagenetic system."""

    # 1. Primary organic carbon degradation constants
    velde = {
        "POC_fast": 1.2e-9,        # Organic carbon degradation rate (s^-1)
        "POC_slow": 3.8e-11,       # Refractory carbon degradation rate (s^-1)
        "TS2_O2": 1.5e-6,          # Sulfide oxidation rate
        "Fe3_hs": 8.0e-7,          # Sulfide-mediated iron reduction rate
        "Fe2_O2": 2.5e-5,          # Ferrous iron oxidation rate
    }

    # 2. Mineral kinetics (Precipitation, Dissolution, Solubility Products)
    velde["FeS_isp"] = 2.0e-3      # FeS apparent precipitation rate
    velde["FeS_isd"] = 1.0e-5      # FeS apparent dissolution rate
    velde["FeS_sp"] = 0.3162       # FeS apparent solubility product

    # 3. Microorganism half-saturation & inhibition thresholds (mmol/L)
    # All half-saturation constants (K_*) belong here in the kinetic constants file:
    k_values = {
        "K_O2": 0.001,             # O2 limitation / inhibition threshold (1 uM)
        "K_O2_TS2": 0.001,         # O2 threshold for sulfide oxidation
        "K_SO4": 0.9,              # Sulfate limitation threshold (0.9 mM)
        "K_Fe3": 1e-3,             # Solid reactive iron threshold
        "K_Fe3_diss_red": 10.4,    # Dissimilatory iron reduction threshold
    }

    return velde, k_values
```



## Step 3: Defining Reactions & Limiters (`fes_reactions.py`)

The reactions file is where you express your geochemical hypotheses. It contains three main sections:

1.  `limiters`: Kinetic switches and fractionating terms.
2.  `reactions`: The list of stoichiometric equations and their regulatory controls.
3.  `reaction_parameters` (or `fractionation_factors`): Isotope fractionation factors ($\alpha$), stoichiometric ratios, and limiter transition constants.



### Understanding Limiters & Inhibitors

In sediment columns, microbes follow the classic redox ladder (Froelich et al., 1979). To model this mathematically:

-   When a microbe **needs** an electron acceptor (like $O_2$ or $SO_4^{2-}$), its rate is multiplied by a **Monod Limiter**: $$ L_{\text{Monod}}(C) = \frac{C}{C + K} $$
-   When a pathway is **suppressed** because a higher-energy oxidant is present, its rate is multiplied by an **Inhibition Term**: $$ I(C) = \frac{K}{C + K} = 1 - \frac{C}{C + K} $$

In MANTIS, you do **not** need to write out explicit fraction formulas for every reaction. You define clean base limiters in a dictionary:

```python
limiters = {
    # Implicit form: 1 / (C + K)
    # MANTIS automatically generates the limiter [C / (C + K)] and inhibitor [K / (C + K)]!
    "O2": "1.0 / (c.O2 + K_O2)",
    "O2_TS2": "1.0 / (c.O2 + mp.K_O2_TS2)",
    "SO4": "1.0 / (c.SO4 + K_SO4)",
    "Fe3_diss_red": "1.0 / (np.maximum(c.Fe3, 0.0) + K_Fe3_diss_red)",

    # Custom isotope explicit fractionation terms
    "SO4_alpha_explicit": "c.SO4 / (c.SO4 + 1e-3)",
    "TS2_alpha_explicit": "c.TS2 / (c.TS2 + 1e-3)",
}
```



### Declaring Chemical Reactions

Reactions are declared as a list of dictionaries with natural chemical formulas. Notice how clean and readable the syntax is:

```python
reactions = [
    # -------------------------------------------------------------------------
    # 1. Aerobic Respiration: POC + O2 -> CO2
    # Limited by O2 availability
    # -------------------------------------------------------------------------
    {
        "reaction_name": "aerobic_respiration_fast",
        "reaction": "POC_fast + O2 -> CO2",
        "k_value_name": "POC_fast",
        "limiters": ["O2"],
    },

    # -------------------------------------------------------------------------
    # 2. Dissimilatory Iron Reduction: POC + 4 Fe3 -> 4 Fe2
    # Limited by reactive iron, INHIBITED by oxygen
    # -------------------------------------------------------------------------
    {
        "reaction_name": "dissimilatory_iron_reduction_fast",
        "reaction": "POC_fast + 4 Fe3 -> 4 Fe2_total",
        "k_value_name": "POC_fast",
        "limiters": ["Fe3_diss_red"],
        "inhibitors": ["O2"],
    },

    # -------------------------------------------------------------------------
    # 3. Sulfate Reduction: 2 POC + SO4 -> TS2 (Total Sulfide)
    # Limited by SO4, INHIBITED by oxygen and reactive iron
    # Includes sulfur isotope fractionation
    # -------------------------------------------------------------------------
    {
        "reaction_name": "sulfate_reduction_fast",
        "reaction": "2 POC_fast + SO4 -> TS2",
        "k_value_name": "POC_fast",
        "limiters": ["SO4"],
        "inhibitors": ["O2", "Fe3_diss_red"],
        "isotope_species": {"SO4": "SO4_32", "TS2": "TS2_32"},
        "isotope_fractionation": {
            "SO4_32 to TS2_32": ["msr_alpha", "SO4_alpha_explicit"]
        },
    },

    # -------------------------------------------------------------------------
    # 4. Secondary Sulfide Oxidation: HS- + 2 O2 -> SO4
    # -------------------------------------------------------------------------
    {
        "reaction_name": "hs_oxidation",
        "reaction": "1 HS + 2 O2 -> 1 SO4",
        "k_value_name": "TS2_O2",
        "limiters": ["O2_TS2"],
        "isotope_species": {"SO4": "SO4_32", "TS2": "TS2_32", "HS": "HS_32"},
    },

    # -------------------------------------------------------------------------
    # 5. Mineral Precipitation & Dissolution (e.g. Iron Monosulfide FeS)
    # -------------------------------------------------------------------------
    {
        "reaction_name": "FeS_precipitation_dissolution",
        "reaction": "Fe2_liq + HS -> FeS",
        "k_value_name": {
            "precipitation": "FeS_isp",
            "dissolution": "FeS_isd",
            "equilibrium": "FeS_sp",
        },
        "isotope_species": {"TS2": "TS2_32", "FeS": "FeS_32", "HS": "HS_32"},
    },
]

# -------------------------------------------------------------------------
# Reaction Parameters & Isotope Fractionation Factors
# -------------------------------------------------------------------------
# Fractionation factors (alpha), stoichiometric ratios, and numerical limiter
# transition factors belong here in the reactions file:
reaction_parameters = {
    # Isotope fractionation factors
    "msr_alpha": 1.07,
    "TS2_O2_alpha": 0.995,
    "S0_O2_alpha": 1.0,
    "dispro_hs_alpha": 0.993,
    "dispro_SO4_hs_split": 3.0,
    "h2s_hs_alpha": 0.99991542,
    # Primary respiration stoichiometry
    "POC_O2_ratio": 1.27,
    # Transition factors for isotope fractionation limiters
    "K_epsilon_msr": 0.2,
    "K_epsilon_TS2_O2": 0.01,
}
fractionation_factors = reaction_parameters
```

### How MANTIS Classifies Reactions (Species-Agnostic Engine)

Historically, reaction-transport codes often relied on hardcoded species names (e.g., checking explicitly for `"POC"` and `"O2"`). In MANTIS, the equation compiler is **fully species-agnostic**. You can invent and model arbitrary reaction pathways (such as denitrification, manganese reduction, carbonate precipitation, or methane oxidation) without modifying any internal solver code.

The compiler determines how to assemble the PDE transport coupling and Jacobian using two primary properties:

1. **Phase Permutations** (derived directly from `species.py`):
   - **`liquid_to_liquid`**: All reactants and products are dissolved porewater solutes (e.g., $\text{HS}^- + 2\,\text{O}_2 \to \text{SO}_4^{2-}$). Reaction rates are scaled by porosity $\phi$.
   - **`solid_to_liquid`**: One or more reactants are solid minerals or particulate organic carbon, and all products are dissolved (e.g., $\text{POC} + \text{O}_2 \to \text{CO}_2$). Reaction rates are scaled by the solid sediment fraction $(1 - \phi)$.
   - **`liquid_to_solid`**: Dissolved porewater reactants precipitate into solid minerals (e.g., $\text{Fe}^{2+} + \text{HS}^- \to \text{FeS}$).
   - **`solid_transformation`**: Solid reactants react to form solid products (e.g., $\text{FeS} + \text{S}^0 \to \text{FeS}_2$ or dissimilatory iron reduction $\text{POC} + 4\,\text{FeOOH} \to 4\,\text{Fe}^{2+}$).

2. **Thermodynamic / Kinetic Mechanisms**:
   - **`general_kinetic`**: Standard single- or multi-reactant irreversible kinetics with arbitrary Monod limiters and inhibition factors. Automatically coupled via `add_coupled_reaction()`.
   - **`reversible_saturation`**: Dynamic mineral precipitation/dissolution driven by the saturation state $\Omega = Q / K_{sp}$. Uses a $C^1$ cubic Hermite spline transition around equilibrium ($\Omega = 1$) to eliminate discontinuities.
   - **`disproportionation`**: A single reactant splits into multiple redox products (e.g., $4\,\text{S}^0 + 4\,\text{H}_2\text{O} \to 3\,\text{H}_2\text{S} + \text{SO}_4^{2-} + 2\,\text{H}^+$).

#### Automated SymPy Chemical Jacobian
Whenever reaction equations are compiled, MANTIS uses `SymPy` to algebraically differentiate the rate laws with respect to every participating state variable:
$$\mathbf{J}_{s, m} = \frac{\partial R_s}{\partial C_m}$$
Because differentiation is performed symbolically:
- Exact machine-precision Jacobians are generated for any arbitrary reaction network.
- Speciation chain rules (such as $\frac{\partial R}{\partial [\text{TS}_2]} = \frac{\partial R}{\partial [\text{HS}^-]} \cdot \alpha_{\text{HS}^-}$) are applied automatically.
- The monolithic Newton solver achieves quadratic convergence without human derivative derivation errors.



## Step 4: The Master Run Script (`run_velde_model_fast.py`)

The driver script is your experiment control panel. It specifies the sediment accumulation rate, bioturbation mixing, domain length, and simulation duration.

```python
#!/usr/bin/env python
from mantis import DiageneticModel
from mantis.diff_lib import get_total_delta

if __name__ == "__main__":
    # 1. Initialize the model
    # MANTIS automatically associates species.py in the same folder
    model = DiageneticModel(
        reactions="fes_reactions.py",
        constants="kinetic_constants_slow.py",
        layout="plot_layout.py",
    )

    # 2. Configure sediment environmental parameters using human-friendly units!
    model.configure(
        max_depth="0.5 m",             # Total sediment column depth (50 cm)
        display_length="0.5 m",        # Depth shown in plots
        t_end="2 kyr",                 # Integrate for 2,000 years
        dt_max="1 year",               # Maximum allowed timestep
        dt_init="1 day",               # Initial startup timestep
        dt_min="1 second",             # Minimum timestep before solver gives up
        w="0.2 cm/yr",                 # Sediment accumulation / burial rate
        phi=0.8,                       # Porosity (80% porewater, 20% solids)
        isotopes=True,                 # Enable 34S/32S isotope tracking

        # Macrofaunal bioturbation mixing (sediment surface down to 6 cm)
        bioturbation={
            "BT0": "4 cm^2/year",      # Surface bio-diffusion coefficient
            "profile": "sigmoidal",    # "sigmoidal" (default) or "exponential"
            "depth": "6 cm",           # Mixed layer depth (for sigmoidal)
            "attenuation": "2 cm",     # Transition attenuation width
        },

        # Porewater bio-irrigation (burrow flushing by worms/crustaceans)
        bioirrigation={
            "BI0": 0.0,                # Set to 0.0 to disable; or e.g. "20 /year"
            "depth": "0 cm",
        },

        # Spatial grid & reaction zone refinement (optional, defaults provided)
        initial_spacing="0.1 mm",      # Spacing at sediment-water interface (z=0)
        reaction_zone=("0.1 cm", "1.0 cm"), # Critical depth window of sharp redox fronts
        reaction_zone_spacing="0.05 mm",    # Dense, uniform spacing in the reaction zone
        max_spacing="1.0 cm",          # Coarser spacing in deep sediment
    )

    # 3. Allow CLI overrides (e.g. running --video from the terminal)
    model.parse_cli_args()

    # 4. Run the simulation
    result = model.run()

    # 5. Save outputs and print summary
    df, csv_path = result.save_csv()
    result.save_state()

    if result.mp.isotopes:
        print(f"Buried Delta34S = {get_total_delta(result.c, result.mp):.2f} mUr")

    # Render final publication plot
    result.plot_summary()
```

### Bioturbation Mixing Profiles: Sigmoidal vs. Exponential

Bioturbation represents the mechanical mixing of sediment particles (solids) by burrowing and deposit-feeding benthic macrofauna. MANTIS provides two standard formulations selectable via the `profile` parameter in the `bioturbation` dictionary:

#### 1. Sigmoidal Profile (`profile="sigmoidal"`, Default)
Based on the formulation of van de Velde and Meysman (2016), this profile represents a distinctly mixed surface layer of thickness $z_{\text{mix}}$ (`depth`) underlain by a smooth transition zone of width $x_{\text{atten}}$ (`attenuation`):

$$D_b(z) = \frac{D_{b0}}{1 + \exp\left(\frac{z - z_{\text{mix}}}{0.25 \, x_{\text{atten}}}\right)}$$

- **When to use**: Ideal for shallow coastal and shelf environments where continuous burrowing activity maintains a uniform, constant mixing intensity in the top few centimeters before tapering off sharply below the mixed layer.
- **Configuration**:
  ```python
  bioturbation={
      "BT0": "4 cm^2/year",      # Surface bio-diffusivity Db0
      "profile": "sigmoidal",    # Default
      "depth": "6 cm",           # Mixed layer thickness z_mix
      "attenuation": "2 cm",     # Transition zone width x_atten
  }
  ```

#### 2. Exponential Profile (`profile="exponential"`)
A classical diagenetic formulation where bio-diffusive mixing intensity is highest at the sediment-water interface ($z = 0$) and continuously decays with depth without a flat shoulder:

$$D_b(z) = D_{b0} \exp\left(-\frac{z}{x_{\text{atten}}}\right)$$

- **When to use**: Widely used in benchmark reactive-transport models (such as SedTrace / SimpleFe; Du 2023) and deep-sea or low-activity sediments where faunal abundance and mixing frequency attenuate exponentially from the sediment surface.
- **Configuration**:
  ```python
  bioturbation={
      "BT0": "10 cm^2/year",     # Surface bio-diffusivity Db0
      "profile": "exponential",
      "attenuation": "3 cm",     # e-folding attenuation depth scale (also accepts "decay")
  }
  ```

Both formulations are automatically converted into SI units ($\text{m}^2/\text{s}$ and $\text{m}$) by Pint, and are coupled into the particulate transport stencils.


### Diffusive Boundary Layer (DBL) Setup

In benthic environments, a thin, hydrodynamic boundary layer of stagnant or laminar water sits directly above the sediment-water interface ($z = 0$). This **Diffusive Boundary Layer (DBL)** has a characteristic thickness $\delta$ (typically $0.2 - 1.0\ \text{mm}$; in many benchmarks such as Du 2023, $\delta = 0.5\ \text{mm} = 0.05\ \text{cm}$).

Within this viscous sublayer, turbulence from overlying bottom water currents is damped, meaning solute transport occurs solely via molecular diffusion through free water ($\phi = 1$, without tortuosity hindrance).

#### 1. Why Model the DBL?
- **Retention vs. Ventilation**: Without a DBL (i.e. pure Dirichlet $C(0) = C_{\mathrm{BW}}$), dissolved metabolic products (such as $\mathrm{Fe}^{2+}$ or $\mathrm{H_2S}$) are instantaneously flushed away into the well-mixed water column. With a DBL, solutes can build up slightly at the sediment-water interface, increasing interfacial concentration gradients and driving authigenic mineral precipitation (e.g., $\mathrm{FeS}$ or carbonate minerals).
- **Benchmark Consistency**: Benchmark models like SedTrace (Du 2023) incorporate an explicit DBL Robin boundary condition ($-\phi_0 D_{\mathrm{sed}} \left.\frac{\partial C}{\partial z}\right|_0 = \frac{D_{\mathrm{free}}}{\delta} (C_0 - C_{\mathrm{BW}})$).

#### 2. Physical & Mathematical Formulation in MANTIS
MANTIS represents the DBL as an exact finite-volume series resistance between the center of the first sediment grid cell ($z_0$) and the overlying well-mixed bottom water ($z = -\delta$):

$$\frac{1}{D_{\mathrm{top}}} = R_{\mathrm{sed}} + R_{\mathrm{DBL}} = \frac{z_0}{\phi_0 D_{\mathrm{sed}}} + \frac{\delta}{D_{\mathrm{free}}}$$

$$\implies D_{\mathrm{top}} = \frac{1}{\frac{z_0}{\phi_0 D_{\mathrm{sed}}} + \frac{\delta}{D_{\mathrm{free}}}}$$

- When $\delta = 0$, $R_{\mathrm{DBL}} = 0$, collapsing smoothly to the standard Dirichlet condition $D_{\mathrm{top}} = \frac{\phi_0 D_{\mathrm{sed}}}{z_0}$.
- Free molecular diffusion $D_{\mathrm{free}}$ is related to the tortuous sediment porewater diffusion coefficient by Boudreau's relation: $D_{\mathrm{free}} = D_{\mathrm{sed}} \cdot (1 - \ln \phi_0^2)$.

#### 3. How to Configure the DBL in MANTIS
You can activate and set the DBL thickness globally across all dissolved solutes in your driver script via `model.configure()`:

```python
model.configure(
    ...,
    diffusive_boundary_layer="0.5 mm",  # Or "0.05 cm"; accepts any Pint length unit
)
```
*(The alias `dbl="0.5 mm"` is also accepted; `"delta"` is strictly reserved for isotopic composition in species definitions).*

Alternatively, you can specify DBL parameters for individual species in `species.py`:
```python
"SO4": {
    "phase": "liquid",
    "bc": {
        "top": 28.0,
        "type": "dissolved",
        "dbl": "0.5 mm",
    },
}
```


### Spatial Grid Configuration & Reaction Zone Resolution

Diagenetic models must resolve spatial scales ranging from sub-millimeter oxygen penetration depths to deep sediment columns spanning decimeters to hundreds of meters. A uniform fine grid across a 50-meter column would require millions of cells and slow down calculations unnecessarily.

MANTIS solves this using an adaptive geometric-mesh builder (`make_grid2`) that refines where steep gradients occur and expands into the deep sediment.

### 1. How Grid Generation Works (`make_grid2`)
The grid structure is divided into four seamless stages from the sediment surface to depth:
1. **Surface Section**: Starts at `initial_spacing` and geometrically refines cells downward toward `reaction_zone[0]`.
2. **Reaction Zone Section**: Maintains a constant, ultra-fine cell size (`reaction_zone_spacing`) across the depth interval `reaction_zone = (z_start, z_end)`.
3. **Coarsening Section**: Below the reaction zone, cell size expands geometrically by factor $r = 1.05$ per cell.
4. **Deep Uniform Section**: Once cell size reaches `max_spacing`, the remainder of the domain down to `max_depth` uses uniform large cells.

#### 2. The Geochemical Importance of Reaction Zone Resolution
Why does grid spacing in the reaction zone matter so much?
Fast bimolecular redox reactions—such as ferrous iron oxidation ($\text{Fe}^{2+} + \text{O}_2 \to \text{Fe}^{3+}$) or sulfide oxidation—occur across exceedingly narrow interfaces. The physical reaction front width $\delta_{\text{phys}}$ is governed by the balance between molecular diffusion and reaction kinetics:
$$\delta_{\text{phys}} \sim \sqrt{\frac{D_{\text{eff}}}{k \cdot C_{\text{scale}}}}$$
For transport-controlled reactions with large kinetic constants ($k = 10^7\ \text{M}^{-1}\text{s}^{-1}$), $\delta_{\text{phys}}$ is often **under $0.05\text{ mm}$**.

If your numerical grid cell size $\Delta z$ is larger than $2 \delta_{\text{phys}}$, the local Damköhler number is large ($Da_{II} \gg 1$). The finite-volume discretization cannot resolve the sharp counter-gradients, producing a numerical artifact called **artificial reactant coexistence**:
- $\text{Fe}^{2+}$ and $\text{O}_2$ artificially overlap in the same cell.
- Reactants diffuse past each other instead of reacting completely at the interface.
- Apparent reaction rates become smeared over several millimeters to centimeters.

#### 3. Automated Coexistence Diagnostic
MANTIS automatically runs a post-simulation front diagnostic (`check_reactant_coexistence`) at the end of every run. If any reaction front is under-resolved, MANTIS outputs a clear warning:
```text
[WARNING: Numerical Reactant Co-existence Detected]
  Reaction: 'Fe2 Oxidation by O2' (fe2 + O2)
  Front location (peak): z = 0.21 cm (depth cell index: 21)
  Local grid resolution: dz = 0.10 mm
  Estimated physical front width: 0.01 mm
  Front cell span (N_cells): 6 cell(s)
  Reasons: Local grid cell dz (0.10 mm) > 2x physical front width (0.01 mm)
  --> Action recommended: In make_grid2(), refine the mesh around z = [0.19, 0.23] cm to reaction_zone_spacing <= 0.01 mm.
```

When you see this warning, simply update your driver script's `model.configure()` with the recommended reaction zone window and spacing:
```python
model.configure(
    reaction_zone=("0.15 cm", "0.30 cm"),
    reaction_zone_spacing="0.01 mm",
)
```


## Step 5: Designing Plot Layouts (`plot_layout.py`)

Geochemists usually want to see:

-   Dissolved porewater concentrations in $\text{mmol/L}$ (or $\mu\text{mol/L}$).
-   Solid minerals and organic carbon in **dry weight percent** ($\text{wt}\%$) or $\mu\text{mol/g}$.
-   Reaction rates and physical transport (e.g., bioturbation coefficient $D_{\text{bio}}$).
-   Stable isotopes in per mil or milli-urey ($\delta^{34}\text{S}$, $\delta^{18}\text{O}$, $\delta^{13}\text{C}$).

MANTIS uses a clean, purely declarative layout schema in `plot_layout.py`. You do not need to write procedural parsing functions or manually convert units—MANTIS automatically converts solid concentrations to $\text{wt}\%$, links color palettes across related elements, and formats axes.


### Declarative Schema Overview

A typical `plot_layout.py` defines specifications for liquid species, solid species, reaction rates, and isotopes, and assigns them flexibly to panels via `plot_panels`:

```python
"""plot_layout.py - Declarative plot configuration."""

# ---- 1. Color Map and Style Setup ----
# Standard colormap name ('default', 'tab10', 'Set1'), a dict of species colors,
# or per-species overrides.
colormap = "default"

# Override styles for specific species (colors, line widths, background fills, etc.)
colormap_override = {
    "D_bio": {"color": "#bbbbbb", "lw": 2, "fill_only": True, "zorder": -10},
}

# ---- 2. Species & Panel Specifications ----
# Liquid dissolved porewater species
liquid_species_1 = liquid_species({
    "scale": "log",                              # 'log' or 'lin' / 'linear'
    "species": ["SO4", "H2S", "Fe2_dis"],
    "lim": (1e-4, 30),                           # Axis limits
    "units": "mmol/l",                           # Axis label
    "position": "left",                          # Default: "left"
})

# Solid minerals and particulates (supports automatic wt% dry sediment conversion)
solid_species_1 = solid_species({
    "scale": "lin",
    "units": "wt%",                              # Automatic wt% conversion
    "relative_to": "Fe",                         # Molecular weight reference (e.g. Fe, S, C; None for no scaling)
    "species": ["FeS"],
    "lim": (0, 2),
    "position": "right",                         # Attach as secondary right axis
})

# Reaction rates panel
# Defaults to symlog scale to handle both positive (production) and negative (consumption) rates
reaction_rates_1 = reaction_rates({
    "scale": "symlog",
    "species": ["POC", "SO4", "H2S", "FeS", "D_bio"],
    "position": "left",
})

# Isotope data panel (Optional)
# Multiple isotope systems (e.g. d34S and d18O) automatically create twin Y-axes.
# 4-element tuples specify broken axes: (min_lower, max_lower, min_upper, max_upper)
isotopes_1 = isotopes({
    "species": ["d34S_SO4", "d34S_FeS"],
    "lim": (-45, -37, 20, 30),                   # Broken y-axis excluding [-37, 20]
    "position": "left",
})

# ---- 3. Flexible Panel Assignment & Ordering ----
# Each sub-list in plot_panels defines a subplot from top to bottom.
# Panels can be combined, reordered, or repeated arbitrarily.
plot_panels = [
    [liquid_species_1, solid_species_1],  # Panel 1: Liquids (left) + Solids (right)
    [reaction_rates_1],                   # Panel 2: Reaction rates
    [isotopes_1],                         # Panel 3: Isotopes
]

# ---- 4. Optional Comparison Against Measured Data ----
# comparison_data = "measured_core_data.csv"
# species_mapping = {
#     "z": ["depth_cm", 1 / 100],                # Map depth column with unit conversion
#     "SO4": ["C_SO4", 1.0],                     # Map porewater sulfate
#     "FeS": ["C_FeS", "wt"],                    # Map solid FeS with automatic wt% conversion
# }
```


### Key Features of the Declarative Layout

1. **Flexible Panel Assignment (`plot_panels`)**:
   Panels can be arranged in any order and composed of any combination of specs. For instance, you can place isotopes on top, create multiple dedicated liquid panels, or place solids on their own separate subplot. The `"position": "left" | "right"` setting directs curves to primary left or secondary twin right axes.

2. **Unified Color Palette Across Panels**:
   Plot elements that relate to the same chemical entity share the exact same color throughout the entire figure. For example, `SO4` in porewater, `SO4` in reaction rates, and `d34S_SO4` in isotopes all automatically share the same color. Colors are assigned generically by cycling through `tab10` (or your chosen `colormap`).

3. **Automatic Solid Weight Percent ($\text{wt}\%$) Conversion**:
   When `"units": ["wt"]` (or `"wt%"`) and `"relative_to": "Fe"` are defined, MANTIS uses sediment porosity ($\phi$) and dry grain density ($\rho_s = 2.6\text{ g/cm}^3$) to convert volumetric concentrations ($\text{mmol/L}$) to dry weight percent:
   $$\text{wt}\% = \frac{C_{\text{solid}} \cdot M_{\text{rel}} \cdot (1 - \phi)}{\rho_s \cdot 10^4}$$
   When `relative_to` is omitted or set to `None`, units like `"pH"` or `"dimensionless"` are formatted as-is without appending any relative element.

4. **Multi-Axis and Broken Axis Support**:
   - Secondary solid mineral groups or secondary isotope systems (e.g. $\delta^{34}\text{S}$ on left and $\delta^{18}\text{O}$ on right) automatically create twin axes.
   - Broken axes can be declared cleanly via a 4-tuple `lim=(min1, max1, min2, max2)` to zoom into separate upper and lower isotope ranges.

5. **Background Fills (`fill_only`)**:
   Physical parameters like bioturbation mixing ($D_{\text{bio}}$) or porosity can be plotted as shaded background areas behind lines by setting `"fill_only": True` (or `"fillonly": True`) and a negative `"zorder"`.



## Step 6: Running the Model & Producing Time-Lapse Videos

### Standard Execution

To run your simulation, simply run the driver script in your terminal:

```sh
python run_velde_model_fast.py
```

During execution, MANTIS outputs progressive milestones indicating physical time reached, current timestep $\Delta t$, and Newton iterations:

```
Grid generated with 1108 points.
Starting Monolithic Block Newton-Raphson Solver. dt_init: 24.00 h
Monolithic solver: Analytical chemical Jacobian enabled.
Step 25, time=8.6 week (0%), dt=36.68 h, iters=6
Step 50, time=20.2 week (0%), dt=5.52 d, iters=5
Step 100, time=26.7 month (0%), dt=25.23 d, iters=6
...
Step 2129, time=2.0 ka (100%), dt=52.18 week, iters=1
Final Report: Maximum steps or end time reached in 2129 steps. Total Wall Time: 26.2s
State saved to run_velde_model_state.npz
Saved final PDF plot: run_velde_model.pdf
d34S = -23.59 mUr
```


### Generating a Time-Lapse Video

To observe the dynamics as reaction fronts establish, add the `--video` CLI flag followed by the recording interval:

```sh
python run_velde_model_fast.py --video "1 year"
```

This spawns an asynchronous background video worker that renders each snapshot into an MP4 file (`run_velde_model.mp4`) without slowing down the core simulation solver!


## Command-Line Options Reference

You can query all options and their configured defaults at any time:

```sh
python run_velde_model_fast.py --help
```

### Interactive Category Queries
You can also inspect specific sub-systems and their current runtime values by passing a category flag with no argument:

```sh
python run_velde_model_fast.py --solver   # Displays all numerical solver parameters
python run_velde_model_fast.py --video    # Displays all video recording options
python run_velde_model_fast.py --run      # Displays time, grid, and domain settings
```

### Complete Options Table

| Category | Flag | Argument | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Solver** | `--solver` | `monolithic` \| `coupled` | `monolithic` | Backend solver architecture |
| | `--solver-tol`, `--newton-tol` | `<float>` | `1e-3` | Newton-Raphson convergence residual tolerance |
| | `--numerical-jacobian` | (flag) | `False` | Use finite differences instead of SymPy analytical Jacobian |
| | `--max-newton-iters` | `<int>` | `12` | Maximum Newton iterations allowed per timestep |
| | `--newton-norm` | `linf` \| `wrms` | `linf` | Residual norm formulation ($L_\infty$ or WRMS) |
| | `--newton-target-optimal`| `<int>` | `4` | Target Newton iterations for timestep expansion |
| | `--newton-max-acceptable`| `<int>` | `8` | Maximum iterations before timestep reduction |
| | `--dt-tolerance` | `<float>` | `1e-5` | Adaptive timestep convergence tolerance |
| | `--force-regen`, `--rebuild-equations` | (flag) | `False` | Force regeneration of equation system & analytical Jacobian |
| **Video** | `--video` | `<interval>` | `None` | Live video recording frame interval (e.g. `'10 year'`, `'none'`) |
| | `--video-dpi` | `<int>` | `120` | Resolution (DPI) of recorded video frames |
| | `--video-queue-size` | `<int>` | `100` | Maximum video frame render queue size |
| | `--video-resume-threshold` | `<int>` | `50` | Video frame resume queue threshold |
| | `--video-workers` | `<int>` | auto | Dedicated background worker count for video rendering |
| **Run & Domain** | `--t-end` | `<duration>` | `1 kyr` | Total duration of simulation (e.g. `'10 kyr'`) |
| | `--dt-max` | `<timestep>` | `1 year` | Maximum adaptive timestep cap (e.g. `'1 year'`) |
| | `--dt-init` | `<timestep>` | `1 day` | Initial startup timestep size (e.g. `'1 day'`) |
| | `--dt-min` | `<timestep>` | `1 minute` | Minimum allowable timestep before aborting |
| | `--max-depth` | `<length>` | `2.0 m` | Total depth of sediment model column (e.g. `'10 m'`) |
| | `--display-length` | `<length>` | `2.0 m` | Plotted column depth for video and PDF summaries |


## Automatic Equation & Jacobian Regeneration

MANTIS automatically keeps your compiled equation system and analytical SymPy Jacobian synchronized with your geochemical input files.

### Do We Need to Regenerate Equations or the Jacobian When Files Change?
**Yes, absolutely.**
- **`species.py`**: Controls the active chemical species vector $\mathbf{c}$, system size $N$, phase partitioning (dissolved vs. solid vs. mixed), and companion isotope pairs ($^{32}\text{S}$ / $^{34}\text{S}$). Adding, removing, or changing species phases alters the dimensions ($N \times N$) of the block tridiagonal Jacobian matrix.
- **`limiters.py` (or `limiters` in `fes_reactions.py`)**: The analytical chemical Jacobian directly embeds the partial derivatives ($\frac{\partial L}{\partial C}$) of all Monod and inhibition terms. Any change in limiter mathematical forms or threshold formulas directly alters the Jacobian entries.
- **`kinetic_constants_*.py`**: Contains rate constants and half-saturation parameters ($K_m$) used in rate equations and validation checks.
- **`fes_reactions.py`**: Defines reaction stoichiometries, active pathways, and isotope fractionation mechanisms.
- **`generate_equation_system.py`**: The code generator itself.

### How MANTIS Handles It
Whenever `model.setup()` or `model.run()` is called:
1. MANTIS checks the file modification timestamps (`st_mtime`) of `fes_reactions.py`, `species.py`, `kinetic_constants_*.py`, `limiters.py`, and `generate_equation_system.py`.
2. If any of these files is **newer** than `<experiment>_equations.py`, MANTIS automatically re-parses the chemical network, derives all symbolic derivatives via SymPy, writes out a freshly updated equations module, invalidates import caches, and loads the new system without requiring manual intervention.
3. You can also force an immediate rebuild at any time using:
   ```sh
   python run_velde_model_fast.py --force-regen
   ```
   or in Python:
   ```python
   model.configure(force_regen=True)
   ```


## Solver Mechanics & Physical Realism (Geochemist's Alert!)

A critical strength of MANTIS's monolithic implicit Block Newton solver is that it is **unconditionally numerically stable**: unlike classical explicit methods (which blow up if the timestep $\Delta t$ exceeds seconds or hours), MANTIS will **not crash** even when taking timesteps of $\Delta t = 1\text{ year}$ or more.

However, **numerical stability is not the same as physical realism**. Here is why this distinction matters:

> **CAUTION: Large $\Delta t$ and Reaction Front Smearing** When a sediment column contains extremely fast, irreversible reactions between opposing fronts—such as the suboxic-anoxic transition zone where dissolved hydrogen sulfide ($\text{HS}^-$) meets oxygen ($\text{O}_2$) or ferric iron ($\text{Fe}^{3+}$)—the natural reaction zone may physically be only a few millimeters thick.
> 
> If the solver takes very large time steps ($\Delta t = 1\text{ year}$), diffusion carries $\text{HS}^-$ and $\text{O}_2$ across several grid cells within a single time step. While the implicit solver computes a mathematically stable solution, it will **artificially broaden (smear)** the reaction front. In the simulation, sulfide and oxygen might appear to coexist over several centimeters, whereas in nature they react immediately to near-zero levels within millimeters.

### How to Ensure Physical Realism:

1.  **Cap $\Delta t_{\text{max}}$**: If you are studying transient front migration or sharp redox boundaries, set `dt_max`"1 month"= or `dt_max`"6 month"= rather than multiple years.
2.  **Refine the Reaction Zone Grid**: Use MANTIS's `make_grid2` with high resolution (e.g. $\Delta z \le 0.5\text{ mm}$) across the critical depth where the redox front sits.
3.  **Check Front Coexistence**: Verify that reactants do not unnaturally overlap by running diagnostic checks across your front.



# Summary & Next Steps

You now have a complete overview of MANTIS's modular modeling system:

1.  **Species** are declared in `species.py` with phases, diffusion laws, and boundary conditions.
2.  **Kinetics** are set in `kinetic_constants_*.py` with half-saturation and rate values.
3.  **Reactions** are written in `fes_reactions.py` with natural chemistry and Monod/inhibition controls.
4.  **Environmental Settings** (depth, burial, time, bioturbation) are configured in your driver script.
5.  **Visuals** are customized in `plot_layout.py` in geochemist-friendly units ($\text{wt}\%$, $\text{mUr}$).

Happy modeling!
