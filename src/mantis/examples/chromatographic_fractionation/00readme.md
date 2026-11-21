# Chromatographic Isotope Fractionation Model

## Overview

This benchmark model demonstrates **chromatographic isotope separation** during 1D advective-dispersive reactive transport in Mantis. It simulates the invasion of dissolved ferrous iron ($\text{Fe}^{2+}$) into an initially iron-free sediment column, accompanied by fast linear equilibrium sorption onto particulate matter:

$$\text{Fe}^{2+}_{\text{liq}} \longleftrightarrow \text{Fe}^{2+}_{\text{sorbed}}$$

Because the sorption equilibrium exhibits an isotope fractionation effect ($\alpha_{\text{eq}} \neq 1$), the companion isotopes ($^{54}\text{Fe}$ and $^{56}\text{Fe}$) experience different partitioning between the mobile porewater and immobile solid matrix. This differential retardation causes the lighter isotope ($^{54}\text{Fe}$) to migrate ahead of the heavier isotope ($^{56}\text{Fe}$), generating a characteristic negative $\delta^{56}\text{Fe}$ excursion at the leading invasion front.

---

## Theory & Mathematics

### 1. Lumped Mixed-Phase Transport

In Mantis, lumped pools (`phase="mixed"`) partition dynamically into dissolved and solid fractions:

$$C_{\text{liq}} = f_{\text{dis}} \cdot C_{\text{total}}, \qquad C_{\text{sorbed}} = (1 - f_{\text{dis}}) \cdot C_{\text{total}}$$

The 1D conservation equation for total iron in porosity $\phi$ with porewater advection $u$ and molecular diffusion $D_{\text{mol}}$ is:

$$\frac{\partial C_{\text{total}}}{\partial t} = \frac{\partial}{\partial z}\left( D_{\text{eff}} \frac{\partial C_{\text{total}}}{\partial z} \right) - \frac{\partial}{\partial z}\left( w_{\text{eff}} C_{\text{total}} \right)$$

where the effective transport parameters are retarded by sorption:
- Effective diffusion: $D_{\text{eff}} = f_{\text{dis}} \cdot D_{\text{mol}}$
- Effective advection: $w_{\text{eff}} = f_{\text{dis}} \cdot u$
- Retardation factor: $R = \frac{1}{f_{\text{dis}}}$

### 2. Isotopic Companion Partitioning & Retardation

For an isotope companion ($^{54}\text{Fe}$), equilibrium fractionation is governed by:

$$\alpha_{\text{eq}} = \frac{(^{54}\text{Fe} / ^{56}\text{Fe})_{\text{sorbed}}}{(^{54}\text{Fe} / ^{56}\text{Fe})_{\text{liq}}}$$

When the heavy isotope $^{56}\text{Fe}$ sorbs preferentially, $\alpha_{\text{eq}} < 1$ for the light isotope $^{54}\text{Fe}$ (e.g. $\alpha_{\text{eq}} = 0.998$). 

Under linear equilibrium, the dissolved fraction of the $^{54}\text{Fe}$ companion ($f_{\text{dis, 54}}$) is:

$$f_{\text{dis, 54}} = \frac{f_{\text{dis}}}{f_{\text{dis}} + \alpha_{\text{eq}} \cdot (1 - f_{\text{dis}})}$$

Since $\alpha_{\text{eq}} < 1$:
$$f_{\text{dis, 54}} > f_{\text{dis}} \implies R_{54} < R_{\text{bulk}}$$

Because the $^{54}\text{Fe}$ companion has a higher dissolved fraction, its effective advection velocity $w_{\text{eff, 54}} = f_{\text{dis, 54}} \cdot u$ and effective diffusivity $D_{\text{eff, 54}} = f_{\text{dis, 54}} \cdot D_{\text{mol}}$ are greater than those of bulk iron ($^{56}\text{Fe}$).

### 3. Chromatographic Separation

As the chemical front invades:
1. $^{54}\text{Fe}$ outpaces $^{56}\text{Fe}$ at the leading edge.
2. The ratio $(^{56}\text{Fe} / ^{54}\text{Fe})$ at the leading edge decreases below the boundary input ratio.
3. This creates a prominent negative $\delta^{56}\text{Fe}$ anomaly at the plume front.

---

## Model Specifications

| Parameter | Symbol | Value | Unit |
| :--- | :--- | :--- | :--- |
| Column length | $L$ | 20 | $\text{cm}$ |
| Grid resolution | $\Delta z$ | 1 | $\text{mm}$ (uniform 200 cells) |
| Porosity | $\phi$ | 0.8 | dimensionless |
| Advection velocity | $u$ | 5.0 | $\text{cm/yr}$ (downward) |
| Sediment burial | $w$ | 0.0 | $\text{cm/yr}$ |
| Temperature | $T$ | 10.0 | $^\circ\text{C}$ |
| Simulation time | $t_{\text{end}}$ | 0.5 | $\text{yr}$ |
| Time stepping | $\Delta t_{\text{init}} / \Delta t_{\text{max}}$ | 0.5 / 2.0 | $\text{day}$ |
| Total Fe boundary | $C_{\text{top}}$ | 1.0 | $\text{mM}$ |
| $\delta^{56}\text{Fe}$ boundary | $\delta_{\text{top}}$ | 0.0 | $\text{mUr}$ |
| Bulk dissolved fraction | $f_{\text{dis}}$ | 0.5 | dimensionless |
| Fractionation factor | $\alpha_{\text{eq}}$ | 0.998 | dimensionless |

---

## Analytical Benchmark Results

After $t = 0.5\text{ yr}$ (95 time steps, ~0.7 s wall time), the simulated profiles closely track the exact semi-infinite analytical Ogata-Banks retarded advection-dispersion solution (`analytical_ogata_banks.csv`):
- Mean grid-weighted normalized error against analytical Ogata-Banks is **0.16%** for concentration profiles (`Fe2_total`, `Fe2_liq`, `Fe2_sorbed`) and **0.36%–0.59%** for $\delta^{56}\text{Fe}$ profiles.
- The front midpoint ($C_{\text{total}} = 0.5\text{ mM}$) advances to $z = 5.55\text{ cm}$ (half the unretarded advection-dispersion front velocity plus effective retarded diffusion).
- Near the upper boundary ($z = 0.05\text{ cm}$), bulk $\delta^{56}\text{Fe} \approx -0.003\text{ mUr}$, while equilibrium partitioning maintains a constant $\Delta^{56}\text{Fe}_{\text{sorbed-liq}} \approx +2.13\text{ mUr}$ offset ($\delta^{56}\text{Fe}_{\text{liq}} = -1.07\text{ mUr}$, $\delta^{56}\text{Fe}_{\text{sorbed}} = +1.06\text{ mUr}$).
- Toward the leading invasion front ($z = 19.95\text{ cm}$), differential retardation ($R_{54} < R_{56}$) drives $\delta^{56}\text{Fe}_{\text{total}}$ down to $-4.28\text{ mUr}$ ($\delta^{56}\text{Fe}_{\text{liq}} = -5.34\text{ mUr}$).

---

## Declarative Configuration

In `species.py`, the partitioning and fractionation are defined declaratively:

```python
speciation = {
    "fe2_partitioning": {
        "type": "partition",
        "primary": "Fe2_total",
        "fractions": {
            "dissolved": "0.5",
            "sorbed": "0.5",
        },
        "equilibrium_fractionation": {
            "standard": "d56Fe",
            "alpha": {"Fe2_sorbed_54_to_Fe2_liq_54": 0.998},
        },
        "species": {
            "Fe2_liq": {
                "phase": "liquid",
                "alias_of": "Fe2_total.dissolved",
                "isotopes": {"Fe2_liq_54": "Fe"},
            },
            "Fe2_sorbed": {
                "phase": "solid",
                "alias_of": "Fe2_total.sorbed",
                "isotopes": {"Fe2_sorbed_54": "Fe"},
            },
        },
    },
}
```

---

## Running the Simulation

Run the simulation from the command line:

```bash
python run_chromatographic_fractionation.py
```

The script outputs:
- `run_chromatographic_fractionation.csv`: Full spatial profiles of concentrations, fluxes, and $\delta^{56}\text{Fe}$ values.
- `run_chromatographic_fractionation_state.npz`: Complete state snapshot for resumption or analysis.
- `run_chromatographic_fractionation.pdf` and `.svg`: Summary visualization plots overlaid against `analytical_ogata_banks.csv`.
